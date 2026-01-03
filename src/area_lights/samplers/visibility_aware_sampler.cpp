//
// Created by ivans on 28/12/2025.
//

#include "visibility_aware_sampler.h"
#include "lights/mesh_area_light.h"
#include "visibility/spatial_visibility_cache.h"
#include "samplers/area_importance_triangle_sampler.h"
#include "mis/bsdf_sampler.h"
#include <cmath>
#include <algorithm>

VisibilityAwareHierarchicalSampler::VisibilityAwareHierarchicalSampler(
    MeshAreaLight* light,
    SpatialVisibilityCache* visCache,
    const Configuration& config)
    : m_light(light)
    , m_visCache(visCache)
    , m_config(config)
{
}

VisibilityAwareHierarchicalSampler::Sample
VisibilityAwareHierarchicalSampler::sampleLight(
    const glm::vec3& shadingPoint,
    const glm::vec3& shadingNormal,
    float u1, float u2, float u3) const {

    Sample result{};
    result.isValid = false;

    if (m_light->m_rootNodeIndex < 0 || m_light->m_bvhNodes.empty()) {
        return result;
    }

    // Use local traversal path (thread-safe)
    std::vector<int> traversalPath;
    traversalPath.reserve(32);  // Pre-allocate to reduce allocations

    // Traverse with adaptive splitting
    TraversalResult traversal = traverseWithAdaptiveSplitting(
        m_light->m_rootNodeIndex,
        shadingPoint,
        shadingNormal,
        u1,
        0,
        traversalPath
    );

    if (traversal.leafNodeIndex < 0 && traversal.sampledLeaves.empty()) {
        return result;
    }

    // Handle single leaf case (no splitting occurred)
    int leafNodeIndex;
    float pathPdf;

    if (!traversal.sampledLeaves.empty()) {
        // Multiple leaves from splitting
        float totalWeight = 0.0f;
        for (const auto& [idx, weight] : traversal.sampledLeaves) {
            totalWeight += weight;
        }

        if (totalWeight <= 0.0f) {
            return result;
        }

        float target = u2 * totalWeight;
        float cumulative = 0.0f;
        leafNodeIndex = traversal.sampledLeaves[0].first;

        float selectedWeight = traversal.sampledLeaves[0].second;

        for (const auto& [idx, weight] : traversal.sampledLeaves) {
            cumulative += weight;
            if (target <= cumulative) {
                leafNodeIndex = idx;
                selectedWeight = weight;
                break;
            }
        }

        // PDF for selecting this specific leaf from the split
        pathPdf = selectedWeight / totalWeight;

    } else {
        leafNodeIndex = traversal.leafNodeIndex;
        pathPdf = traversal.pathPdf;
    }

    if (leafNodeIndex < 0 || leafNodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return result;
    }

    const auto& leaf = m_light->m_bvhNodes[leafNodeIndex];
    if (!leaf.isLeaf || leaf.triangleIndices.empty()) {
        return result;
    }

    m_light->setLastSampledNode(leafNodeIndex);

    // Select triangle within leaf (weighted by flux)
    float fluxSum = 0.0f;
    std::vector<float> triangleFluxes;
    triangleFluxes.reserve(leaf.triangleIndices.size());

    for (int triIdx : leaf.triangleIndices) {
        if (triIdx < 0 || triIdx >= static_cast<int>(m_light->m_triangles.size())) {
            continue;
        }
        const auto& tri = m_light->m_triangles[triIdx];
        float flux = tri.intensity * tri.area *
                    ((m_light->m_emission.r + m_light->m_emission.g + m_light->m_emission.b) / 3.0f);
        triangleFluxes.push_back(flux);
        fluxSum += flux;
    }

    if (fluxSum <= 0.0f || triangleFluxes.empty()) {
        return result;
    }

    // Select triangle
    float adjustedU2 = traversal.sampledLeaves.empty() ? u2 : std::fmod(u2 * 7.13f, 1.0f);
    float target = adjustedU2 * fluxSum;
    float cumulative = 0.0f;
    size_t selectedIdx = 0;

    for (size_t i = 0; i < triangleFluxes.size(); ++i) {
        cumulative += triangleFluxes[i];
        if (target <= cumulative) {
            selectedIdx = i;
            break;
        }
    }

    if (selectedIdx >= leaf.triangleIndices.size()) {
        return result;
    }

    int triangleIndex = leaf.triangleIndices[selectedIdx];
    if (triangleIndex < 0 || triangleIndex >= static_cast<int>(m_light->m_triangles.size())) {
        return result;
    }

    const auto& triangle = m_light->m_triangles[triangleIndex];

    // Sample point on triangle
    auto triangleSample = triangle.areaImportanceSampler->sample(
        shadingPoint, u3, std::fmod(u1 + u2 + u3, 1.0f)
    );

    if (triangleSample.pdf <= 0.0f) {
        return result;
    }

    // Calculate combined PDF
    float triangleSelectionProb = triangleFluxes[selectedIdx] / fluxSum;
    float lightSamplingPdf = pathPdf * triangleSelectionProb * triangleSample.pdf;

    if (lightSamplingPdf <= 0.0f) {
        return result;
    }

    // Calculate MIS weight
    float misWeight = 1.0f;
    if (m_config.enableMIS) {
        glm::vec3 toLight = triangleSample.position - shadingPoint;
        float distanceSq = glm::dot(toLight, toLight);

        if (distanceSq > 1e-8f) {
            glm::vec3 lightDir = toLight / std::sqrt(distanceSq);
            float bsdfPdf = BSDFSampler::pdfDiffuse(shadingNormal, lightDir);

            float cosLight = std::max(glm::dot(triangleSample.normal, -lightDir), 1e-4f);
            float lightPdfSolidAngle = lightSamplingPdf * distanceSq / cosLight;

            misWeight = MISWeightCalculator::calculateWeight(
                lightPdfSolidAngle,
                bsdfPdf,
                m_config.misHeuristic
            );
        }
    }

    // Fill result
    result.position = triangleSample.position;
    result.normal = triangleSample.normal;
    result.radiance = triangleSample.radiance;
    result.pdf = lightSamplingPdf;
    result.misWeight = misWeight;
    result.isValid = true;

    // Store traversal path in result for later visibility recording
    result.traversalPath = std::move(traversalPath);

    return result;
}

float VisibilityAwareHierarchicalSampler::calculateNodeImportance(
    int nodeIndex,
    const glm::vec3& shadingPoint,
    const glm::vec3& shadingNormal) const {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return 0.0f;
    }

    const auto& node = m_light->m_bvhNodes[nodeIndex];

    // Vector from shading point to cluster center
    glm::vec3 toNode = node.getCentroid() - shadingPoint;
    float dist = glm::length(toNode);
    float distSq = std::max(dist * dist, 1e-4f);

    // Clamp distance to half bounding sphere radius when inside/near cluster
    float boundingRadius = node.getBoundingSphereRadius();
    if (dist < boundingRadius) {
        dist = std::max(dist, boundingRadius * 0.5f);
        distSq = dist * dist;
    }

    glm::vec3 dirToNode = (dist > 1e-4f) ? toNode / dist : glm::vec3(0.0f, 1.0f, 0.0f);

    // theta_u: angle of cone covering bounding box from shading point
    float thetaU = (dist > 1e-4f) ? std::atan(boundingRadius / dist) : static_cast<float>(M_PI);
    thetaU = std::min(thetaU, static_cast<float>(M_PI));

    // theta_i: incident angle from shading point normal to cluster center
    float cosThetaI = glm::dot(shadingNormal, dirToNode);
    float thetaI = std::acos(glm::clamp(cosThetaI, -1.0f, 1.0f));

    // theta_i' = max(theta_i - theta_u, 0)
    float thetaIPrime = std::max(thetaI - thetaU, 0.0f);

    // theta: angle between cluster orientation axis and direction to shading point
    float cosTheta = glm::dot(node.orientationBounds.axis, -dirToNode);
    float theta = std::acos(glm::clamp(cosTheta, -1.0f, 1.0f));

    // theta' = max(theta - theta_o - theta_u, 0)
    float thetaPrime = std::max(
        theta - node.orientationBounds.thetaO - thetaU,
        0.0f
    );

    // Check if emission is within profile
    float orientationFactor = 0.0f;
    if (thetaPrime < node.orientationBounds.thetaE) {
        orientationFactor = std::cos(thetaPrime);
    }

    // Incident angle factor
    float incidentFactor = std::abs(std::cos(thetaIPrime));

    // Base geometric importance (Equation 3)
    float geometricImportance = (incidentFactor * orientationFactor * node.totalFlux) / distSq;

    // Query visibility probability
    float visibilityProb = 1.0f;
    if (m_config.enableVisibilityLearning && m_visCache) {
        visibilityProb = m_visCache->queryVisibility(shadingPoint, nodeIndex);

        // Confidence-based blending
        if (auto* cell = m_visCache->getCell(shadingPoint)) {
            if (const auto* nodeVis = cell->getNodeStats(nodeIndex)) {
                uint32_t samples = nodeVis->totalSamples.load(std::memory_order_relaxed);

                if (samples < static_cast<uint32_t>(m_config.minSamplesForLearning)) {
                    float observedProb = nodeVis->getVisibilityProbability();
                    float uncertainty = 1.0f / std::sqrt(1.0f + static_cast<float>(samples));
                    visibilityProb = std::min(observedProb + uncertainty, 1.0f);
                }
            }
        }
    }

    // Combine geometric and visibility importance
    float importance = geometricImportance *
                      glm::mix(1.0f, visibilityProb, m_config.visibilityWeight);

    return std::max(importance, 0.0f);
}

bool VisibilityAwareHierarchicalSampler::shouldSplit(
    int nodeIndex,
    const glm::vec3& shadingPoint,
    int currentDepth) const {

    if (!m_config.enableAdaptiveSplitting) {
        return false;
    }

    if (currentDepth >= m_config.maxSplitDepth) {
        return false;
    }

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return false;
    }

    const auto& node = m_light->m_bvhNodes[nodeIndex];
    if (node.isLeaf) {
        return false;
    }

    float varianceScore = computeSplitVarianceScore(nodeIndex, shadingPoint);

    // Remap to [0,1] using fourth root as in paper
    float normalizedScore = std::pow(1.0f / (1.0f + varianceScore), 0.25f);

    // slit if variance is high (normalized score is low)
    return normalizedScore < (1.0f - m_config.splitThreshold);
}

float VisibilityAwareHierarchicalSampler::computeSplitVarianceScore(
    int nodeIndex,
    const glm::vec3& shadingPoint) const {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return 0.0f;
    }

    const auto& node = m_light->m_bvhNodes[nodeIndex];

    // distance range to cluster (a, b)
    glm::vec3 toCenter = node.getCentroid() - shadingPoint;
    float centerDist = glm::length(toCenter);
    float radius = node.getBoundingSphereRadius();

    float a = std::max(centerDist - radius, 1e-4f);
    float b = centerDist + radius;

    // variance of geometric term 1/d^2
    float E_g = 1.0f / (a * b);

    float range = b - a;
    float E_g2 = (range > 1e-6f)
        ? (b*b*b - a*a*a) / (3.0f * range * a*a*a * b*b*b)
        : E_g * E_g;

    float V_g = std::max(E_g2 - E_g * E_g, 0.0f);

    // flux variance
    float V_e = node.fluxVariance;
    int numTris = node.isLeaf ? static_cast<int>(node.triangleIndices.size()) : 1;
    float E_e = node.totalFlux / std::max(static_cast<float>(numTris), 1.0f);

    // total variance (Equation 10)
    int N = std::max(numTris, 1);
    float variance = (V_e * V_g + V_e * E_g * E_g + E_e * E_e * V_g) * static_cast<float>(N * N);

    return variance;
}

VisibilityAwareHierarchicalSampler::TraversalResult
VisibilityAwareHierarchicalSampler::traverseWithAdaptiveSplitting(
    int nodeIndex,
    const glm::vec3& shadingPoint,
    const glm::vec3& shadingNormal,
    float u,
    int depth,
    std::vector<int>& outTraversalPath) const {

    TraversalResult result;

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return result;
    }

    // record this node in traversal path
    outTraversalPath.push_back(nodeIndex);

    const auto& node = m_light->m_bvhNodes[nodeIndex];

    // base case: leaf node
    if (node.isLeaf) {
        result.leafNodeIndex = nodeIndex;
        result.pathPdf = 1.0f;
        return result;
    }

    // validate children
    if (node.leftChild < 0 || node.rightChild < 0 ||
        node.leftChild >= static_cast<int>(m_light->m_bvhNodes.size()) ||
        node.rightChild >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return result;
    }

    // check if we should split
    if (shouldSplit(nodeIndex, shadingPoint, depth)) {
        // sample BOTH children
        std::vector<int> leftPath, rightPath;
        leftPath.reserve(16);
        rightPath.reserve(16);

        auto leftResult = traverseWithAdaptiveSplitting(
            node.leftChild, shadingPoint, shadingNormal, u, depth + 1, leftPath);
        auto rightResult = traverseWithAdaptiveSplitting(
            node.rightChild, shadingPoint, shadingNormal,
            std::fmod(u + 0.5f, 1.0f), depth + 1, rightPath);

        // merge paths into main path
        outTraversalPath.insert(outTraversalPath.end(), leftPath.begin(), leftPath.end());
        outTraversalPath.insert(outTraversalPath.end(), rightPath.begin(), rightPath.end());

        // collect all leaves from both branches
        float leftImportance = calculateNodeImportance(node.leftChild, shadingPoint, shadingNormal);
        float rightImportance = calculateNodeImportance(node.rightChild, shadingPoint, shadingNormal);
        float totalImportance = leftImportance + rightImportance;

        if (totalImportance <= 0.0f) {
            return result;
        }

        if (leftResult.leafNodeIndex >= 0) {
            result.sampledLeaves.emplace_back(
                leftResult.leafNodeIndex,
                leftImportance
            );
        }
        for (const auto& leaf : leftResult.sampledLeaves) {
            result.sampledLeaves.emplace_back(leaf.first, leftImportance * leaf.second);
        }

        if (rightResult.leafNodeIndex >= 0) {
            result.sampledLeaves.emplace_back(
                rightResult.leafNodeIndex,
                rightImportance
            );
        }
        for (const auto& leaf : rightResult.sampledLeaves) {
            result.sampledLeaves.emplace_back(leaf.first, rightImportance * leaf.second);
        }

        return result;
    }

    // Stochastic traversal (no split)
    float importanceLeft = calculateNodeImportance(node.leftChild, shadingPoint, shadingNormal);
    float importanceRight = calculateNodeImportance(node.rightChild, shadingPoint, shadingNormal);

    float totalImportance = importanceLeft + importanceRight;
    if (totalImportance <= 0.0f) {
        return result;
    }

    float probLeft = importanceLeft / totalImportance;

    if (u < probLeft) {
        result = traverseWithAdaptiveSplitting(
            node.leftChild, shadingPoint, shadingNormal,
            (probLeft > 0.0f) ? u / probLeft : 0.0f,
            depth + 1, outTraversalPath
        );
        result.pathPdf *= probLeft;
    } else {
        float probRight = 1.0f - probLeft;
        result = traverseWithAdaptiveSplitting(
            node.rightChild, shadingPoint, shadingNormal,
            (probRight > 0.0f) ? (u - probLeft) / probRight : 0.0f,
            depth + 1, outTraversalPath
        );
        result.pathPdf *= probRight;
    }

    return result;
}

int VisibilityAwareHierarchicalSampler::selectBVHNodeStochastic(
    int nodeIndex,
    const glm::vec3& shadingPoint,
    const glm::vec3& shadingNormal,
    float u,
    float& outPdf,
    std::vector<int>& outTraversalPath) const {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return -1;
    }

    outTraversalPath.push_back(nodeIndex);

    const auto& node = m_light->m_bvhNodes[nodeIndex];

    if (node.isLeaf) {
        return nodeIndex;
    }

    if (node.leftChild < 0 || node.rightChild < 0) {
        return -1;
    }

    float importanceLeft = calculateNodeImportance(node.leftChild, shadingPoint, shadingNormal);
    float importanceRight = calculateNodeImportance(node.rightChild, shadingPoint, shadingNormal);

    float totalImportance = importanceLeft + importanceRight;
    if (totalImportance <= 0.0f) {
        return -1;
    }

    float probLeft = importanceLeft / totalImportance;

    if (u < probLeft) {
        outPdf *= probLeft;
        float newU = (probLeft > 0.0f) ? u / probLeft : 0.0f;
        return selectBVHNodeStochastic(
            node.leftChild, shadingPoint, shadingNormal, newU, outPdf, outTraversalPath);
    } else {
        float probRight = 1.0f - probLeft;
        outPdf *= probRight;
        float newU = (probRight > 0.0f) ? (u - probLeft) / probRight : 0.0f;
        return selectBVHNodeStochastic(
            node.rightChild, shadingPoint, shadingNormal, newU, outPdf, outTraversalPath);
    }
}

void VisibilityAwareHierarchicalSampler::recordVisibilitySample(
    const glm::vec3& shadingPoint,
    int nodeIndex,
    bool wasVisible) {

    if (m_config.enableVisibilityLearning && m_visCache) {
        m_visCache->recordVisibility(shadingPoint, nodeIndex, wasVisible);
    }
}

void VisibilityAwareHierarchicalSampler::recordTraversalVisibility(
    const glm::vec3& shadingPoint,
    const std::vector<int>& traversalPath,
    bool wasVisible) {

    if (!m_config.enableVisibilityLearning || !m_visCache) {
        return;
    }

    // record visibility for ALL nodes in the traversal path
    for (int nodeIndex : traversalPath) {
        m_visCache->recordVisibility(shadingPoint, nodeIndex, wasVisible);
    }
}

float VisibilityAwareHierarchicalSampler::evaluatePdf(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    // Find which triangle contains the point
    size_t triangleIndex = m_light->m_triangles.size();
    for (size_t i = 0; i < m_light->m_triangles.size(); ++i) {
        const auto& tri = m_light->m_triangles[i];
        float u, v, w;
        if (MeshAreaLight::isPointInTriangle(lightPoint, tri.v0, tri.v1, tri.v2, u, v, w)) {
            triangleIndex = i;
            break;
        }
    }

    if (triangleIndex >= m_light->m_triangles.size()) {
        return 0.0f;
    }

    return m_light->computeHierarchicalPdfForTriangle(
        triangleIndex, shadingPoint, lightPoint
    );
}