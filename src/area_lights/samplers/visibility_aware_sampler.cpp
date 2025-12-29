//
// Created by ivans on 28/12/2025.
//

#include "visibility_aware_sampler.h"

#include <iostream>
#include <ostream>

#include "area_importance_triangle_sampler.h"
#include "mis/bsdf_sampler.h"
#include "visibility/spatial_visibility_cache.h"

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
VisibilityAwareHierarchicalSampler::sampleLight(const glm::vec3& shadingPoint,
                                               const glm::vec3& normal,
                                               float u1, float u2, float u3) const {
    Sample result{};
    result.isValid = false;
    result.misWeight = 1.0f;

    // get root node
    if (m_light->m_rootNodeIndex < 0 ||
        m_light->m_bvhNodes.empty()) {
        return result;
    }

    // 1. select leaf node using visibility-aware hierarchical sampling
    float pathPdf = 1.0f;
    int leafNodeIndex = selectBVHNodeWithVisibility(
        m_light->m_rootNodeIndex,
        shadingPoint,
        u1,
        pathPdf
    );
    if (leafNodeIndex < 0) {
        return result;
    }

    const auto& leaf = m_light->m_bvhNodes[leafNodeIndex];
    if (!leaf.isLeaf || leaf.triangleIndices.empty()) {
        return result;
    }

    m_light->setLastSampledNode(leafNodeIndex);

    // 2. select triangle within leaf (weighted by flux)
    float fluxSum = 0.0f;
    std::vector<float> triangleFluxes;

    for (int triIdx : leaf.triangleIndices) {
        const auto& tri = m_light->m_triangles[triIdx];
        float flux = tri.intensity * tri.area *
                    ((m_light->m_emission.r + m_light->m_emission.g + m_light->m_emission.b) / 3.0f);
        triangleFluxes.push_back(flux);
        fluxSum += flux;
    }

    if (fluxSum <= 0.0f) {
        return result;
    }

    // select triangle
    float target = u2 * fluxSum;
    float cumulative = 0.0f;
    int selectedIdx = 0;

    for (size_t i = 0; i < leaf.triangleIndices.size(); ++i) {
        cumulative += triangleFluxes[i];
        if (target <= cumulative) {
            selectedIdx = static_cast<int>(i);
            break;
        }
    }

    int triangleIndex = leaf.triangleIndices[selectedIdx];
    const auto& triangle = m_light->m_triangles[triangleIndex];

    // 3. sample point on triangle
    auto triangleSample = triangle.areaImportanceSampler->sample(
        shadingPoint, u3, std::fmod(u1 + u2, 1.0f)
    );

    if (triangleSample.pdf <= 0.0f) {
        return result;
    }

    // 4. calculate combined PDF
    float triangleSelectionProb = triangleFluxes[selectedIdx] / fluxSum;
    float lightSamplingPdf = pathPdf * triangleSelectionProb * triangleSample.pdf;

    // 5. calculate MIS weight if enabled
    float misWeight = 1.0f;
    if (m_config.enableMIS) {
        glm::vec3 lightDir = glm::normalize(triangleSample.position - shadingPoint);
        float bsdfPdf = BSDFSampler::pdfDiffuse(normal, lightDir);

        float distanceSq = glm::dot(triangleSample.position - shadingPoint,
                                  triangleSample.position - shadingPoint);
        float cosLight = glm::dot(triangleSample.normal, -lightDir);
        cosLight = std::max(cosLight, 1e-4f);

        float lightPdfSolidAngle = lightSamplingPdf * cosLight / distanceSq;

        misWeight = MISWeightCalculator::calculateWeight(
            lightPdfSolidAngle,
            bsdfPdf,
            m_config.misHeuristic
        );
    }

    // Fill result
    result.position = triangleSample.position;
    result.normal = triangleSample.normal;
    result.radiance = triangleSample.radiance;
    result.pdf = lightSamplingPdf;
    result.misWeight = misWeight;
    result.isValid = true;

    return result;
}

int VisibilityAwareHierarchicalSampler::selectBVHNodeWithVisibility(
    int nodeIndex,
    const glm::vec3& shadingPoint,
    float u,
    float& outPdf) const {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return -1;
    }

    const auto& node = m_light->m_bvhNodes[nodeIndex];

    // base case: leaf node
    if (node.isLeaf) {
        return nodeIndex;
    }

    // calculate importance for each child
    const float importanceLeft = calculateNodeImportance(node.leftChild, shadingPoint);
    const float importanceRight = calculateNodeImportance(node.rightChild, shadingPoint);

    const float totalImportance = importanceLeft + importanceRight;
    if (totalImportance <= 0.0f) {
        return -1;
    }

    // calculate probabilities

    // traverse and update PDF
    if (float probLeft = importanceLeft / totalImportance; u < probLeft) {
        outPdf *= probLeft;
        float newU = (probLeft > 0.0f) ? u / probLeft : 0.0f;
        return selectBVHNodeWithVisibility(node.leftChild, shadingPoint, newU, outPdf);
    } else {
        float probRight = 1.0f - probLeft;
        outPdf *= probRight;
        float newU = (probRight > 0.0f) ? (u - probLeft) / probRight : 0.0f;
        return selectBVHNodeWithVisibility(node.rightChild, shadingPoint, newU, outPdf);
    }
}

float VisibilityAwareHierarchicalSampler::calculateNodeImportance(
    int nodeIndex,
    const glm::vec3& shadingPoint) const {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(m_light->m_bvhNodes.size())) {
        return 0.0f;
    }

    const auto& node = m_light->m_bvhNodes[nodeIndex];

    // 1. calculate geometric importance (flux / distance^2)
    const glm::vec3 toNode = node.getCentroid() - shadingPoint;
    float distSq = glm::dot(toNode, toNode);
    distSq = std::max(distSq, 1e-4f); // Avoid division by zero

    float geometricImportance = node.totalFlux / distSq;

    // 2. query visibility probability if enabled
    float visibilityProb = 1.0f;
    if (m_config.enableVisibilityLearning && m_visCache) {
        visibilityProb = m_visCache->queryVisibility(shadingPoint, nodeIndex);

        // only use visibility if we have enough samples
        if (auto* cell = m_visCache->getCell(shadingPoint)) {
            if (const auto* nodeVis = cell->getNodeStats(nodeIndex)) {
                if (const uint32_t samples = nodeVis->totalSamples.load(std::memory_order_relaxed); samples < static_cast<uint32_t>(m_config.minSamplesForLearning)) {
                    // blend toward optimistic estimate
                    const float confidence = static_cast<float>(samples) /
                                     static_cast<float>(m_config.minSamplesForLearning);
                    visibilityProb = visibilityProb * confidence + 1.0f * (1.0f - confidence);
                }
            }
        }
    }

    // 3. combine geometric and visibility importance
    float importance = geometricImportance *
                      glm::mix(1.0f, visibilityProb, m_config.visibilityWeight);

    return importance;
}

float VisibilityAwareHierarchicalSampler::evaluatePdf(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    // find which triangle contains the point
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

    // traverse BVH to compute hierarchical PDF
    return m_light->computeHierarchicalPdfForTriangle(
        triangleIndex, shadingPoint, lightPoint
    );
}

void VisibilityAwareHierarchicalSampler::recordVisibilitySample(
    const glm::vec3& shadingPoint,
    int nodeIndex,
    bool wasVisible) {

    if (m_config.enableVisibilityLearning && m_visCache) {
        m_visCache->recordVisibility(shadingPoint, nodeIndex, wasVisible);
    }
}
