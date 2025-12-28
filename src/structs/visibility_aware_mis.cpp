//
// Created by ivans on 28/12/2025.
//

#include "visibility_aware_mis.h"

#include "area_light.h"
#include <algorithm>
#include <cmath>
#include <iostream>

/// Spatial Visibility Cache Implementation
SpatialVisibilityCache::SpatialVisibilityCache(const glm::vec3& sceneMin,
                                               const glm::vec3& sceneMax,
                                               int resolution)
    : m_sceneMin(sceneMin)
    , m_sceneMax(sceneMax)
{
    glm::vec3 sceneExtent = sceneMax - sceneMin;

    float maxExtent = std::max({sceneExtent.x, sceneExtent.y, sceneExtent.z});
    float cellSizeTarget = maxExtent / static_cast<float>(resolution);

    m_resolution.x = std::max(1, static_cast<int>(sceneExtent.x / cellSizeTarget));
    m_resolution.y = std::max(1, static_cast<int>(sceneExtent.y / cellSizeTarget));
    m_resolution.z = std::max(1, static_cast<int>(sceneExtent.z / cellSizeTarget));

    m_cellSize = sceneExtent / glm::vec3(m_resolution);

    // allocate cells
    int totalCells = m_resolution.x * m_resolution.y * m_resolution.z;
    m_cells.resize(totalCells);
    // construct cells in place
    for (int z = 0; z < m_resolution.z; ++z) {
        for (int y = 0; y < m_resolution.y; ++y) {
            for (int x = 0; x < m_resolution.x; ++x) {
                int index = x + y * m_resolution.x + z * m_resolution.x * m_resolution.y;
                m_cells[index] = VisibilityCacheCell();

                glm::ivec3 coord(x, y, z);

                glm::vec3 cellMin = m_sceneMin + glm::vec3(coord) * m_cellSize;
                glm::vec3 cellMax = cellMin + m_cellSize;

                // construct cell in place
                VisibilityCacheCell cell;
                cell.center = (cellMin + cellMax) * 0.5f;
                cell.radius = glm::length(m_cellSize) * 0.5f;

                m_cells.push_back(std::move(cell));
            }
        }
    }

    std::cout << "Visibility cache initialized: "
              << m_resolution.x << "x" << m_resolution.y << "x" << m_resolution.z
              << " = " << totalCells << " cells" << std::endl;
}

glm::ivec3 SpatialVisibilityCache::pointToCellCoord(const glm::vec3& point) const {
    glm::vec3 normalized = (point - m_sceneMin) / (m_sceneMax - m_sceneMin);
    glm::vec3 coord = normalized * glm::vec3(m_resolution);

    return {
        glm::clamp(static_cast<int>(coord.x), 0, m_resolution.x - 1),
        glm::clamp(static_cast<int>(coord.y), 0, m_resolution.y - 1),
        glm::clamp(static_cast<int>(coord.z), 0, m_resolution.z - 1)
    };
}

int SpatialVisibilityCache::coordToCellIndex(const glm::ivec3& coord) const {
    return coord.x + coord.y * m_resolution.x + coord.z * m_resolution.x * m_resolution.y;
}

bool SpatialVisibilityCache::isValidCoord(const glm::ivec3& coord) const {
    return coord.x >= 0 && coord.x < m_resolution.x &&
           coord.y >= 0 && coord.y < m_resolution.y &&
           coord.z >= 0 && coord.z < m_resolution.z;
}

int SpatialVisibilityCache::getCellIndex(const glm::vec3& point) const {
    glm::ivec3 coord = pointToCellCoord(point);
    return coordToCellIndex(coord);
}

VisibilityCacheCell* SpatialVisibilityCache::getCell(const glm::vec3& point) {
    if (int idx = getCellIndex(point); idx >= 0 && idx < static_cast<int>(m_cells.size())) {
        return &m_cells[idx];
    }
    return nullptr;
}

void SpatialVisibilityCache::recordVisibility(const glm::vec3& shadingPoint,
                                              int nodeIndex,
                                              bool wasVisible) {
    if (auto* cell = getCell(shadingPoint)) {
        std::lock_guard<std::mutex> lock(cell->cellMutex);
        auto* nodeVis = cell->getOrCreateNodeVisibility(nodeIndex);
        nodeVis->recordSample(wasVisible);
    }
}

float SpatialVisibilityCache::queryVisibility(const glm::vec3& shadingPoint,
                                              int nodeIndex) const {
    if (int idx = getCellIndex(shadingPoint); idx >= 0 && idx < static_cast<int>(m_cells.size())) {
        for (const auto& cell = m_cells[idx]; const auto& nv : cell.nodeVisibilities) {
            if (nv.nodeIndex == nodeIndex) {
                return nv.getVisibilityProbability();
            }
        }
    }

    // default: optimistically assume visible
    return 1.0f;
}

void SpatialVisibilityCache::getStatistics(int& totalCells, int& activeCells,
                                          int& totalSamples) const {
    totalCells = static_cast<int>(m_cells.size());
    activeCells = 0;
    totalSamples = 0;

    for (const auto& cell : m_cells) {
        if (!cell.nodeVisibilities.empty()) {
            activeCells++;
            for (const auto& nv : cell.nodeVisibilities) {
                totalSamples += nv.totalSamples.load(std::memory_order_relaxed);
            }
        }
    }
}

void SpatialVisibilityCache::reset() {
    for (auto& cell : m_cells) {
        std::lock_guard<std::mutex> lock(cell.cellMutex);
        cell.nodeVisibilities.clear();
    }
}

/// MIS Weight Calculator Implementation
float MISWeightCalculator::calculateWeight(float thisPdf, float otherPdf,
                                          Heuristic heuristic) {
    if (thisPdf <= 0.0f) return 0.0f;
    if (otherPdf <= 0.0f) return 1.0f;

    switch (heuristic) {
        case Heuristic::Balance:
            return balanceHeuristic(thisPdf, otherPdf);
        case Heuristic::Power:
            return powerHeuristic(thisPdf, otherPdf);
        case Heuristic::Maximum:
            return thisPdf > otherPdf ? 1.0f : 0.0f;
        default:
            return balanceHeuristic(thisPdf, otherPdf);
    }
}

float MISWeightCalculator::balanceHeuristic(float thisPdf, float otherPdf) {
    return thisPdf / (thisPdf + otherPdf);
}

float MISWeightCalculator::powerHeuristic(float thisPdf, float otherPdf, float beta) {
    float thisWeight = std::pow(thisPdf, beta);
    float otherWeight = std::pow(otherPdf, beta);
    return thisWeight / (thisWeight + otherWeight);
}


/// Visibility-Aware Hierarchical Sampler Implementation
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
        // calculate BSDF sampling PDF for this direction
        glm::vec3 lightDir = glm::normalize(triangleSample.position - shadingPoint);
        float bsdfPdf = BSDFSampler::pdfDiffuse(normal, lightDir);

        misWeight = MISWeightCalculator::calculateWeight(
            lightSamplingPdf,
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

    // get children
    const auto& left = m_light->m_bvhNodes[node.leftChild];
    const auto& right = m_light->m_bvhNodes[node.rightChild];

    // calculate importance for each child
    float importanceLeft = calculateNodeImportance(node.leftChild, shadingPoint);
    float importanceRight = calculateNodeImportance(node.rightChild, shadingPoint);

    float totalImportance = importanceLeft + importanceRight;
    if (totalImportance <= 0.0f) {
        return -1;
    }

    // calculate probabilities
    float probLeft = importanceLeft / totalImportance;

    // traverse and update PDF
    if (u < probLeft) {
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
            if (auto* nodeVis = cell->getNodeVisibility(nodeIndex)) {
                if (uint32_t samples = nodeVis->totalSamples.load(std::memory_order_relaxed); samples < static_cast<uint32_t>(m_config.minSamplesForLearning)) {
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

    // Find which triangle contains the point
    size_t triangleIndex = m_light->m_triangles.size();
    for (size_t i = 0; i < m_light->m_triangles.size(); ++i) {
        const auto& tri = m_light->m_triangles[i];
        float u, v, w;
        if (m_light->isPointInTriangle(lightPoint, tri.v0, tri.v1, tri.v2, u, v, w)) {
            triangleIndex = i;
            break;
        }
    }

    if (triangleIndex >= m_light->m_triangles.size()) {
        return 0.0f;
    }

    // Traverse BVH to compute hierarchical PDF
    // This mirrors the sampling procedure
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

// ============================================================================
// BSDF Sampler Implementation
// ============================================================================

BSDFSampler::BSDFSample BSDFSampler::sampleDiffuse(const glm::vec3& normal,
                                                   float u1, float u2) {
    BSDFSample sample{};

    // cosine-weighted hemisphere sampling
    glm::vec3 localDir = VisibilityAwareHelpers::cosineSampleHemisphere(u1, u2);
    sample.direction = VisibilityAwareHelpers::alignHemisphereWithNormal(localDir, normal);

    // PDF for cosine-weighted sampling: cos(theta) / pi
    float cosTheta = std::max(0.0f, glm::dot(sample.direction, normal));
    sample.pdf = cosTheta / glm::pi<float>();

    // diffuse BRDF is 1/pi (albedo assumed to be 1)
    sample.reflectance = glm::vec3(1.0f / glm::pi<float>());

    return sample;
}

float BSDFSampler::pdfDiffuse(const glm::vec3& normal, const glm::vec3& direction) {
    float cosTheta = std::max(0.0f, glm::dot(direction, normal));
    return cosTheta / glm::pi<float>();
}