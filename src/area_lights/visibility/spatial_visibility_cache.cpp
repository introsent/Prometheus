//
// Created by ivans on 28/12/2025.
//

#include "spatial_visibility_cache.h"
#include <algorithm>
#include <iostream>
#include <glm/glm.hpp>

SpatialVisibilityCache::SpatialVisibilityCache(
    const glm::vec3& sceneMin,
    const glm::vec3& sceneMax,
    int resolution)
    : m_sceneMin(sceneMin)
    , m_sceneMax(sceneMax) {

    const glm::vec3 sceneExtent = sceneMax - sceneMin;

    // compute cell size based on target resolution
    const float maxExtent = std::max({sceneExtent.x, sceneExtent.y, sceneExtent.z});
    const float cellSizeTarget = maxExtent / static_cast<float>(resolution);

    // compute actual resolution for each axis
    m_resolution.x = std::max(1, static_cast<int>(sceneExtent.x / cellSizeTarget));
    m_resolution.y = std::max(1, static_cast<int>(sceneExtent.y / cellSizeTarget));
    m_resolution.z = std::max(1, static_cast<int>(sceneExtent.z / cellSizeTarget));

    m_cellSize = sceneExtent / glm::vec3(m_resolution);

    // allocate cells
    const int totalCells = m_resolution.x * m_resolution.y * m_resolution.z;
    m_cells.reserve(totalCells);

    // initialize cell centers and radii
    for (int z = 0; z < m_resolution.z; ++z) {
        for (int y = 0; y < m_resolution.y; ++y) {
            for (int x = 0; x < m_resolution.x; ++x) {
                const glm::vec3 cellMin = m_sceneMin + glm::vec3(x, y, z) * m_cellSize;
                const glm::vec3 cellMax = cellMin + m_cellSize;

                VisibilityCacheCell cell;
                cell.center = (cellMin + cellMax) * 0.5f;
                cell.radius = glm::length(m_cellSize) * 0.5f;

                m_cells.push_back(std::move(cell));
            }
        }
    }

    std::cout << "Visibility cache initialized: "
              << m_resolution.x << "x" << m_resolution.y << "x" << m_resolution.z
              << " = " << totalCells << " cells\n";
}

glm::ivec3 SpatialVisibilityCache::pointToCellCoord(const glm::vec3& point) const {
    const glm::vec3 normalized = (point - m_sceneMin) / (m_sceneMax - m_sceneMin);
    const glm::vec3 coord = normalized * glm::vec3(m_resolution);

    return {
        glm::clamp(static_cast<int>(coord.x), 0, m_resolution.x - 1),
        glm::clamp(static_cast<int>(coord.y), 0, m_resolution.y - 1),
        glm::clamp(static_cast<int>(coord.z), 0, m_resolution.z - 1)
    };
}

int SpatialVisibilityCache::coordToCellIndex(const glm::ivec3& coord) const {
    return coord.x +
           coord.y * m_resolution.x +
           coord.z * m_resolution.x * m_resolution.y;
}

bool SpatialVisibilityCache::isValidCoord(const glm::ivec3& coord) const {
    return coord.x >= 0 && coord.x < m_resolution.x &&
           coord.y >= 0 && coord.y < m_resolution.y &&
           coord.z >= 0 && coord.z < m_resolution.z;
}

int SpatialVisibilityCache::getCellIndex(const glm::vec3& point) const {
    const glm::ivec3 coord = pointToCellCoord(point);
    return coordToCellIndex(coord);
}

VisibilityCacheCell* SpatialVisibilityCache::getCell(const glm::vec3& point) {
    const int idx = getCellIndex(point);
    if (idx >= 0 && idx < static_cast<int>(m_cells.size())) {
        return &m_cells[idx];
    }
    return nullptr;
}

void SpatialVisibilityCache::recordVisibility(
    const glm::vec3& shadingPoint,
    int nodeIndex,
    bool wasVisible) {

    if (auto* cell = getCell(shadingPoint)) {
        auto* stats = cell->getOrCreateNodeStats(nodeIndex);
        stats->recordSample(wasVisible);
    }
}

float SpatialVisibilityCache::queryVisibility(
    const glm::vec3& shadingPoint,
    int nodeIndex) const {
    if (const int idx = getCellIndex(shadingPoint); idx >= 0 && idx < static_cast<int>(m_cells.size())) {
        for (const auto& cell = m_cells[idx]; const auto& stats : cell.nodeStats) {
            if (stats.nodeIndex == nodeIndex) {
                return stats.getVisibilityProbability();
            }
        }
    }

    // default: optimistically assume visible
    return 1.0f;
}

void SpatialVisibilityCache::getStatistics(
    int& totalCells,
    int& activeCells,
    int& totalSamples) const {

    totalCells = static_cast<int>(m_cells.size());
    activeCells = 0;
    totalSamples = 0;

    for (const auto& cell : m_cells) {
        if (!cell.nodeStats.empty()) {
            activeCells++;
            for (const auto& stats : cell.nodeStats) {
                totalSamples += stats.totalSamples.load(std::memory_order_relaxed);
            }
        }
    }
}

void SpatialVisibilityCache::reset() {
    for (auto& cell : m_cells) {
        std::lock_guard<std::mutex> lock(cell.cellMutex);
        cell.nodeStats.clear();
    }
}