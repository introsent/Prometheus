//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_SPATIAL_VISIBILITY_CACHE_H
#define PROMETHEUS_SPATIAL_VISIBILITY_CACHE_H
#include <glm/vec3.hpp>

#include "visibility_cache_cell.h"

/// Spatial visibility cache
// divides scene into uniform grid cells, each tracking visibility
class SpatialVisibilityCache {
public:
    SpatialVisibilityCache(
        const glm::vec3& sceneMin,
        const glm::vec3& sceneMax,
        int resolution = 16);

    // query visibility probability for a BVH node from a point
    // returns 1.0 if no data available (optimistic)
    [[nodiscard]] float queryVisibility(
        const glm::vec3& shadingPoint,
        int nodeIndex) const;

    // record a visibility sample (thread-safe)
    void recordVisibility(
        const glm::vec3& shadingPoint,
        int nodeIndex,
        bool wasVisible);

    // get cache statistics
    void getStatistics(
        int& totalCells,
        int& activeCells,
        int& totalSamples) const;

    // clear all cached data
    void reset();

    // get cell for a point (for direct access)
    VisibilityCacheCell* getCell(const glm::vec3& point);

private:
    glm::vec3 m_sceneMin;
    glm::vec3 m_sceneMax;
    glm::vec3 m_cellSize{};
    glm::ivec3 m_resolution{};
    std::vector<VisibilityCacheCell> m_cells;

    // convert point to cell coordinates
    [[nodiscard]] glm::ivec3 pointToCellCoord(const glm::vec3& point) const;

    // convert cell coordinates to linear index
    [[nodiscard]] int coordToCellIndex(const glm::ivec3& coord) const;

    // check if cell coordinates are valid
    [[nodiscard]] bool isValidCoord(const glm::ivec3& coord) const;

    // get cell index for a point
    [[nodiscard]] int getCellIndex(const glm::vec3& point) const;
};

#endif //PROMETHEUS_SPATIAL_VISIBILITY_CACHE_H