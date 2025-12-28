//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_VISIBILITY_CACHE_CELL_H
#define PROMETHEUS_VISIBILITY_CACHE_CELL_H
#include <mutex>
#include <vector>
#include <glm/vec3.hpp>
#include "node_visibility_stats.h"

/// Cache cell
// stores visibility statistics for all BVH nodes from this spatial region
struct VisibilityCacheCell {
    glm::vec3 center{0.0f};
    float radius{0.0f};
    std::vector<NodeVisibilityStats> nodeStats;
    std::mutex cellMutex;  // protects nodeStats vector

    VisibilityCacheCell() = default;

    // mutexes are not copyable/movable
    VisibilityCacheCell(const VisibilityCacheCell&) = delete;
    VisibilityCacheCell& operator=(const VisibilityCacheCell&) = delete;

    VisibilityCacheCell(VisibilityCacheCell&& other) noexcept
        : center(other.center)
        , radius(other.radius)
        , nodeStats(std::move(other.nodeStats))
        , cellMutex()  // construct new mutex
    {}

    VisibilityCacheCell& operator=(VisibilityCacheCell&& other) noexcept {
        if (this != &other) {
            center = other.center;
            radius = other.radius;
            nodeStats = std::move(other.nodeStats);
        }
        return *this;
    }

    // find existing node statistics (read-only)
    [[nodiscard]] NodeVisibilityStats* getNodeStats(int nodeIndex) {
        for (auto& stats : nodeStats) {
            if (stats.nodeIndex == nodeIndex) {
                return &stats;
            }
        }
        return nullptr;
    }

    // get or create node statistics (thread-safe)
    NodeVisibilityStats* getOrCreateNodeStats(int nodeIndex) {
        std::lock_guard<std::mutex> lock(cellMutex);

        // check if already exists
        for (auto& stats : nodeStats) {
            if (stats.nodeIndex == nodeIndex) {
                return &stats;
            }
        }

        // create new entry
        nodeStats.emplace_back();
        nodeStats.back().nodeIndex = nodeIndex;
        return &nodeStats.back();
    }
};

#endif //PROMETHEUS_VISIBILITY_CACHE_CELL_H