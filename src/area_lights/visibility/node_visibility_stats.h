//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_NODE_VISIBILITY_STATS_H
#define PROMETHEUS_NODE_VISIBILITY_STATS_H
#include <atomic>


/// Visibility statistics for a single BVH node
// tracks number of visible vs total samples
struct NodeVisibilityStats {
    int nodeIndex{-1};
    std::atomic<uint32_t> visibleSamples{0};
    std::atomic<uint32_t> totalSamples{0};

    NodeVisibilityStats() = default;

    // atomics are not copyable, only movable
    NodeVisibilityStats(const NodeVisibilityStats&) = delete;
    NodeVisibilityStats& operator=(const NodeVisibilityStats&) = delete;

    NodeVisibilityStats(NodeVisibilityStats&& other) noexcept
        : nodeIndex(other.nodeIndex)
        , visibleSamples(other.visibleSamples.load(std::memory_order_relaxed))
        , totalSamples(other.totalSamples.load(std::memory_order_relaxed))
    {}

    NodeVisibilityStats& operator=(NodeVisibilityStats&& other) noexcept {
        if (this != &other) {
            nodeIndex = other.nodeIndex;
            visibleSamples.store(
                other.visibleSamples.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
            totalSamples.store(
                other.totalSamples.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
        }
        return *this;
    }

    // get estimated visibility probability
    // returns 1.0 if no samples recorded yet (optimistic default)
    [[nodiscard]] float getVisibilityProbability() const {
        const uint32_t total = totalSamples.load(std::memory_order_relaxed);
        if (total == 0) return 1.0f;

        const uint32_t visible = visibleSamples.load(std::memory_order_relaxed);
        return static_cast<float>(visible) / static_cast<float>(total);
    }

    // record a new sample
    void recordSample(bool wasVisible) {
        totalSamples.fetch_add(1, std::memory_order_relaxed);
        if (wasVisible) {
            visibleSamples.fetch_add(1, std::memory_order_relaxed);
        }
    }
};

#endif //PROMETHEUS_NODE_VISIBILITY_STATS_H