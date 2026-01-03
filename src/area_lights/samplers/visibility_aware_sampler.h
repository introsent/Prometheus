//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_VISIBILITY_AWARE_SAMPLER_H
#define PROMETHEUS_VISIBILITY_AWARE_SAMPLER_H

#include <vector>
#include <glm/glm.hpp>
#include "acceleration/bvh_node.h"
#include "mis/mis_weights.h"

class MeshAreaLight;
class SpatialVisibilityCache;

/// Visibility-aware hierarchical light sampler
/// Implements techniques from:
/// ref: Conty & Kulla "Importance Sampling of Many Lights with Adaptive Tree Splitting" (HPG 2018)
/// ref: Vévoda et al. "Bayesian online regression for adaptive direct illumination sampling" (SIGGRAPH 2018)
class VisibilityAwareHierarchicalSampler {
public:
    struct Configuration {
        bool enableVisibilityLearning{true};
        bool enableAdaptiveSplitting{true};
        bool enableMIS{true};

        float visibilityWeight{0.7f};
        float splitThreshold{0.02f};
        int minSamplesForLearning{8};
        int maxSplitDepth{4};

        MISHeuristic misHeuristic{MISHeuristic::Balance};
    };

    struct Sample {
        glm::vec3 position{0.0f};
        glm::vec3 normal{0.0f};
        glm::vec3 radiance{0.0f};
        float pdf{0.0f};
        float misWeight{1.0f};
        bool isValid{false};

        // store traversal path in sample result for thread-safe visibility recording
        std::vector<int> traversalPath;
    };

    VisibilityAwareHierarchicalSampler(
        MeshAreaLight* light,
        SpatialVisibilityCache* visCache,
        const Configuration& config = Configuration());

    Sample sampleLight(
        const glm::vec3& shadingPoint,
        const glm::vec3& normal,
        float u1, float u2, float u3) const;

    float evaluatePdf(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    void recordVisibilitySample(
        const glm::vec3& shadingPoint,
        int nodeIndex,
        bool wasVisible);

    // Record visibility for entire traversal path (thread-safe version)
    void recordTraversalVisibility(
        const glm::vec3& shadingPoint,
        const std::vector<int>& traversalPath,
        bool wasVisible);

private:
    MeshAreaLight* m_light;
    SpatialVisibilityCache* m_visCache;
    Configuration m_config;

    /// Importance calculation (Equation 3 from paper)
    float calculateNodeImportance(
        int nodeIndex,
        const glm::vec3& shadingPoint,
        const glm::vec3& shadingNormal) const;

    /// Adaptive splitting decision (Section 5.4)
    bool shouldSplit(
        int nodeIndex,
        const glm::vec3& shadingPoint,
        int currentDepth) const;

    /// Compute split variance score (Equation 10 from paper)
    float computeSplitVarianceScore(
        int nodeIndex,
        const glm::vec3& shadingPoint) const;

    /// Traversal result structure
    struct TraversalResult {
        int leafNodeIndex{-1};
        float pathPdf{1.0f};
        std::vector<std::pair<int, float>> sampledLeaves;
    };

    /// Main traversal with adaptive splitting
    TraversalResult traverseWithAdaptiveSplitting(
        int nodeIndex,
        const glm::vec3& shadingPoint,
        const glm::vec3& shadingNormal,
        float u,
        int depth,
        std::vector<int>& outTraversalPath) const;

    /// Single-path stochastic traversal (fallback)
    int selectBVHNodeStochastic(
        int nodeIndex,
        const glm::vec3& shadingPoint,
        const glm::vec3& shadingNormal,
        float u,
        float& outPdf,
        std::vector<int>& outTraversalPath) const;
};

#endif // PROMETHEUS_VISIBILITY_AWARE_SAMPLER_H