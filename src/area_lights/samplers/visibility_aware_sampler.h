//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_VISIBILITY_AWARE_SAMPLER_H
#define PROMETHEUS_VISIBILITY_AWARE_SAMPLER_H
#include "lights/mesh_area_light.h"
#include "mis/mis_weights.h"

/// Visibility-aware hierarchical sampler
///
// combines:
// - hierarchical importance sampling (flux / distance^2)
// - learned visibility information
// - multiple importance sampling with BSDF

// ref: "Bayesian Online Regression for Adaptive Direct Illumination Sampling" (Vévoda et al., 2018)
class VisibilityAwareHierarchicalSampler {
public:
    struct Configuration {
        bool enableVisibilityLearning{true};
        bool enableMIS{true};
        float visibilityWeight{0.5f};        // 0=pure flux, 1=pure visibility
        int minSamplesForLearning{4};        // min samples before using visibility
        MISHeuristic misHeuristic{MISHeuristic::Balance};
    };

    struct Sample {
        glm::vec3 position;
        glm::vec3 normal;
        glm::vec3 radiance;
        float pdf;
        float misWeight;
        bool isValid;
    };

    VisibilityAwareHierarchicalSampler(
        MeshAreaLight* light,
        SpatialVisibilityCache* visCache,
        const Configuration& config = Configuration());

    // sample light using visibility-aware importance
    [[nodiscard]] Sample sampleLight(
        const glm::vec3& shadingPoint,
        const glm::vec3& normal,
        float u1, float u2, float u3) const;

    // evaluate pdf for a given light point
    [[nodiscard]] float evaluatePdf(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    // record visibility for learning
    void recordVisibilitySample(
        const glm::vec3& shadingPoint,
        int nodeIndex,
        bool wasVisible);

    // configuration access
    [[nodiscard]] const Configuration& getConfig() const { return m_config; }
    void setConfig(const Configuration& config) { m_config = config; }

private:
    MeshAreaLight* m_light;
    SpatialVisibilityCache* m_visCache;
    Configuration m_config;

    // hierarchical selection with visibility weighting
    int selectBVHNodeWithVisibility(
        int nodeIndex,
        const glm::vec3& shadingPoint,
        float u,
        float& outPdf) const;

    // calculate node importance combining flux, distance, and visibility
    [[nodiscard]] float calculateNodeImportance(
        int nodeIndex,
        const glm::vec3& shadingPoint) const;
};
#endif //PROMETHEUS_VISIBILITY_AWARE_SAMPLER_H