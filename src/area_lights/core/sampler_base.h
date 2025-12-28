//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_SAMPLER_BASE_H
#define PROMETHEUS_SAMPLER_BASE_H

#include <glm/glm.hpp>

/// Area light sample structure
// contains all information about a sampled point on an area light
struct AreaLightSample {
    glm::vec3 position;     // sampled position on light surface
    glm::vec3 normal;       // surface normal at sample point
    float pdf;              // probability density function value
    float misWeight;        // multiple importance sampling weight (1.0 if not using MIS)
    glm::vec3 radiance;     // emitted radiance (emission * intensity)
    float area;             // total area of the light source
};


/// Base sampler interface
// all sampling strategies implement this interface
// ref: "Importance Sampling of Many Lights" (Estevez & Kulla, 2018)
class IAreaLightSampler {
public:
    virtual ~IAreaLightSampler() = default;

    // sample a point on the light surface given a shading point
    // u1, u2: uniform random numbers in [0,1]
    // returns: sample information including position, normal, and pdf
    [[nodiscard]] virtual AreaLightSample sample(
        const glm::vec3& shadingPoint,
        float u1, float u2) const = 0;

    // evaluate pdf for a given light point
    // used for multiple importance sampling
    [[nodiscard]] virtual float pdf(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const = 0;

    // update light intensity (for dynamic lights)
    virtual void setIntensity(float intensity) = 0;

    // compute total flux emitted by the light
    // flux = average(emission) * intensity * area * pi (for lambertian)
    [[nodiscard]] virtual float getTotalFlux() const = 0;
};


/// Sampling strategy enumeration
enum class SamplingStrategy {
    Uniform,                        // uniform area sampling
    AreaImportance,                 // solid angle importance sampling
    HierarchicalFlux,               // BVH-based flux importance sampling
    VisibilityAwareHierarchical     // visibility-learning hierarchical sampling
};

#endif //PROMETHEUS_SAMPLER_BASE_H