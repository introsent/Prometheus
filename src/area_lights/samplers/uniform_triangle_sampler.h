//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_UNIFORM_TRIANGLE_SAMPLER_H
#define PROMETHEUS_UNIFORM_TRIANGLE_SAMPLER_H

#include "sampler_base.h"
#include "math_helpers.h"


/// Uniform triangle sampler
// samples points uniformly distributed across triangle surface
// uses square-root parameterization for area-preserving mapping

// ref: "Generating Random Points in Triangles" (Turk, Graphics Gems 1990)
// ref: "A Low-Distortion Map Between Triangle and Square" (Heitz, 2019)
class UniformTriangleSampler final : public IAreaLightSampler {
public:
    UniformTriangleSampler(
        const glm::vec3& v0,
        const glm::vec3& v1,
        const glm::vec3& v2,
        const glm::vec3& normal,
        float area,
        const glm::vec3& emission,
        float intensity);

    [[nodiscard]] AreaLightSample sample(
        const glm::vec3& shadingPoint,
        float u1, float u2) const override;

    [[nodiscard]] float pdf(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const override;

    void setIntensity(float intensity) override;
    [[nodiscard]] float getTotalFlux() const override;

private:
    glm::vec3 m_v0, m_v1, m_v2;  // triangle vertices
    glm::vec3 m_normal;           // surface normal
    float m_area;                 // triangle area
    glm::vec3 m_emission;         // emissive color
    float m_intensity;            // emission multiplier
};


#endif //PROMETHEUS_UNIFORM_TRIANGLE_SAMPLER_H