//
// Created by ivans on 28/12/2025.
//

#include "uniform_triangle_sampler.h"
#include <cmath>

UniformTriangleSampler::UniformTriangleSampler(
    const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
    const glm::vec3& normal, float area,
    const glm::vec3& emission, float intensity)
    : m_v0(v0), m_v1(v1), m_v2(v2)
    , m_normal(normal)
    , m_area(area)
    , m_emission(emission)
    , m_intensity(intensity)
{}

AreaLightSample UniformTriangleSampler::sample(
    const glm::vec3& shadingPoint,
    float u1, float u2) const {

    // Square-root parameterization for uniform sampling
    // maps unit square [0,1]^2 to triangle with area-preserving transformation
    // the sqrt(u1) term stretches the parameterization to maintain uniform density
    const float sqrtU1 = std::sqrt(u1);
    const float b0 = 1.0f - sqrtU1;
    const float b1 = sqrtU1 * (1.0f - u2);
    const float b2 = sqrtU1 * u2;

    const glm::vec3 position = b0 * m_v0 + b1 * m_v1 + b2 * m_v2;

    // for uniform area sampling, pdf is constant: 1/area
    const float pdf = 1.0f / m_area;

    const glm::vec3 radiance = m_emission * m_intensity;

    return AreaLightSample{
        position,
        m_normal,
        pdf,
        1.0f,  // no MIS weight for uniform sampling
        radiance,
        m_area
    };
}

float UniformTriangleSampler::pdf(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    // uniform sampling has constant pdf = 1/area
    return 1.0f / m_area;
}

void UniformTriangleSampler::setIntensity(float intensity) {
    m_intensity = intensity;
}

float UniformTriangleSampler::getTotalFlux() const {
    return MathHelpers::computeTriangleFlux(m_emission, m_intensity, m_area);
}