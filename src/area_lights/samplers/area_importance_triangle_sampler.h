//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_AREA_IMPORTANCE_TRIANGLE_SAMPLER_H
#define PROMETHEUS_AREA_IMPORTANCE_TRIANGLE_SAMPLER_H

#include "sampler_base.h"
#include "math_helpers.h"

/// Area importance triangle sampler
// samples triangle proportional to subtended solid angle
// reduces variance when light covers small solid angle from shading point

// ref: "Stratified Sampling of Spherical Triangles" (Arvo, 1995)
// ref: "An Area-Preserving Parametrization for Spherical Rectangles" (Ureña et al., 2013)

// solid angle formula (van Oosterom & Strackee, 1983):
// omega = 2 * arctan( |a dot (b cross c)| / (1 + a dot b + b dot c + c dot a) )
// where a, b, c are unit vectors from shading point to triangle vertices
class AreaImportanceTriangleSampler final : public IAreaLightSampler {
public:
    AreaImportanceTriangleSampler(
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

    // calculate solid angle subtended by triangle as seen from point p
    [[nodiscard]] float calculateSolidAngle(const glm::vec3& p) const;

private:
    // sample direction to spherical triangle using Arvo's method
    // returns direction from shading point, pdf in solid angle measure
    [[nodiscard]] glm::vec3 sampleSphericalTriangle(
        const glm::vec3& shadingPoint,
        float u1, float u2,
        float* outPdfSolidAngle = nullptr) const;

    // intersect ray with triangle plane to get exact sample position
    [[nodiscard]] std::pair<float, glm::vec3> intersectTrianglePlane(
        const glm::vec3& rayOrigin,
        const glm::vec3& rayDirection) const;

    glm::vec3 m_v0, m_v1, m_v2;
    glm::vec3 m_normal;
    float m_area;
    glm::vec3 m_emission;
    float m_intensity;
};


#endif //PROMETHEUS_AREA_IMPORTANCE_TRIANGLE_SAMPLER_H