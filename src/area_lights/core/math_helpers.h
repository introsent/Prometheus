//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_MATH_HELPERS_H
#define PROMETHEUS_MATH_HELPERS_H

#include <glm/glm.hpp>
#include <algorithm>
#include <cmath>

namespace MathHelpers {
    /// Safe mathematical operations
    // safe arccosine with clamping to avoid numerical issues
    inline float safeAcos(const float x) {
        return std::acos(std::clamp(x, -1.0f, 1.0f));
    }

    // safe square root, returns 0 for negative inputs
    inline float safeSqrt(const float x) {
        return std::sqrt(std::max(0.0f, x));
    }

    /// Gram-schmidt orthonormalization
    // creates a unit vector orthogonal to 'normal' from input vector 'v'
    // ref: spherical triangle sampling (Arvo, 1995)
    inline glm::vec3 gramSchmidtNormalize(const glm::vec3& v, const glm::vec3& normal) {
        const glm::vec3 orthogonal = v - glm::dot(v, normal) * normal;
        const float length = glm::length(orthogonal);
        return length > 1e-6f ? orthogonal / length : glm::vec3(0.0f);
    }

    /// Bounding box operations
    inline void computeTriangleBoundingBox(
        const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
        glm::vec3& bboxMin, glm::vec3& bboxMax) {

        bboxMin = glm::min(glm::min(v0, v1), v2);
        bboxMax = glm::max(glm::max(v0, v1), v2);
    }

    inline void expandBoundingBox(
        glm::vec3& bboxMin, glm::vec3& bboxMax,
        const glm::vec3& point) {

        bboxMin = glm::min(bboxMin, point);
        bboxMax = glm::max(bboxMax, point);
    }

    /// Geometric computations
    inline glm::vec3 computeTriangleCentroid(
        const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2) {

        return (v0 + v1 + v2) / 3.0f;
    }

    // compute total flux emitted by a triangle
    // flux = emission * intensity * area
    // for lambertian emitters, multiply by π for hemispherical emission
    inline float computeTriangleFlux(
        const glm::vec3& emission,
        const float intensity,
        const float area) {

        // use average of RGB components as scalar approximation
        const float avgEmission = (emission.r + emission.g + emission.b) / 3.0f;
        return avgEmission * intensity * area;
    }

    /// Barycentric coordinates
    // for point-in-triangle tests and interpolation
    struct BarycentricCoords {
        float u, v, w;

        [[nodiscard]] bool isInside() const {
            constexpr float eps = 1e-4f;
            return (u >= -eps) && (v >= -eps) && (w >= -eps);
        }
    };

    // compute barycentric coordinates of point p with respect to triangle (v0, v1, v2)
    inline BarycentricCoords computeBarycentricCoordinates(
        const glm::vec3& p,
        const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2) {

        const glm::vec3 v0v1 = v1 - v0;
        const glm::vec3 v0v2 = v2 - v0;
        const glm::vec3 v0p = p - v0;

        const float d00 = glm::dot(v0v1, v0v1);
        const float d01 = glm::dot(v0v1, v0v2);
        const float d11 = glm::dot(v0v2, v0v2);
        const float d20 = glm::dot(v0p, v0v1);
        const float d21 = glm::dot(v0p, v0v2);

        float denominator = d00 * d11 - d01 * d01;
        if (std::abs(denominator) < 1e-8f) {
            return {0.0f, 0.0f, 0.0f};
        }

        BarycentricCoords result{};
        result.v = (d11 * d20 - d01 * d21) / denominator;
        result.w = (d00 * d21 - d01 * d20) / denominator;
        result.u = 1.0f - result.v - result.w;

        return result;
    }
}

#endif //PROMETHEUS_MATH_HELPERS_H