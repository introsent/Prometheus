//
// Created by ivans on 25/12/2025.
//

#ifndef PROMETHEUS_MATH_HELPERS_H
#define PROMETHEUS_MATH_HELPERS_H

#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

namespace MathHelpers
{
    // clamp-safe acos
    inline float safeAcos(float x)
    {
        return std::acos(std::clamp(x, -1.0f, 1.0f));
    }

    // clamp-safe sqrt
    inline float safeSqrt(float x)
    {
        return std::sqrt(std::max(0.0f, x));
    }

    // Gram-Schmidt: return (v - proj_a(v)) normalized
    // if near zero, return any vector perpendicular to a
    inline glm::vec3 gramSchmidtNormalize(const glm::vec3& v, const glm::vec3& a)
    {
        const glm::vec3 proj = a * glm::dot(v, a);
        glm::vec3 w = v - proj;

        float len2 = glm::dot(w, w);
        if (len2 <= 1e-10f)
        {
            // fallback: find any vector perpendicular to a
            const glm::vec3 absA = glm::abs(a);
            const glm::vec3 other =
                (absA.x <= absA.y && absA.x <= absA.z) ? glm::vec3(1, 0, 0) :
                (absA.y <= absA.x && absA.y <= absA.z) ? glm::vec3(0, 1, 0) :
                                                        glm::vec3(0, 0, 1);

            w = glm::cross(a, other);
            len2 = glm::dot(w, w);

            // last resort (almost parallel to all axes)
            if (len2 <= 1e-10f)
                return glm::normalize(glm::cross(a, glm::vec3(1, 0, 0) + 0.0001f));
        }

        return glm::normalize(w);
    }

    // compute triangle axis-aligned bounding box
    inline void computeTriangleBoundingBox(
        const glm::vec3& v0,
        const glm::vec3& v1,
        const glm::vec3& v2,
        glm::vec3& bboxMin,
        glm::vec3& bboxMax)
    {
        bboxMin = glm::min(v0, glm::min(v1, v2));
        bboxMax = glm::max(v0, glm::max(v1, v2));
    }

    // compute triangle centroid
    inline glm::vec3 computeTriangleCentroid(
        const glm::vec3& v0,
        const glm::vec3& v1,
        const glm::vec3& v2)
    {
        return (v0 + v1 + v2) / 3.0f;
    }

    // compute triangle flux (power)
    inline float computeTriangleFlux(
        const glm::vec3& emission,
        float intensity,
        float area)
    {
        const float avgEmission =
            (emission.r + emission.g + emission.b) / 3.0f;

        return avgEmission * intensity * area;
    }

    // expand bounding box to include a point
    inline void expandBoundingBox(
        glm::vec3& bboxMin,
        glm::vec3& bboxMax,
        const glm::vec3& point)
    {
        bboxMin = glm::min(bboxMin, point);
        bboxMax = glm::max(bboxMax, point);
    }
}

#endif //PROMETHEUS_MATH_HELPERS_H