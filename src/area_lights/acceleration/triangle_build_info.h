//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_TRIANGLE_BUILD_INFO_H
#define PROMETHEUS_TRIANGLE_BUILD_INFO_H
#include <glm/vec3.hpp>

/// Triangle information for BVH construction
struct TriangleBuildInfo {
    glm::vec3 centroid;      // triangle centroid for spatial partitioning
    glm::vec3 bboxMin;       // triangle bounding box
    glm::vec3 bboxMax;
    float flux;              // emitted flux
    float area;              // surface area
    int originalIndex;       // index in original triangle array

    TriangleBuildInfo(
        const glm::vec3& cent,
        const glm::vec3& min,
        const glm::vec3& max,
        float f, float a, int idx)
        : centroid(cent)
        , bboxMin(min)
        , bboxMax(max)
        , flux(f)
        , area(a)
        , originalIndex(idx)
    {}
};

#endif //PROMETHEUS_TRIANGLE_BUILD_INFO_H