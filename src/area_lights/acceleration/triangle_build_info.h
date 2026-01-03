//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_TRIANGLE_BUILD_INFO_H
#define PROMETHEUS_TRIANGLE_BUILD_INFO_H
#include <glm/vec3.hpp>

struct TriangleBuildInfo {
    glm::vec3 centroid;
    glm::vec3 bboxMin;
    glm::vec3 bboxMax;
    glm::vec3 normal;
    float flux;
    float area;
    int originalIndex;

    TriangleBuildInfo(
        const glm::vec3& cent,
        const glm::vec3& min,
        const glm::vec3& max,
        const glm::vec3& norm,  // NEW
        float f, float a, int idx)
        : centroid(cent)
        , bboxMin(min)
        , bboxMax(max)
        , normal(norm)
        , flux(f)
        , area(a)
        , originalIndex(idx)
    {}
};

#endif //PROMETHEUS_TRIANGLE_BUILD_INFO_H