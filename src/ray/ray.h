//
// Created by minaj on 10/27/2025.
//

#ifndef RAY_H
#define RAY_H
#include "glm/vec3.hpp"

struct Ray {
    glm::vec3 origin {0.f, 0.f, -1.f};
    glm::vec3 direction {0.f, 0.f, 1.f};
    float tNear = 0.f;
    float tFar = std::numeric_limits<float>::infinity();
    unsigned int mask = 0xFFFFFFFF;
};

#endif //RAY_H
