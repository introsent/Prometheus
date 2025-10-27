//
// Created by minaj on 10/27/2025.
//

#ifndef LIGHT_H
#define LIGHT_H
#include "glm/vec3.hpp"

enum class LightType
{
    Point,
    Directional
};

struct Light
{
    glm::vec3 origin{};
    glm::vec3 direction{};
    glm::vec3 color{};
    float intensity{};
    LightType type{};
};
#endif //LIGHT_H
