//
// Created by minaj on 10/27/2025.
//

#ifndef LIGHT_H
#define LIGHT_H
#include "glm/vec3.hpp"

class TriangleAreaLight;
class MeshAreaLight;

enum class LightType {
    Point,
    Directional,
    TriangleArea,
    MeshArea
};

/// Base light structure for simple lights
struct Light {
    glm::vec3 origin{};
    glm::vec3 direction{};
    glm::vec3 color{};
    float intensity{};
    LightType type{};

    // for Area Lights store pointer to area light object
    TriangleAreaLight* triangleAreaLight{nullptr};
    MeshAreaLight* meshAreaLight{nullptr};

    // check if this is an area light
    [[nodiscard]] bool isAreaLight() const {
        return type == LightType::TriangleArea || type == LightType::MeshArea;
    }
};
#endif //LIGHT_H
