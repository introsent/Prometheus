//
// Created by minaj on 10/27/2025.
//

#ifndef CAMERA_H
#define CAMERA_H
#include "glm/vec3.hpp"
#include "../ray/ray.h"

class Camera {
public:
    Camera(const glm::vec3& position, float fovDegrees, float aspectRatio);

    [[nodiscard]] Ray generateRay(float u, float v) const;

    void setPosition(const glm::vec3& pos) { m_position = pos; }
    void setFOV(float fovDegrees);

private:
    glm::vec3 m_position;
    glm::vec3 m_forward;
    glm::vec3 m_right;
    glm::vec3 m_up;
    float m_aspectRatio;
    float m_viewportHeight {};
    float m_viewportWidth  {};
};




#endif //CAMERA_H
