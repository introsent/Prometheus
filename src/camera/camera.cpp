//
// Created by minaj on 10/27/2025.
//

#include "camera.h"
#include <cmath>

#include "glm/ext/quaternion_geometric.hpp"

Camera::Camera(const glm::vec3& position, float fovDegrees, float aspectRatio)
    : m_position(position), m_aspectRatio(aspectRatio),
      m_forward(0.f, 0.f, 1.f), m_right(1.f, 0.f, 0.f), m_up(0.f, 1.f, 0.f) {
    setFOV(fovDegrees);
}

void Camera::setFOV(const float fovDegrees) {
    const float fovRadians = fovDegrees * 3.14159265f / 180.f;
    m_viewportHeight = 2.f * std::tan(fovRadians / 2.f);
    m_viewportWidth = m_viewportHeight * m_aspectRatio;
}

Ray Camera::generateRay(const float u, const float v) const {
    const float x = (u - 0.5f) * m_viewportWidth;
    const float y = (0.5f - v) * m_viewportHeight; // Flip Y for screen coordinates

    const glm::vec3 direction = glm::normalize(
        m_forward + x * m_right + y * m_up
    );

    return Ray{m_position, direction};
}