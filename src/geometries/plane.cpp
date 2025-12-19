//
// Created by minaj on 10/27/2025.
//

#include "plane.h"
#include <limits>
#include "glm/ext/quaternion_geometric.hpp"

Plane::Plane(SceneManager* scene, const glm::vec3& origin, const glm::vec3& normal, const EmbreeDevice* devicePtr)
    : Geometry(scene, RTC_GEOMETRY_TYPE_QUAD, devicePtr->handle()),
      m_origin(origin),
      m_normal(glm::normalize(normal))
{
    // Create a large quad to represent the plane (ideally infinite, but we simulate it)
    constexpr float planeSize = 10000.f;

    // Calculate perpendicular vectors for the plane
    glm::vec3 u;
    if (std::abs(normal.x) > 0.1f) {
        u = glm::normalize(glm::cross(glm::vec3(0, 1, 0), normal));
    } else {
        u = glm::normalize(glm::cross(glm::vec3(1, 0, 0), normal));
    }
    const glm::vec3 v = glm::normalize(glm::cross(normal, u));

    // Create quad vertices
    m_vertexBuffer = static_cast<float*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), 4));

    // Store original vertices
    m_originalVertices[0] = origin - u * planeSize - v * planeSize;
    m_originalVertices[1] = origin + u * planeSize - v * planeSize;
    m_originalVertices[2] = origin + u * planeSize + v * planeSize;
    m_originalVertices[3] = origin - u * planeSize + v * planeSize;

    // Fill vertex buffer
    for (int i = 0; i < 4; ++i) {
        m_vertexBuffer[i * 3 + 0] = m_originalVertices[i].x;
        m_vertexBuffer[i * 3 + 1] = m_originalVertices[i].y;
        m_vertexBuffer[i * 3 + 2] = m_originalVertices[i].z;
    }

    // Create quad indices
    auto* indices = static_cast<unsigned*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4,
        4 * sizeof(unsigned), 1));

    indices[0] = 0;
    indices[1] = 1;
    indices[2] = 2;
    indices[3] = 3;

    // Initialize AABB
    Plane::updateAABB();
    updateTransformedAABB();
}

void Plane::updateAABB() {
    if (!m_vertexBuffer) return;

    m_minAABB = m_originalVertices[0];
    m_maxAABB = m_originalVertices[0];

    for (int i = 1; i < 4; ++i) {
        m_minAABB = glm::min(m_minAABB, m_originalVertices[i]);
        m_maxAABB = glm::max(m_maxAABB, m_originalVertices[i]);
    }
}

void Plane::applyTransforms() {
    if (!m_vertexBuffer) return;

    const glm::mat4 finalTransform = getFinalTransform();

    // Apply transformation to all 4 vertices
    for (int i = 0; i < 4; ++i) {
        const glm::vec3 transformedPos = transformPoint(m_originalVertices[i], finalTransform);
        m_vertexBuffer[i * 3 + 0] = transformedPos.x;
        m_vertexBuffer[i * 3 + 1] = transformedPos.y;
        m_vertexBuffer[i * 3 + 2] = transformedPos.z;
    }
}