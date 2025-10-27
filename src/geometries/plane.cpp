//
// Created by minaj on 10/27/2025.
//

#include "plane.h"
#include <limits>

#include "glm/ext/quaternion_geometric.hpp"

Plane::Plane(const glm::vec3& origin, const glm::vec3& normal, const EmbreeDevice* devicePtr)
    : Geometry(RTC_GEOMETRY_TYPE_QUAD, devicePtr->getDevice()) {

    // Create a large quad to represent the plane (ideally infinite, but we simulate it)
    const float planeSize = 10000.f;

    // Calculate perpendicular vectors for the plane
    glm::vec3 u, v;
    if (std::abs(normal.x) > 0.1f) {
        u = glm::normalize(glm::cross(glm::vec3(0, 1, 0), normal));
    } else {
        u = glm::normalize(glm::cross(glm::vec3(1, 0, 0), normal));
    }
    v = glm::normalize(glm::cross(normal, u));

    // Create quad vertices
    auto* vertices = static_cast<float*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), 4));

    glm::vec3 p0 = origin - u * planeSize - v * planeSize;
    glm::vec3 p1 = origin + u * planeSize - v * planeSize;
    glm::vec3 p2 = origin + u * planeSize + v * planeSize;
    glm::vec3 p3 = origin - u * planeSize + v * planeSize;

    vertices[0] = p0.x; vertices[1] = p0.y; vertices[2] = p0.z;
    vertices[3] = p1.x; vertices[4] = p1.y; vertices[5] = p1.z;
    vertices[6] = p2.x; vertices[7] = p2.y; vertices[8] = p2.z;
    vertices[9] = p3.x; vertices[10] = p3.y; vertices[11] = p3.z;

    // Create quad indices
    auto* indices = static_cast<unsigned*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4,
        4 * sizeof(unsigned), 1));

    indices[0] = 0;
    indices[1] = 1;
    indices[2] = 2;
    indices[3] = 3;
}
