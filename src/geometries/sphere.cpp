//
// Created by minaj on 10/27/2025.
//

#include "sphere.h"
#include "init/embree_device.h"

Sphere::Sphere(const glm::vec3& center, float radius, const EmbreeDevice* devicePtr)
    : Geometry(RTC_GEOMETRY_TYPE_SPHERE_POINT, devicePtr->handle()),
      m_center(center),
      m_radius(radius)
{
    struct Sphere4 {
        float x, y, z, r;
    };

    auto* spheres = static_cast<Sphere4*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4,
        sizeof(Sphere4), 1));

    spheres[0].x = center.x;
    spheres[0].y = center.y;
    spheres[0].z = center.z;
    spheres[0].r = radius;

    // Initialize AABB
    Sphere::updateAABB();
    updateTransformedAABB();
}

void Sphere::updateAABB() {
    // Sphere AABB is simply center +/- radius in all directions
    m_minAABB = m_center - glm::vec3(m_radius);
    m_maxAABB = m_center + glm::vec3(m_radius);
}

void Sphere::applyTransforms() {
    const glm::mat4 finalTransform = getFinalTransform();

    // Transform the center
    const glm::vec3 transformedCenter = transformPoint(m_center, finalTransform);

    // Extract scale from transformation matrix to adjust radius
    // Use the maximum scale component since sphere should scale uniformly
    const glm::vec3 scaleVec = glm::vec3(
        glm::length(glm::vec3(finalTransform[0])),
        glm::length(glm::vec3(finalTransform[1])),
        glm::length(glm::vec3(finalTransform[2]))
    );
    const float maxScale = glm::max(glm::max(scaleVec.x, scaleVec.y), scaleVec.z);
    const float transformedRadius = m_radius * maxScale;

    // Update Embree buffer
    struct Sphere4 {
        float x, y, z, r;
    };

    auto* spheres = static_cast<Sphere4*>(rtcGetGeometryBufferData(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0));

    if (spheres) {
        spheres[0].x = transformedCenter.x;
        spheres[0].y = transformedCenter.y;
        spheres[0].z = transformedCenter.z;
        spheres[0].r = transformedRadius;
    }
}