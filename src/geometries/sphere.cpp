//
// Created by minaj on 10/27/2025.
//

#include "sphere.h"

Sphere::Sphere(const glm::vec3& center, float radius, const EmbreeDevice* devicePtr)
    : Geometry(RTC_GEOMETRY_TYPE_SPHERE_POINT, devicePtr->getDevice()) {

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
}