//
// Created by minaj on 10/4/2025.
//

#include "geometry.h"

Geometry::Geometry(RTCGeometryType geometryType, RTCDevice device)
    : m_type(geometryType)
{
    m_geometry = rtcNewGeometry(device, geometryType);
}

RTCGeometryType Geometry::getType() const {
    return m_type;
}

RTCGeometry Geometry::getGeometry() const {
    return m_geometry;
}

void Geometry::commit() const {
    rtcCommitGeometry(m_geometry);
}

unsigned Geometry::attach(const EmbreeScene* embreeScenePtr) const {
    return rtcAttachGeometry(embreeScenePtr->handle(), m_geometry);
}

void Geometry::release() const {
    rtcReleaseGeometry(m_geometry);
}

void Geometry::translate(const glm::vec3& translation) {
    m_translationMatrix = glm::translate(glm::mat4(1.0f), translation);
}

void Geometry::rotateY(float yawRadians) {
    m_rotationMatrix = glm::rotate(glm::mat4(1.0f), yawRadians, glm::vec3(0.f, 1.f, 0.f));
}

void Geometry::rotate(float angle, const glm::vec3& axis) {
    m_rotationMatrix = glm::rotate(glm::mat4(1.0f), angle, axis);
}

void Geometry::scale(const glm::vec3& scaleVec) {
    m_scaleMatrix = glm::scale(glm::mat4(1.0f), scaleVec);
}

void Geometry::updateTransforms() {
    applyTransforms();
    updateTransformedAABB();

    rtcUpdateGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_VERTEX, 0);
    rtcCommitGeometry(m_geometry);
}

void Geometry::updateTransformedAABB() {
    const glm::mat4 finalTransform = getFinalTransform();

    // Transform all 8 corners of the AABB
    const glm::vec3 corners[8] = {
        glm::vec3(m_minAABB.x, m_minAABB.y, m_minAABB.z),
        glm::vec3(m_maxAABB.x, m_minAABB.y, m_minAABB.z),
        glm::vec3(m_maxAABB.x, m_minAABB.y, m_maxAABB.z),
        glm::vec3(m_minAABB.x, m_minAABB.y, m_maxAABB.z),
        glm::vec3(m_minAABB.x, m_maxAABB.y, m_minAABB.z),
        glm::vec3(m_maxAABB.x, m_maxAABB.y, m_minAABB.z),
        glm::vec3(m_maxAABB.x, m_maxAABB.y, m_maxAABB.z),
        glm::vec3(m_minAABB.x, m_maxAABB.y, m_maxAABB.z)
    };

    m_transformedMinAABB = transformPoint(corners[0], finalTransform);
    m_transformedMaxAABB = m_transformedMinAABB;

    for (int i = 1; i < 8; ++i) {
        const glm::vec3 transformedCorner = transformPoint(corners[i], finalTransform);
        m_transformedMinAABB = glm::min(m_transformedMinAABB, transformedCorner);
        m_transformedMaxAABB = glm::max(m_transformedMaxAABB, transformedCorner);
    }
}

glm::vec3 Geometry::transformPoint(const glm::vec3& point, const glm::mat4& matrix) {
    const glm::vec4 homogeneous(point, 1.0f);
    const glm::vec4 transformed = matrix * homogeneous;
    return glm::vec3(transformed) / transformed.w;
}

glm::mat4 Geometry::getFinalTransform() const {
    return m_translationMatrix * m_rotationMatrix * m_scaleMatrix;
}