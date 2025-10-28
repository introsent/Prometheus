//
// Created by minaj on 10/4/2025.
//
#include "triangle.h"

Triangle::Triangle(const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr,
                   CullingMode culling)
    : Geometry(RTC_GEOMETRY_TYPE_TRIANGLE, devicePtr->handle())
    , m_cullingMode(culling) {

    m_vertexBuffer = static_cast<float*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), vertices.size()));
    fillVertexBuffer(vertices);

    const unsigned triCount = vertices.size() / 3;
    m_indexBuffer = static_cast<unsigned*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        3 * sizeof(unsigned), triCount));
    fillIndexBuffer(vertices);

    m_originalVertices = vertices;

    // Initialize AABB
    updateAABB();
    updateTransformedAABB();
}

void Triangle::fillVertexBuffer(const std::vector<Vertex>& vertices) const {
    if (!m_vertexBuffer) return;
    for (size_t i = 0; i < vertices.size(); ++i) {
        m_vertexBuffer[i * 3 + 0] = vertices[i].x;
        m_vertexBuffer[i * 3 + 1] = vertices[i].y;
        m_vertexBuffer[i * 3 + 2] = vertices[i].z;
    }
}

void Triangle::fillIndexBuffer(const std::vector<Vertex>& vertices) const {
    const auto triCount = static_cast<unsigned>(vertices.size() / 3);
    for (unsigned t = 0; t < triCount; ++t) {
        m_indexBuffer[t * 3 + 0] = t * 3 + 0;
        m_indexBuffer[t * 3 + 1] = t * 3 + 1;
        m_indexBuffer[t * 3 + 2] = t * 3 + 2;
    }
}

CullingMode Triangle::getCullingMode() const {
    return m_cullingMode;
}

void Triangle::translate(const glm::vec3 &translation) {
    m_translationMatrix = glm::translate(glm::mat4(1.0f), translation);
}

void Triangle::rotateY(float yawRadians) {
    m_rotationMatrix = glm::rotate(glm::mat4(1.0f), yawRadians, glm::vec3(0.f, 1.f, 0.f));
}

void Triangle::scale(const glm::vec3 &scaleVec) {
    m_scaleMatrix = glm::scale(glm::mat4(1.0f), scaleVec);
}

void Triangle::updateTransforms() {
    applyTransforms();
    updateTransformedAABB();

    rtcUpdateGeometryBuffer(m_geometry, RTC_BUFFER_TYPE_VERTEX, 0);
    rtcCommitGeometry(m_geometry);
}

void Triangle::updateAABB() {
    if (m_originalVertices.empty()) return;

    m_minAABB = glm::vec3(m_originalVertices[0].x,
                          m_originalVertices[0].y,
                          m_originalVertices[0].z);
    m_maxAABB = m_minAABB;

    for (const auto& vertex : m_originalVertices) {
        const glm::vec3 pos(vertex.x, vertex.y, vertex.z);
        m_minAABB = glm::min(m_minAABB, pos);
        m_maxAABB = glm::max(m_maxAABB, pos);
    }
}

void Triangle::updateTransformedAABB() {
    const glm::mat4 finalTransform = m_translationMatrix * m_rotationMatrix * m_scaleMatrix;

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

void Triangle::applyTransforms() const {
    // Combine transformations: Scale -> Rotate -> Translate
    const glm::mat4 finalTransform = m_translationMatrix * m_rotationMatrix * m_scaleMatrix;

    // Apply transformation to all vertices
    for (size_t i = 0; i < m_originalVertices.size(); ++i) {
        const glm::vec3 originalPos(m_originalVertices[i].x,
                                   m_originalVertices[i].y,
                                   m_originalVertices[i].z);
        const glm::vec3 transformedPos = transformPoint(originalPos, finalTransform);

        m_vertexBuffer[i * 3 + 0] = transformedPos.x;
        m_vertexBuffer[i * 3 + 1] = transformedPos.y;
        m_vertexBuffer[i * 3 + 2] = transformedPos.z;
    }
}

glm::vec3 Triangle::transformPoint(const glm::vec3 &point, const glm::mat4 &matrix) {
    const glm::vec4 homogeneous(point, 1.0f);
    const glm::vec4 transformed = matrix * homogeneous;
    return glm::vec3(transformed) / transformed.w;
}
