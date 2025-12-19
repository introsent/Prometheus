//
// Created by minaj on 10/4/2025.
//
#include "triangle.h"

Triangle::Triangle(SceneManager* scene, const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr)
    : Geometry(scene, RTC_GEOMETRY_TYPE_TRIANGLE, devicePtr->handle()),
      m_originalVertices(vertices)
{
    m_vertexBuffer = static_cast<float*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), vertices.size()));
    fillVertexBuffer(vertices);

    const unsigned triCount = vertices.size() / 3;
    m_indexBuffer = static_cast<unsigned*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        3 * sizeof(unsigned), triCount));
    fillIndexBuffer(vertices);

    // Initialize AABB
    Triangle::updateAABB();
    updateTransformedAABB();
}

void Triangle::fillVertexBuffer(const std::vector<Vertex>& vertices) const {
    if (!m_vertexBuffer) return;
    for (size_t i = 0; i < vertices.size(); ++i) {
        m_vertexBuffer[i * 3 + 0] = vertices[i].position.x;
        m_vertexBuffer[i * 3 + 1] = vertices[i].position.y;
        m_vertexBuffer[i * 3 + 2] = vertices[i].position.z;
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

void Triangle::updateAABB() {
    if (m_originalVertices.empty()) return;

    m_minAABB = glm::vec3(m_originalVertices[0].position.x,
                          m_originalVertices[0].position.y,
                          m_originalVertices[0].position.z);
    m_maxAABB = m_minAABB;

    for (const auto& vertex : m_originalVertices) {
        const glm::vec3 pos(vertex.position.x, vertex.position.y, vertex.position.z);
        m_minAABB = glm::min(m_minAABB, pos);
        m_maxAABB = glm::max(m_maxAABB, pos);
    }
}

void Triangle::applyTransforms() {
    if (!m_vertexBuffer) return;

    const glm::mat4 finalTransform = getFinalTransform();

    // Apply transformation to all vertices
    for (size_t i = 0; i < m_originalVertices.size(); ++i) {
        const glm::vec3 originalPos(m_originalVertices[i].position.x,
                                   m_originalVertices[i].position.y,
                                   m_originalVertices[i].position.z);
        const glm::vec3 transformedPos = transformPoint(originalPos, finalTransform);

        m_vertexBuffer[i * 3 + 0] = transformedPos.x;
        m_vertexBuffer[i * 3 + 1] = transformedPos.y;
        m_vertexBuffer[i * 3 + 2] = transformedPos.z;
    }
}