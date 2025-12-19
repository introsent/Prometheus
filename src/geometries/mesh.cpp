//
// Created by ivans on 19/12/2025.
//

#include "mesh.h"
#include <embree4/rtcore.h>

Mesh::Mesh(const std::vector<Vertex>& vertices,
           const std::vector<uint32_t>& indices,
           const EmbreeDevice* device)
    : Geometry(RTC_GEOMETRY_TYPE_TRIANGLE, device->handle()),
      m_originalVertices(vertices),
      m_indices(indices)
{
    uploadToEmbree(vertices, indices);

    // Initialize AABB
    updateAABB();
    updateTransformedAABB();
}

void Mesh::uploadToEmbree(const std::vector<Vertex>& vertices,
                          const std::vector<uint32_t>& indices)
{
    // Vertex buffer
    m_vertexBuffer = static_cast<float*>(
        rtcSetNewGeometryBuffer(
            m_geometry,
            RTC_BUFFER_TYPE_VERTEX,
            0,
            RTC_FORMAT_FLOAT3,
            3 * sizeof(float),
            vertices.size()
        )
    );

    for (size_t i = 0; i < vertices.size(); ++i) {
        m_vertexBuffer[i * 3 + 0] = vertices[i].position.x;
        m_vertexBuffer[i * 3 + 1] = vertices[i].position.y;
        m_vertexBuffer[i * 3 + 2] = vertices[i].position.z;
    }

    // Index buffer
    auto* idx = static_cast<uint32_t*>(
        rtcSetNewGeometryBuffer(
            m_geometry,
            RTC_BUFFER_TYPE_INDEX,
            0,
            RTC_FORMAT_UINT3,
            sizeof(uint32_t) * 3,
            indices.size() / 3
        )
    );

    std::memcpy(idx, indices.data(), indices.size() * sizeof(uint32_t));

    rtcCommitGeometry(m_geometry);
}

void Mesh::updateAABB() {
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

void Mesh::applyTransforms() {
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