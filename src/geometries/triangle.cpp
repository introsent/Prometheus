//
// Created by minaj on 10/4/2025.
//

#include "triangle.h"

Triangle::Triangle(const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr)
    : Geometry(RTC_GEOMETRY_TYPE_TRIANGLE, devicePtr->handle()) {

    m_vertexBuffer = static_cast<float*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), vertices.size()));
    FillVertexBuffer(vertices);

    const unsigned triCount = vertices.size() / 3;
    m_indexBuffer = static_cast<unsigned*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        3 * sizeof(unsigned), triCount));
    FillIndexBuffer(vertices);
}

void Triangle::FillVertexBuffer(const std::vector<Vertex>& vertices) {
    if (!m_vertexBuffer) return;
    for (size_t i = 0; i < vertices.size(); ++i) {
        m_vertexBuffer[i * 3 + 0] = vertices[i].x;
        m_vertexBuffer[i * 3 + 1] = vertices[i].y;
        m_vertexBuffer[i * 3 + 2] = vertices[i].z;
    }
}

void Triangle::FillIndexBuffer(const std::vector<Vertex>& vertices) {
    const unsigned triCount = static_cast<unsigned>(vertices.size() / 3);
    for (unsigned t = 0; t < triCount; ++t) {
        m_indexBuffer[t * 3 + 0] = t * 3 + 0;
        m_indexBuffer[t * 3 + 1] = t * 3 + 1;
        m_indexBuffer[t * 3 + 2] = t * 3 + 2;
    }
}
