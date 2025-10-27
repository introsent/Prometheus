//
// Created by minaj on 10/4/2025.
//

#include "triangle.h"

Triangle::Triangle(const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr) : Geometry(RTC_GEOMETRY_TYPE_TRIANGLE, devicePtr->handle()) {

    m_vertexBuffer = static_cast<float*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        vertices.size() * sizeof(float), vertices.size()));
    FillVertexBuffer(vertices);

    m_indexBuffer = static_cast<unsigned*>(rtcSetNewGeometryBuffer(m_geometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        vertices.size() * sizeof(unsigned), vertices.size()));
    FillIndexBuffer(vertices);
}

void Triangle::FillVertexBuffer(const std::vector<Vertex> &vertices) {
    if (m_vertexBuffer == nullptr) return;

    for (int vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
    {
        m_vertexBuffer[vertexIndex * vertices.size()] = vertices[vertexIndex].x;
        m_vertexBuffer[vertexIndex * vertices.size() + 1] = vertices[vertexIndex].y;
        m_vertexBuffer[vertexIndex * vertices.size() + 2] = vertices[vertexIndex].z;
    }
}

void Triangle::FillIndexBuffer(const std::vector<Vertex> &vertices) {
    for (int i = 0; i < vertices.size(); ++i)
    {
        m_indexBuffer[i] = i;
    }
}
