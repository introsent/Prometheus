//
// Created by minaj on 10/4/2025.
//

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "structs/vertex.h"
#include "geometries/geometry.h"

enum class CullingMode {
    NONE = 0,           // No culling (default)
    BACK_FACE = 1,      // Cull back faces
    FRONT_FACE = 2      // Cull front faces
};

class Triangle : public Geometry {
public:
    Triangle(const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr,
             CullingMode culling = CullingMode::NONE);

    [[nodiscard]] CullingMode getCullingMode() const;

private:
    void fillVertexBuffer(const std::vector<Vertex>& vertices) const;
    void fillIndexBuffer(const std::vector<Vertex>& vertices) const;

    float* m_vertexBuffer = nullptr;
    unsigned* m_indexBuffer = nullptr;
    CullingMode m_cullingMode = CullingMode::NONE;
};



#endif //TRIANGLE_H
