//
// Created by minaj on 10/4/2025.
//

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "structs/vertex.h"
#include "geometries/geometry.h"

class Triangle : public Geometry {
public:
    Triangle(SceneManager* scene, const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr);

    [[nodiscard]] const std::vector<Vertex>& getOriginalVertices() const { return m_originalVertices; }
    void updateAABB() override;

private:
    float* m_vertexBuffer = nullptr;
    unsigned* m_indexBuffer = nullptr;
    std::vector<Vertex> m_originalVertices;

    void fillVertexBuffer(const std::vector<Vertex>& vertices) const;
    void fillIndexBuffer(const std::vector<Vertex>& vertices) const;
    void applyTransforms() override;
};

#endif // TRIANGLE_H