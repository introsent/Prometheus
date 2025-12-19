//
// Created by ivans on 19/12/2025.
//

#ifndef PROMETHEUS_MESH_H
#define PROMETHEUS_MESH_H

#include "geometries/geometry.h"
#include "structs/vertex.h"
#include <vector>

class Mesh final : public Geometry {
public:
    Mesh(SceneManager* scene, const std::vector<Vertex>& vertices,
         const std::vector<uint32_t>& indices,
         const EmbreeDevice* device);

    void updateAABB() override;

private:
    std::vector<Vertex> m_originalVertices;
    std::vector<uint32_t> m_indices;
    float* m_vertexBuffer = nullptr;

    void uploadToEmbree(const std::vector<Vertex>& vertices,
                        const std::vector<uint32_t>& indices);
    void applyTransforms() override;
};

#endif //PROMETHEUS_MESH_H