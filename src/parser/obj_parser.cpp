//
// Created by ivans on 19/12/2025.
//

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include "obj_parser.h"
#include <iostream>

bool ParseOBJ(
    const std::string& filename,
    std::vector<Vertex>& outVertices,
    std::vector<uint32_t>& outIndices)
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &err, filename.c_str())) {
        std::cerr << err << '\n';
        return false;
    }

    for (const auto& shape : shapes) {
        size_t indexOffset = 0;

        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
            int fv = shape.mesh.num_face_vertices[f];

            for (int v = 0; v < fv; v++) {
                const auto& idx = shape.mesh.indices[indexOffset + v];

                Vertex vert{};
                vert.position = {
                    attrib.vertices[3 * idx.vertex_index + 0],
                    attrib.vertices[3 * idx.vertex_index + 1],
                    attrib.vertices[3 * idx.vertex_index + 2]
                };

                if (!attrib.normals.empty() && idx.normal_index >= 0) {
                    vert.normal = {
                        attrib.normals[3 * idx.normal_index + 0],
                        attrib.normals[3 * idx.normal_index + 1],
                        attrib.normals[3 * idx.normal_index + 2]
                    };
                }

                outVertices.push_back(vert);
                outIndices.push_back(static_cast<uint32_t>(outVertices.size() - 1));
            }

            indexOffset += fv;
        }
    }

    return true;
}


