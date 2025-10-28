//
// Created by minaj on 10/4/2025.
//

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "structs/vertex.h"
#include "geometries/geometry.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

enum class CullingMode {
    NONE = 0,           // No culling (default)
    BACK_FACE = 1,      // Cull back faces
    FRONT_FACE = 2      // Cull front faces
};

class Triangle : public Geometry {
public:
    Triangle(const std::vector<Vertex>& vertices, const EmbreeDevice* devicePtr,
             CullingMode culling = CullingMode::BACK_FACE);

    void fillVertexBuffer(const std::vector<Vertex>& vertices) const;
    void fillIndexBuffer(const std::vector<Vertex>& vertices) const;

    [[nodiscard]] CullingMode getCullingMode() const;

    // Transform operations
    void translate(const glm::vec3& translation);
    void rotateY(float yawRadians);
    void scale(const glm::vec3& scaleVec);
    void updateTransforms();

    // AABB operations
    void updateAABB();
    void updateTransformedAABB();
    [[nodiscard]] glm::vec3 getTransformedMinAABB() const { return m_transformedMinAABB; }
    [[nodiscard]] glm::vec3 getTransformedMaxAABB() const { return m_transformedMaxAABB; }

private:
    float* m_vertexBuffer = nullptr;
    unsigned* m_indexBuffer = nullptr;
    CullingMode m_cullingMode;

    // Original vertices (untransformed)
    std::vector<Vertex> m_originalVertices;

    // Transform matrices
    glm::mat4 m_translationMatrix{1.0f};
    glm::mat4 m_rotationMatrix{1.0f};
    glm::mat4 m_scaleMatrix{1.0f};

    // AABB bounds
    glm::vec3 m_minAABB{0.f};
    glm::vec3 m_maxAABB{0.f};
    glm::vec3 m_transformedMinAABB{0.f};
    glm::vec3 m_transformedMaxAABB{0.f};

    void applyTransforms() const;
    [[nodiscard]] static glm::vec3 transformPoint(const glm::vec3& point, const glm::mat4& matrix) ;
};

#endif // TRIANGLE_H