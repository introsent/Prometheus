//
// Created by minaj on 10/4/2025.
//

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <functional>
#include <vector>
#include <glm/vec3.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "init/embree_scene.h"
#include "embree4/rtcore.h"

struct Vertex;
using UpdateFunc = std::function<void(float dt, float totalTime)>;
class SceneManager;
class Geometry {
public:
    explicit Geometry(SceneManager* scene, RTCGeometryType geometryType, RTCDevice device);

    [[nodiscard]] RTCGeometryType getType() const;
    [[nodiscard]] RTCGeometry getGeometry() const;

    virtual ~Geometry() = default;

    void commit() const;
    unsigned attach(const EmbreeScene* embreeScenePtr) const;
    void release() const;

    // Transform operations (now available for all geometries)
    void translate(const glm::vec3& translation);
    void rotateY(float yawRadians);
    void rotate(float angle, const glm::vec3& axis); // More general rotation
    void scale(const glm::vec3& scaleVec);
    virtual void updateTransforms(); // Virtual so subclasses can override

    // AABB operations
    virtual void updateAABB() = 0; // Pure virtual - each geometry computes differently
    void updateTransformedAABB();
    [[nodiscard]] glm::vec3 getTransformedMinAABB() const { return m_transformedMinAABB; }
    [[nodiscard]] glm::vec3 getTransformedMaxAABB() const { return m_transformedMaxAABB; }
    [[nodiscard]] glm::vec3 getMinAABB() const { return m_minAABB; }
    [[nodiscard]] glm::vec3 getMaxAABB() const { return m_maxAABB; }

    // Update Function
    void setUpdateFunc(UpdateFunc func) { m_updateFunc = std::move(func); }
    void update(float dt, float totalTime);

protected:
    RTCGeometryType m_type;
    RTCGeometry m_geometry = nullptr;
    SceneManager* m_scene = nullptr;

    // Transform matrices
    glm::mat4 m_translationMatrix{1.0f};
    glm::mat4 m_rotationMatrix{1.0f};
    glm::mat4 m_scaleMatrix{1.0f};

    // AABB bounds
    glm::vec3 m_minAABB{0.f};
    glm::vec3 m_maxAABB{0.f};
    glm::vec3 m_transformedMinAABB{0.f};
    glm::vec3 m_transformedMaxAABB{0.f};

    // Helper methods
    virtual void applyTransforms() = 0; // Pure virtual - each geometry applies differently
    [[nodiscard]] static glm::vec3 transformPoint(const glm::vec3& point, const glm::mat4& matrix);
    [[nodiscard]] glm::mat4 getFinalTransform() const;

    // Update Function
    UpdateFunc m_updateFunc;
};

#endif //GEOMETRY_H