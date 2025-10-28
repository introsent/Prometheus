//
// Created by minaj on 10/27/2025.
//

#include "scene_manager.h"

#include "global_indices.h"
#include "plane.h"
#include "sphere.h"
#include "triangle.h"

SceneManager::SceneManager() {
    m_devicePtr = std::make_unique<EmbreeDevice>();
    m_scenePtr = std::make_unique<EmbreeScene>(m_devicePtr.get());

    // Default red material at index 0
    m_materials.push_back(std::make_unique<Material_SolidColor>(colors::red));
}

SceneManager::~SceneManager() {
    m_scenePtr->release();
    m_devicePtr->release();
}

unsigned char SceneManager::addMaterial(Material* material) {
    m_materials.push_back(std::unique_ptr<Material>(material));
    return static_cast<unsigned char>(m_materials.size() - 1);
}

void SceneManager::addSphere(const glm::vec3& center, float radius, unsigned char materialId) {
    auto sphere = std::make_unique<Sphere>(center, radius, m_devicePtr.get());
    sphere->commit();

    // Get the geomID when attaching
    unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), sphere->getGeometry());

    // Store the geometry type
    m_geometryTypes[geomID] = sphere->getType();

    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(sphere));
}

void SceneManager::addPlane(const glm::vec3& origin, const glm::vec3& normal, unsigned char materialId) {
    auto plane = std::make_unique<Plane>(origin, normal, m_devicePtr.get());
    plane->commit();

    // Get the geomID when attaching
    unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), plane->getGeometry());

    // Store the geometry type
    m_geometryTypes[geomID] = plane->getType();

    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(plane));
}

unsigned int SceneManager::addTriangle(const std::vector<Vertex>& vertices, unsigned char materialId,
                                       CullingMode culling) {
    auto triangle = std::make_unique<Triangle>(vertices, m_devicePtr.get(), culling);
    triangle->commit();

    const unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), triangle->getGeometry());

    m_geometryTypes[geomID] = triangle->getType();
    m_cullingModes[geomID] = culling;
    m_geometryMaterials.push_back(materialId);

    // Store the index in the geometries vector for later retrieval
    const size_t geometryIndex = m_geometries.size();
    m_triangleIndices.push_back(geometryIndex);

    m_geometries.push_back(std::move(triangle));

    // Return the triangle index (not the geometry index)
    return static_cast<unsigned int>(m_triangleIndices.size() - 1);
}

void SceneManager::addLight(Light* light) {
    m_lights.push_back(std::unique_ptr<Light>(light));
}

void SceneManager::commit() {
    m_scenePtr->commit();
}

void SceneManager::update(const Timer* pTimer) {
    // Calculate rotation angle (oscillates between 0 and 2*PI)
    const float yawAngle = (std::cos(pTimer->getTotal()) + 1.f) / 2.f  * 6.283185307179586476925f;

    // Rotate all triangles
    for (const unsigned int triIdx : g_triangleIndices) {
        if (Triangle* pTriangle = getTriangle(triIdx)) {
            pTriangle->rotateY(yawAngle);
            pTriangle->updateAABB();
            pTriangle->updateTransforms();
        }
    }

    // Recommit the scene after transformations
    commit();
}

const Material* SceneManager::getMaterial(unsigned char materialId) const {
    if (materialId < m_materials.size()) {
        return m_materials[materialId].get();
    }
    return m_materials[0].get(); // Default to red
}

const Light* SceneManager::getLight(unsigned char lightId) const {
    if (lightId < m_lights.size()) {
        return m_lights[lightId].get();
    }
    return m_lights[0].get();
}

unsigned char SceneManager::getGeometryMaterial(unsigned int geomID) const {
    if (geomID < m_geometryMaterials.size()) {
        return m_geometryMaterials[geomID];
    }
    return 0; // Default material
}

RTCGeometryType SceneManager::getGeometryType(unsigned int geomID) const {
    if (const auto it = m_geometryTypes.find(geomID); it != m_geometryTypes.end()) {
        return it->second;
    }
    return RTC_GEOMETRY_TYPE_TRIANGLE; // Default fallback
}

CullingMode SceneManager::getCullingMode(unsigned geomID) const {
    const auto it = m_cullingModes.find(geomID);
    return (it != m_cullingModes.end()) ? it->second : CullingMode::NONE;
}

Triangle *SceneManager::getTriangle(unsigned int triangleIndex) const {
    if (triangleIndex < m_triangleIndices.size()) {
        if (const size_t geometryIndex = m_triangleIndices[triangleIndex]; geometryIndex < m_geometries.size()) {
            return dynamic_cast<Triangle*>(m_geometries[geometryIndex].get());
        }
    }
    return nullptr;
}

size_t SceneManager::getTriangleCount() const {
    return m_triangleIndices.size();
}