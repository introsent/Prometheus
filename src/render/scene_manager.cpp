//
// Created by minaj on 10/27/2025.
//

#include "scene_manager.h"
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

void SceneManager::addTriangle(const std::vector<Vertex> &vertices, unsigned char materialId, CullingMode culling) {
    auto triangle = std::make_unique<Triangle>(vertices, m_devicePtr.get(), culling);
    triangle->commit();

    unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), triangle->getGeometry());

    m_geometryTypes[geomID] = triangle->getType();
    m_cullingModes[geomID] = culling;  // Store culling mode
    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(triangle));
}


void SceneManager::addLight(Light* light) {
    m_lights.push_back(std::unique_ptr<Light>(light));
}

void SceneManager::commit() {
    m_scenePtr->commit();
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
