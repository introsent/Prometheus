//
// Created by minaj on 10/27/2025.
//

#include "scene_manager.h"

#include "plane.h"
#include "sphere.h"

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
    sphere->attach(m_scenePtr.get());

    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(sphere));
}

void SceneManager::addPlane(const glm::vec3& origin, const glm::vec3& normal, unsigned char materialId) {
    auto plane = std::make_unique<Plane>(origin, normal, m_devicePtr.get());
    plane->commit();
    plane->attach(m_scenePtr.get());

    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(plane));
}

void SceneManager::addLight(Light* light)
{
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

const Light * SceneManager::getLight(unsigned char lightId) const {
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

