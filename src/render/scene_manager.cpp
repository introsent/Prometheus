//
// Optimized scene manager - reduce unnecessary updates
//

#include "scene_manager.h"
#include "global_indices.h"
#include "plane.h"
#include "sphere.h"
#include "triangle.h"
#include "mesh.h"

SceneManager::SceneManager() {
    m_devicePtr = std::make_unique<EmbreeDevice>();
    m_scenePtr = std::make_unique<EmbreeScene>(m_devicePtr.get());
    m_materials.push_back(std::make_unique<Material_SolidColor>(colors::red));
    m_needsCommit = false;
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
    auto sphere = std::make_unique<Sphere>(this, center, radius, m_devicePtr.get());
    sphere->commit();

    unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), sphere->getGeometry());
    m_geometryTypes[geomID] = sphere->getType();
    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(sphere));
    m_needsCommit = true;
}

void SceneManager::addPlane(const glm::vec3& origin, const glm::vec3& normal, unsigned char materialId) {
    auto plane = std::make_unique<Plane>(this, origin, normal, m_devicePtr.get());
    plane->commit();

    unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), plane->getGeometry());
    m_geometryTypes[geomID] = plane->getType();
    m_geometryMaterials.push_back(materialId);
    m_geometries.push_back(std::move(plane));
    m_needsCommit = true;
}

unsigned int SceneManager::addTriangle(const std::vector<Vertex>& vertices, unsigned char materialId) {
    auto triangle = std::make_unique<Triangle>(this, vertices, m_devicePtr.get());
    triangle->commit();

    const unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), triangle->getGeometry());

    m_geometryTypes[geomID] = triangle->getType();
    m_geometryMaterials.push_back(materialId);

    const size_t geometryIndex = m_geometries.size();
    m_triangleIndices.push_back(geometryIndex);
    m_geometries.push_back(std::move(triangle));
    m_needsCommit = true;

    return static_cast<unsigned int>(m_triangleIndices.size() - 1);
}

unsigned int SceneManager::addMesh(const std::vector<Vertex>& vertices,
                                    const std::vector<uint32_t>& indices,
                                    unsigned char materialId) {
    auto mesh = std::make_unique<Mesh>(this, vertices, indices, m_devicePtr.get());
    mesh->commit();

    const unsigned geomID = rtcAttachGeometry(m_scenePtr->handle(), mesh->getGeometry());

    m_geometryTypes[geomID] = mesh->getType();
    m_geometryMaterials.push_back(materialId);

    const size_t geometryIndex = m_geometries.size();
    m_meshIndices.push_back(geometryIndex);
    m_geometries.push_back(std::move(mesh));
    m_needsCommit = true;

    return static_cast<unsigned int>(m_meshIndices.size() - 1);
}

void SceneManager::addLight(Light* light) {
    m_lights.push_back(std::unique_ptr<Light>(light));
}

void SceneManager::commit() {
    if (m_needsCommit) {
        m_scenePtr->commit();
        m_needsCommit = false;
    }
}

void SceneManager::update(const Timer* pTimer) {
    const float dt = pTimer->getElapsed();
    const float total = pTimer->getTotal();

    for (auto& geometry : m_geometries)
    {
        geometry->update(dt, total);
    }

    if (m_needsCommit)
        commit();
}

const Material* SceneManager::getMaterial(unsigned char materialId) const {
    if (materialId < m_materials.size()) {
        return m_materials[materialId].get();
    }
    return m_materials[0].get();
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
    return 0;
}

RTCGeometryType SceneManager::getGeometryType(unsigned int geomID) const {
    if (const auto it = m_geometryTypes.find(geomID); it != m_geometryTypes.end()) {
        return it->second;
    }
    return RTC_GEOMETRY_TYPE_TRIANGLE;
}

Triangle* SceneManager::getTriangle(unsigned int triangleIndex) const {
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

Mesh* SceneManager::getMesh(unsigned int meshIndex) const {
    if (meshIndex < m_meshIndices.size()) {
        if (const size_t geometryIndex = m_meshIndices[meshIndex]; geometryIndex < m_geometries.size()) {
            return dynamic_cast<Mesh*>(m_geometries[geometryIndex].get());
        }
    }
    return nullptr;
}

size_t SceneManager::getMeshCount() const {
    return m_meshIndices.size();
}