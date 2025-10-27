//
// Created by minaj on 10/27/2025.
//

#ifndef SCENE_MANAGER_H
#define SCENE_MANAGER_H
#include <memory>
#include <vector>
#include "embree_scene.h"
#include "material.h"
#include "embree4/rtcore_device.h"

class EmbreeScene;
class Geometry;
class EmbreeDevice;

class SceneManager {
public:
    SceneManager();
    ~SceneManager();

    unsigned char addMaterial(Material* material);
    void addSphere(const glm::vec3& center, float radius, unsigned char materialId);
    void addPlane(const glm::vec3& origin, const glm::vec3& normal, unsigned char materialId);

    void commit();

    RTCScene getScene() const { return m_scenePtr->handle(); }
    const Material* getMaterial(unsigned char materialId) const;
    unsigned char getGeometryMaterial(unsigned int geomID) const;

private:
    std::unique_ptr<EmbreeDevice> m_devicePtr;
    std::unique_ptr<EmbreeScene> m_scenePtr;
    std::vector<std::unique_ptr<Material>> m_materials;
    std::vector<std::unique_ptr<Geometry>> m_geometries;
    std::vector<unsigned char> m_geometryMaterials;
};

#endif //SCENE_MANAGER_H
