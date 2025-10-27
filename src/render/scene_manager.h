//
// Created by minaj on 10/27/2025.
//

#ifndef SCENE_MANAGER_H
#define SCENE_MANAGER_H
#include <memory>
#include <vector>
#include "embree_scene.h"
#include "light.h"
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
    void addLight(Light* light);

    void commit();

    [[nodiscard]] RTCScene getScene() const { return m_scenePtr->handle(); }
    [[nodiscard]] const Material* getMaterial(unsigned char materialId) const;
    [[nodiscard]] const Light* getLight(unsigned char lightId) const;
    [[nodiscard]] unsigned char getGeometryMaterial(unsigned int geomID) const;

    [[nodiscard]] const std::vector<std::unique_ptr<Light>>& getLights() const noexcept { return m_lights; }
private:
    std::unique_ptr<EmbreeDevice> m_devicePtr;
    std::unique_ptr<EmbreeScene> m_scenePtr;
    std::vector<std::unique_ptr<Material>> m_materials;
    std::vector<std::unique_ptr<Geometry>> m_geometries;
    std::vector<std::unique_ptr<Light>> m_lights;
    std::vector<unsigned char> m_geometryMaterials;
};

#endif //SCENE_MANAGER_H
