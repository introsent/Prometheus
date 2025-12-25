//
// Created by minaj on 10/27/2025.
//

#ifndef SCENE_MANAGER_H
#define SCENE_MANAGER_H
#include <memory>
#include <unordered_map>
#include <vector>

#include "area_light.h"
#include "embree_scene.h"
#include "light.h"
#include "material.h"
#include "timer.h"
#include "triangle.h"
#include "mesh.h"
#include "vertex.h"
#include "embree4/rtcore_device.h"

enum class CullingMode;
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
    unsigned addTriangle(const std::vector<Vertex>& vertices, unsigned char materialId);
    unsigned addMesh(const std::vector<Vertex>& vertices, const std::vector<uint32_t>& indices, unsigned char materialId);
    void addLight(Light* light);

    unsigned int addTriangleAreaLight(const std::vector<Vertex>& vertices,
                                 const glm::vec3& emission,
                                 float intensity,
                                 unsigned char materialId = 0);

    unsigned int addMeshAreaLight(const std::vector<Vertex>& vertices,
                                 const std::vector<uint32_t>& indices,
                                 const glm::vec3& albedo,
                                 const glm::vec3& emission,
                                 float intensity,
                                 unsigned char materialId = 0);

    void setAreaLightSamplingStrategy(SamplingStrategy strategy) const;

    void commit();
    void markDirty() { m_needsCommit = true; }

    void update(const Timer* pTimer);

    [[nodiscard]] RTCScene getScene() const { return m_scenePtr->handle(); }
    [[nodiscard]] const Material* getMaterial(unsigned char materialId) const;
    [[nodiscard]] const Light* getLight(unsigned char lightId) const;
    [[nodiscard]] unsigned char getGeometryMaterial(unsigned int geomID) const;
    [[nodiscard]] RTCGeometryType getGeometryType(unsigned int geomID) const;
    [[nodiscard]] const std::vector<std::unique_ptr<Light>>& getLights() const noexcept { return m_lights; }

    // Triangle access for transformations
    [[nodiscard]] Triangle* getTriangle(unsigned int triangleIndex) const;
    [[nodiscard]] size_t getTriangleCount() const;

    // Mesh access for transformations
    [[nodiscard]] Mesh* getMesh(unsigned int meshIndex) const;
    [[nodiscard]] size_t getMeshCount() const;

    // Get area lights
    [[nodiscard]] const std::vector<std::unique_ptr<TriangleAreaLight>>& getTriangleAreaLights() const {
        return m_triangleAreaLights;
    }
    [[nodiscard]] const std::vector<std::unique_ptr<MeshAreaLight>>& getMeshAreaLights() const {
        return m_meshAreaLights;
    }


private:
    std::unique_ptr<EmbreeDevice> m_devicePtr;
    std::unique_ptr<EmbreeScene> m_scenePtr;
    std::vector<std::unique_ptr<Material>> m_materials;
    std::vector<std::unique_ptr<Geometry>> m_geometries;
    std::vector<std::unique_ptr<Light>> m_lights;
    std::vector<unsigned char> m_geometryMaterials;
    std::unordered_map<unsigned, RTCGeometryType> m_geometryTypes;

    std::vector<std::unique_ptr<TriangleAreaLight>> m_triangleAreaLights;
    std::vector<std::unique_ptr<MeshAreaLight>> m_meshAreaLights;

    std::vector<size_t> m_triangleIndices;
    std::vector<size_t> m_meshIndices;

    bool m_needsCommit;
};

#endif //SCENE_MANAGER_H