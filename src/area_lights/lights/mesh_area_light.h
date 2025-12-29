//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_MESH_AREA_LIGHT_H
#define PROMETHEUS_MESH_AREA_LIGHT_H
#include "sampler_base.h"
#include <glm/glm.hpp>
#include <vector>
#include <memory>
#include "acceleration/bvh_node.h"
#include "samplers/area_importance_triangle_sampler.h"
#include "samplers/uniform_triangle_sampler.h"
#include "visibility/spatial_visibility_cache.h"

class SceneManager;
class VisibilityAwareHierarchicalSampler;

// ============================================================================
// mesh area light
// manages collection of emitting triangles with multiple sampling strategies
// ============================================================================
class MeshAreaLight {
public:
    MeshAreaLight(
        unsigned int meshIndex,
        const glm::vec3& emission,
        float intensity,
        SceneManager* scene);

    ~MeshAreaLight();

    // main sampling interface (dispatches to current strategy)
    [[nodiscard]] AreaLightSample sample(
        const glm::vec3& shadingPoint,
        float u1, float u2, float u3) const;
    [[nodiscard]] AreaLightSample sample(const glm::vec3& shadingPoint,
                          const glm::vec3& surfaceNormal,
                          float u1, float u2, float u3) const;


    // pdf evaluation (dispatches to current strategy)
    [[nodiscard]] float pdf(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    // strategy-specific methods
    [[nodiscard]] AreaLightSample sampleUniform(
        const glm::vec3& shadingPoint,
        float u1, float u2, float u3) const;

    [[nodiscard]] AreaLightSample sampleAreaImportance(
        const glm::vec3& shadingPoint,
        float u1, float u2, float u3) const;

    [[nodiscard]] AreaLightSample sampleHierarchicalFlux(
        const glm::vec3& shadingPoint,
        float u1, float u2, float u3) const;

    [[nodiscard]] AreaLightSample sampleVisibilityAware(
        const glm::vec3& shadingPoint,
        float u1, float u2, float u3) const;
    [[nodiscard]] AreaLightSample sampleVisibilityAware(const glm::vec3& shadingPoint,
                                         const glm::vec3& surfaceNormal,
                                         float u1, float u2, float u3) const;


    // pdf methods
    [[nodiscard]] float pdfUniform(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    [[nodiscard]] float pdfAreaImportance(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    [[nodiscard]] float pdfHierarchicalFlux(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    // configuration
    void setSamplingStrategy(SamplingStrategy strategy);
    void setTriangleIntensity(size_t triangleIndex, float intensity);
    void updateBVH();

    // accessors
    [[nodiscard]] float getTotalPower() const;
    [[nodiscard]] SamplingStrategy getSamplingStrategy() const { return m_strategy; }
    [[nodiscard]] size_t getTriangleCount() const { return m_triangles.size(); }
    [[nodiscard]] float getTotalArea() const { return m_totalArea; }

    int getLastSampledNode() const { return m_lastSampledNode; }
    void setLastSampledNode(int leafNodeIndex) { m_lastSampledNode = leafNodeIndex; }


    // static visibility cache (shared across all mesh lights)
    static void initializeVisibilityCache(
        const glm::vec3& sceneMin,
        const glm::vec3& sceneMax,
        int resolution = 16);

    static void getVisibilityStats(
        int& totalCells,
        int& activeCells,
        int& totalSamples);

    static void clearVisibilityCache();
    VisibilityAwareHierarchicalSampler* getVisibilitySampler() const {
        return m_visAwareSampler.get();
    };

    [[nodiscard]] bool containsPoint(const glm::vec3& point, float epsilon = 1e-4f) const;
    [[nodiscard]] glm::vec3 getEmissionAt(const glm::vec3& point) const;

private:
    friend class VisibilityAwareHierarchicalSampler;

    // triangle data with pre-computed samplers
    struct TriangleData {
        glm::vec3 v0, v1, v2;
        glm::vec3 normal;
        float area;
        float intensity;
        std::unique_ptr<UniformTriangleSampler> uniformSampler;
        std::unique_ptr<AreaImportanceTriangleSampler> areaImportanceSampler;
    };

    // initialization helpers
    void extractTrianglesFromMesh();
    void buildAreaCDF();

    // BVH management
    void buildBVH();
    void calculateLeafFlux();
    void buildTriangleToNodeMapping();

    // hierarchical sampling helpers
    [[nodiscard]] int selectBVHLeafNode(
        const glm::vec3& shadingPoint,
        float u,
        float& outPathPdf) const;

    [[nodiscard]] size_t selectTriangleFromLeaf(
        int leafNodeIndex,
        float u,
        float& outSelectionProb) const;

    [[nodiscard]] float computeHierarchicalPdfForTriangle(
        size_t triangleIndex,
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    [[nodiscard]] float computePathPdfToLeaf(
        int leafNodeIndex,
        const glm::vec3& shadingPoint) const;

    // geometry queries
    [[nodiscard]] size_t findTriangleContainingPoint(
        const glm::vec3& point) const;

    static bool isPointInTriangle(
        const glm::vec3& point,
        const glm::vec3& v0,
        const glm::vec3& v1,
        const glm::vec3& v2,
        float& u, float& v, float& w) ;

    // member data
    unsigned int m_meshIndex;
    SceneManager* m_scene;
    glm::vec3 m_emission;
    float m_intensity;
    SamplingStrategy m_strategy;

    // triangle storage
    std::vector<TriangleData> m_triangles;
    float m_totalArea;

    // uniform sampling
    std::vector<float> m_areaCDF;

    // hierarchical sampling
    std::vector<BVHNode> m_bvhNodes;
    int m_rootNodeIndex;
    std::vector<int> m_triangleToLeafNode;

    // visibility-aware sampling
    static std::unique_ptr<SpatialVisibilityCache> s_visibilityCache;
    std::unique_ptr<VisibilityAwareHierarchicalSampler> m_visAwareSampler;
    mutable int m_lastSampledNode;
};
#endif //PROMETHEUS_MESH_AREA_LIGHT_H