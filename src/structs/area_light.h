//
// Created by ivans on 20/12/2025.
//

#ifndef PROMETHEUS_AREA_LIGHT_H
#define PROMETHEUS_AREA_LIGHT_H

#include "glm/vec3.hpp"
#include "glm/geometric.hpp"
#include <vector>
#include <memory>
#include <random>

class Triangle;
class SceneManager;

/// Sampling strategy
enum class SamplingStrategy
{
    Uniform,
    AreaImportance,
    HierarchicalFlux,
    VisibilityAwareHierarchical
};

/// Sample result from area light
struct AreaLightSample
{
    glm::vec3 position;
    glm::vec3 normal;
    float pdf; // probability density function (the larger the light, the brighter should it be)
    glm::vec3 radiance;
    float area;
};


/// Base class for area light sampling strategies
class AreaLightSampler
{
public:
    virtual ~AreaLightSampler() = default;

    // sample a point of the area light
    [[nodiscard]] virtual AreaLightSample sample(const glm::vec3& shadingPoint,
                                  float u1, float u2) const = 0;

    // calculate pdf for the point
    [[nodiscard]] virtual float pdf(const glm::vec3& shadingPoint,
                    const glm::vec3& lightPoint) const = 0;

    // get total flux (for importance sampling)
    [[nodiscard]] virtual float getTotalFlux() const = 0;
};


/// S1: Uniform sampling over triangle surface
class UniformTriangleSampler : public AreaLightSampler {
public:
    UniformTriangleSampler(const glm::vec3& v0,
                          const glm::vec3& v1,
                          const glm::vec3& v2,
                          const glm::vec3& normal,
                          float area,
                          const glm::vec3& emission,
                          float intensity);

    [[nodiscard]] AreaLightSample sample(const glm::vec3& shadingPoint,
                          float u1, float u2) const override;

    [[nodiscard]] float pdf(const glm::vec3& shadingPoint,
             const glm::vec3& lightPoint) const override;

    [[nodiscard]] float getTotalFlux() const override;

private:
    glm::vec3 m_v0, m_v1, m_v2;
    glm::vec3 m_normal;
    float m_area;
    glm::vec3 m_emission;
    float m_intensity;
};

// tiny helpers

static float safeAcos(float x) { return std::acos(std::clamp(x, -1.0f, 1.0f)); }
static float safeSqrt(float x) { return std::sqrt(std::max(0.0f, x)); }

// Gram-Schmidt: return (v - proj_a(v)) normalized; if near zero, return any perpendicular to a
static glm::vec3 gramSchmidtNormalize(const glm::vec3 &v, const glm::vec3 &a) {
    const glm::vec3 proj = a * glm::dot(v, a);
    glm::vec3 w = v - proj;
    if (float len2 = glm::dot(w, w); len2 <= 1e-10f) {
        // fallback: find any vector perpendicular to a
        const glm::vec3 absA = glm::abs(a);
        const glm::vec3 other = (absA.x <= absA.y && absA.x <= absA.z) ? glm::vec3(1,0,0)
                         : (absA.y <= absA.x && absA.y <= absA.z) ? glm::vec3(0,1,0)
                                                                 : glm::vec3(0,0,1);
        w = glm::cross(a, other);
        len2 = glm::dot(w, w);
        if (len2 <= 1e-10f) return glm::normalize(glm::cross(a, glm::vec3(1,0,0) + 0.0001f)); // last resort
    }
    return glm::normalize(w);
}

class AreaImportanceTriangleSampler : public AreaLightSampler {
public:
    AreaImportanceTriangleSampler(const glm::vec3& v0,
                         const glm::vec3& v1,
                         const glm::vec3& v2,
                         const glm::vec3& normal,
                         float area,
                         const glm::vec3& emission,
                         float intensity);


    [[nodiscard]] AreaLightSample sample(const glm::vec3& shadingPoint,
                                         float u1, float u2) const override;

    [[nodiscard]] float pdf(const glm::vec3& shadingPoint,
             const glm::vec3& lightPoint) const override;

    [[nodiscard]] float getTotalFlux() const override;
    [[nodiscard]] float calculateSolidAngle(const glm::vec3& p) const;
private:

    glm::vec3 sampleSphericalTriangle(const glm::vec3& p, float u1, float u2, float* outPdf = nullptr) const;
    [[nodiscard]] std::pair<float, glm::vec3> getRayPlaneIntersection(const glm::vec3& shadingPoint, const glm::vec3& direction) const;
    void recomputePositionUsingBarycentricCoordinates(float t, glm::vec3& position) const;


    glm::vec3 m_v0, m_v1, m_v2;
    glm::vec3 m_normal;
    float m_area;
    glm::vec3 m_emission;
    float m_intensity;
};

class HierarchicalFluxTriangleSampler : public AreaLightSampler {
public:
    HierarchicalFluxTriangleSampler(const glm::vec3& v0,
                         const glm::vec3& v1,
                         const glm::vec3& v2,
                         const glm::vec3& normal,
                         float area,
                         const glm::vec3& emission,
                         float intensity);


    [[nodiscard]] AreaLightSample sample(const glm::vec3& shadingPoint,
                                         float u1, float u2) const override;

    [[nodiscard]] float pdf(const glm::vec3& shadingPoint,
             const glm::vec3& lightPoint) const override;

    [[nodiscard]] float getTotalFlux() const override;

private:
    glm::vec3 m_v0, m_v1, m_v2;
    glm::vec3 m_normal;
    float m_area;
    glm::vec3 m_emission;
    float m_intensity;
};

/// Triangle area light
class TriangleAreaLight
{
public:
    TriangleAreaLight(unsigned int triangleIndex,
                     const glm::vec3& emission,
                     float intensity,
                     SceneManager* scene);
    ~TriangleAreaLight() = default;

    [[nodiscard]] AreaLightSample sample(const glm::vec3& shadingPoint,
                         float u1, float u2) const;

    [[nodiscard]] float pdf(const glm::vec3& shadingPoint,
             const glm::vec3& lightPoint) const;

    [[nodiscard]] glm::vec3 getNormal() const;
    [[nodiscard]] float getArea() const { return m_area; }
    [[nodiscard]] glm::vec3 getEmission() const { return m_emission; }
    [[nodiscard]] float getIntensity() const { return m_intensity; }
    [[nodiscard]] unsigned int getTriangleIndex() const { return m_triangleIndex; }

    // Strategy management
    void setSamplingStrategy(SamplingStrategy strategy);
    [[nodiscard]] SamplingStrategy getSamplingStrategy() const { return m_strategy; }

    // Power calculation (emission * area * intensity)
    [[nodiscard]] float getPower() const;

    // For mesh sampling
    [[nodiscard]] float calculateSolidAngle(const glm::vec3& p) const;

private:
    void initializeSampler();
    void cacheGeometry();

    unsigned int m_triangleIndex;
    SceneManager* m_scene;

    glm::vec3 m_emission;   // emissive color
    float m_intensity;      // emission intensity

    // cached geometry
    glm::vec3 m_v0{}, m_v1{}, m_v2{};   // triangle vertices
    glm::vec3 m_normal{};           // triangle normal
    float m_area{};                 // triangle area

    // sampling
    SamplingStrategy m_strategy;
    std::unique_ptr<AreaLightSampler> m_sampler;
};

/// Mesh area light (for now collection of triangles)
class MeshAreaLight {
public:
    MeshAreaLight(unsigned int meshIndex,
                 const glm::vec3& emission,
                 float intensity,
                 SceneManager* scene);
    ~MeshAreaLight() = default;

    // sampling methods with explicit strategy
    [[nodiscard]] AreaLightSample sampleUniform(const glm::vec3& shadingPoint,
                          float u1, float u2, float u3) const;

    [[nodiscard]] AreaLightSample sampleAreaImportance(const glm::vec3& shadingPoint,
                          float u1, float u2, float u3) const;

    // convenience method that uses current strategy
    [[nodiscard]] AreaLightSample sample(const glm::vec3& shadingPoint,
                          float u1, float u2, float u3) const;

    // PDF methods
    [[nodiscard]] float pdfUniform(const glm::vec3& shadingPoint,
                     const glm::vec3& lightPoint) const;

    [[nodiscard]] float pdfAreaImportance(const glm::vec3& shadingPoint,
                         const glm::vec3& lightPoint) const;

    [[nodiscard]] float pdf(const glm::vec3& shadingPoint,
             const glm::vec3& lightPoint) const;

    [[nodiscard]] float getTotalPower() const;

    void setSamplingStrategy(SamplingStrategy strategy);
    [[nodiscard]] SamplingStrategy getSamplingStrategy() const { return m_strategy; }

private:
    struct TriangleData {
        glm::vec3 v0, v1, v2;
        glm::vec3 normal;
        float area;
        std::unique_ptr<UniformTriangleSampler> uniformSampler;
        std::unique_ptr<AreaImportanceTriangleSampler> areaImportanceSampler;

        TriangleData(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
                    const glm::vec3& normal, float area,
                    const glm::vec3& emission, float intensity)
            : v0(v0), v1(v1), v2(v2), normal(normal), area(area)
        {
            uniformSampler = std::make_unique<UniformTriangleSampler>(
                v0, v1, v2, normal, area, emission, intensity);
            areaImportanceSampler = std::make_unique<AreaImportanceTriangleSampler>(
                v0, v1, v2, normal, area, emission, intensity);
        }
    };

    // Check if point is inside triangle
    [[nodiscard]] bool isPointInTriangle(const glm::vec3& p,
                                        const glm::vec3& v0,
                                        const glm::vec3& v1,
                                        const glm::vec3& v2,
                                        float& u, float& v, float& w) const;

    unsigned int m_meshIndex;
    SceneManager* m_scene;
    glm::vec3 m_emission;
    float m_intensity;
    SamplingStrategy m_strategy;

    std::vector<TriangleData> m_triangles;
    std::vector<float> m_areaCDF;      // For uniform sampling
    float m_totalArea;
};


struct BVHNode {
    BVHNode();

    // for nodes
    BVHNode(int left, int right, float flux, float area, bool leaf);

    void initializeInternalNode(int left, int right, float flux, float a);

    void initializeLeafNode(int start, int end, float flux, float a);

    // for internal nodes
    int leftChild{};
    int rightChild{};

    // for leaf nodes
    int startTri{};
    int endTri{};

    // for bounding box
    glm::vec3 bboxMin{};
    glm::vec3 bboxMax{};

    // total flux in current BVH node
    float totalFlux{};

    // total area in this node
    float totalArea{};

    // whether this node is a leaf or no
    bool isLeaf;

    // centroid (for sorting)
    [[nodiscard]] glm::vec3 getCentroid() const { return (bboxMin + bboxMax) * 0.5f; }
};


#endif //PROMETHEUS_AREA_LIGHT_H