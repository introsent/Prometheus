//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_VISIBILITY_AWARE_MIS_H
#define PROMETHEUS_VISIBILITY_AWARE_MIS_H

#include <glm/glm.hpp>
#include <vector>
#include <memory>
#include <atomic>
#include <mutex>
#include <glm/ext/scalar_constants.hpp>
class MeshAreaLight;
class SceneManager;
class Ray;

/// Visibility Cache Cell
struct VisibilityCacheCell {
    struct NodeVisibility {
        int nodeIndex = -1;
        std::atomic<uint32_t> visibleSamples{0};
        std::atomic<uint32_t> totalSamples{0};

        NodeVisibility() = default;

        // Delete copy (atomics can't be copied)
        NodeVisibility(const NodeVisibility&) = delete;
        NodeVisibility& operator=(const NodeVisibility&) = delete;

        // Implement move by loading/storing atomic values
        NodeVisibility(NodeVisibility&& other) noexcept
            : nodeIndex(other.nodeIndex)
            , visibleSamples(other.visibleSamples.load(std::memory_order_relaxed))
            , totalSamples(other.totalSamples.load(std::memory_order_relaxed))
        {}

        NodeVisibility& operator=(NodeVisibility&& other) noexcept {
            if (this != &other) {
                nodeIndex = other.nodeIndex;
                visibleSamples.store(other.visibleSamples.load(std::memory_order_relaxed),
                                    std::memory_order_relaxed);
                totalSamples.store(other.totalSamples.load(std::memory_order_relaxed),
                                  std::memory_order_relaxed);
            }
            return *this;
        }

        float getVisibilityProbability() const {
            uint32_t total = totalSamples.load(std::memory_order_relaxed);
            if (total == 0) return 1.0f;
            uint32_t visible = visibleSamples.load(std::memory_order_relaxed);
            return static_cast<float>(visible) / static_cast<float>(total);
        }

        void recordSample(bool wasVisible) {
            totalSamples.fetch_add(1, std::memory_order_relaxed);
            if (wasVisible) {
                visibleSamples.fetch_add(1, std::memory_order_relaxed);
            }
        }
    };

    glm::vec3 center;
    float radius;
    std::vector<NodeVisibility> nodeVisibilities;
    std::mutex cellMutex;

    // Default constructor
    VisibilityCacheCell() : center(0.0f), radius(0.0f) {}

    // Delete copy (can't copy mutex)
    VisibilityCacheCell(const VisibilityCacheCell&) = delete;
    VisibilityCacheCell& operator=(const VisibilityCacheCell&) = delete;

    // Custom move constructor (move everything except mutex)
    VisibilityCacheCell(VisibilityCacheCell&& other) noexcept
        : center(other.center)
        , radius(other.radius)
        , nodeVisibilities(std::move(other.nodeVisibilities))
        , cellMutex()  // Construct new mutex (can't move)
    {}

    // Custom move assignment
    VisibilityCacheCell& operator=(VisibilityCacheCell&& other) noexcept {
        if (this != &other) {
            center = other.center;
            radius = other.radius;
            nodeVisibilities = std::move(other.nodeVisibilities);
            // cellMutex stays as-is (can't move)
        }
        return *this;
    }

    NodeVisibility* getNodeVisibility(int nodeIndex) {
        for (auto& nv : nodeVisibilities) {
            if (nv.nodeIndex == nodeIndex) {
                return &nv;
            }
        }
        return nullptr;
    }

    NodeVisibility* getOrCreateNodeVisibility(int nodeIndex) {
        std::lock_guard<std::mutex> lock(cellMutex);

        for (auto& nv : nodeVisibilities) {
            if (nv.nodeIndex == nodeIndex) {
                return &nv;
            }
        }

        nodeVisibilities.emplace_back();
        nodeVisibilities.back().nodeIndex = nodeIndex;
        return &nodeVisibilities.back();
    }
};


/// Spatial Visibility Cache
class SpatialVisibilityCache {
public:
    SpatialVisibilityCache(const glm::vec3& sceneMin, const glm::vec3& sceneMax,
                          int resolution = 16);

    // get the cell index for a point
    int getCellIndex(const glm::vec3& point) const;

    // get the cell for a point
    VisibilityCacheCell* getCell(const glm::vec3& point);

    // record a visibility sample
    void recordVisibility(const glm::vec3& shadingPoint, int nodeIndex, bool wasVisible);

    // query visibility probability
    [[nodiscard]] float queryVisibility(const glm::vec3& shadingPoint, int nodeIndex) const;

    // get statistics
    void getStatistics(int& totalCells, int& activeCells, int& totalSamples) const;

    // reset the cache (for adaptive refinement)
    void reset();

private:
    glm::vec3 m_sceneMin;
    glm::vec3 m_sceneMax;
    glm::vec3 m_cellSize{};
    glm::ivec3 m_resolution{};
    std::vector<VisibilityCacheCell> m_cells;

    [[nodiscard]] glm::ivec3 pointToCellCoord(const glm::vec3& point) const;
    [[nodiscard]] int coordToCellIndex(const glm::ivec3& coord) const;
    [[nodiscard]] bool isValidCoord(const glm::ivec3& coord) const;
};

/// MIS Weight Calculator
class MISWeightCalculator {
public:
    enum class Heuristic {
        Balance,    // w_i = pdf_i / sum(pdf_j)
        Power,      // w_i = pdf_i^2 / sum(pdf_j^2)
        Maximum     // w_i = 1 if i = argmax(pdf_j), else 0
    };

    static float calculateWeight(float thisPdf, float otherPdf,
                                Heuristic heuristic = Heuristic::Balance);

    static float balanceHeuristic(float thisPdf, float otherPdf);
    static float powerHeuristic(float thisPdf, float otherPdf, float beta = 2.0f);
};


/// Visibility-Aware Hierarchical Sampler with MIS
class VisibilityAwareHierarchicalSampler {
public:
    struct Sample {
        glm::vec3 position;
        glm::vec3 normal;
        glm::vec3 radiance;
        float pdf;
        float misWeight;
        bool isValid;
    };

    struct Configuration {
        bool enableVisibilityLearning = true;
        bool enableMIS = true;
        float visibilityWeight = 0.5f;      // 0 = pure flux, 1 = pure visibility
        int maxBVHDepth = 20;
        int minSamplesForLearning = 4;      // minimum samples before trusting visibility
        MISWeightCalculator::Heuristic misHeuristic = MISWeightCalculator::Heuristic::Balance;
    };

    struct X : Configuration {
    };

    VisibilityAwareHierarchicalSampler(MeshAreaLight* light,
                                       SpatialVisibilityCache* visCache,
                                       const Configuration& config = Configuration());

    // main sampling function
    [[nodiscard]] Sample sampleLight(const glm::vec3& shadingPoint,
                      const glm::vec3& normal,
                      float u1, float u2, float u3) const;

    // PDF evaluation
    [[nodiscard]] float evaluatePdf(const glm::vec3& shadingPoint,
                     const glm::vec3& lightPoint) const;

    // record visibility for learning
    void recordVisibilitySample(const glm::vec3& shadingPoint,
                               int nodeIndex,
                               bool wasVisible);

    // get configuration
    [[nodiscard]] const Configuration& getConfig() const { return m_config; }
    void setConfig(const Configuration& config) { m_config = config; }

private:
    MeshAreaLight* m_light;
    SpatialVisibilityCache* m_visCache;
    Configuration m_config;

    // hierarchical selection with visibility awareness
    int selectBVHNodeWithVisibility(int nodeIndex,
                                    const glm::vec3& shadingPoint,
                                    float u,
                                    float& outPdf) const;

    // calculate importance with visibility
    [[nodiscard]] float calculateNodeImportance(int nodeIndex,
                                 const glm::vec3& shadingPoint) const;

    // MIS combination
    [[nodiscard]] Sample combineWithMIS(const Sample& lightSample,
                         const glm::vec3& shadingPoint,
                         const glm::vec3& normal) const;
};


/// BSDF Sampler (for MIS)
class BSDFSampler {
public:
    struct BSDFSample {
        glm::vec3 direction;
        float pdf;
        glm::vec3 reflectance;
    };

    // sample diffuse BRDF (cosine-weighted hemisphere)
    static BSDFSample sampleDiffuse(const glm::vec3& normal, float u1, float u2);

    // evaluate diffuse BRDF PDF
    static float pdfDiffuse(const glm::vec3& normal, const glm::vec3& direction);

    // TODO: Add specular, glossy, etc.
};


/// Helper Functions
namespace VisibilityAwareHelpers {
    // convert between coordinate systems
    inline glm::vec3 uniformSampleHemisphere(float u1, float u2) {
        float z = u1;
        const float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
        const float phi = 2.0f * glm::pi<float>() * u2;
        return {r * std::cos(phi), r * std::sin(phi), z};
    }

    inline glm::vec3 cosineSampleHemisphere(float u1, float u2) {
        float z = std::sqrt(u1);
        const float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
        const float phi = 2.0f * glm::pi<float>() * u2;
        return {r * std::cos(phi), r * std::sin(phi), z};
    }

    inline glm::vec3 alignHemisphereWithNormal(const glm::vec3& sample,
                                               const glm::vec3& normal) {
        const glm::vec3 up = std::abs(normal.y) < 0.999f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
        const glm::vec3 tangent = glm::normalize(glm::cross(up, normal));
        const glm::vec3 bitangent = glm::cross(normal, tangent);
        return tangent * sample.x + bitangent * sample.y + normal * sample.z;
    }

    inline float powerHeuristic(float fPdf, float gPdf, float beta = 2.0f) {
        const float f = std::pow(fPdf, beta);
        const float g = std::pow(gPdf, beta);
        return f / (f + g);
    }
}


#endif //PROMETHEUS_VISIBILITY_AWARE_MIS_H