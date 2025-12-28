//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_TRIANGLE_AREA_LIGHT_H
#define PROMETHEUS_TRIANGLE_AREA_LIGHT_H
#include <glm/vec3.hpp>
#include "sampler_base.h"
#include <memory>

class SceneManager;
class Triangle;
class Mesh;

/// Triangle area light
// single triangle emitting light uniformly or with importance sampling
class TriangleAreaLight {
public:
    TriangleAreaLight(
        unsigned int triangleIndex,
        const glm::vec3& emission,
        float intensity,
        SceneManager* scene);

    // sample a point on the light surface
    [[nodiscard]] AreaLightSample sample(
        const glm::vec3& shadingPoint,
        float u1, float u2) const;

    // evaluate pdf for a given light point
    [[nodiscard]] float pdf(
        const glm::vec3& shadingPoint,
        const glm::vec3& lightPoint) const;

    // accessors
    [[nodiscard]] glm::vec3 getNormal() const { return m_normal; }
    [[nodiscard]] float getArea() const { return m_area; }
    [[nodiscard]] glm::vec3 getEmission() const { return m_emission; }
    [[nodiscard]] float getIntensity() const { return m_intensity; }
    [[nodiscard]] float getPower() const;
    [[nodiscard]] unsigned int getTriangleIndex() const { return m_triangleIndex; }

    // strategy management
    void setSamplingStrategy(SamplingStrategy strategy);
    [[nodiscard]] SamplingStrategy getSamplingStrategy() const { return m_strategy; }

    // for hierarchical sampling in mesh lights
    [[nodiscard]] float calculateSolidAngle(const glm::vec3& p) const;

private:
    void initializeSampler();
    void cacheGeometry();

    unsigned int m_triangleIndex;
    SceneManager* m_scene;

    // light properties
    glm::vec3 m_emission;
    float m_intensity;

    // cached geometry
    glm::vec3 m_v0{}, m_v1{}, m_v2{};
    glm::vec3 m_normal{};
    float m_area{};

    // sampling strategy
    SamplingStrategy m_strategy{};
    std::unique_ptr<IAreaLightSampler> m_sampler{};
};

#endif //PROMETHEUS_TRIANGLE_AREA_LIGHT_H