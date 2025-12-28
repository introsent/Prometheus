//
// Created by ivans on 28/12/2025.
//

#include "triangle_area_light.h"

#include "sampler_base.h"
#include "render/scene_manager.h"
#include "samplers/area_importance_triangle_sampler.h"
#include "samplers/uniform_triangle_sampler.h"

TriangleAreaLight::TriangleAreaLight(unsigned int triangleIndex,
                                     const glm::vec3& emission,
                                     float intensity,
                                     SceneManager* scene)
    : m_triangleIndex(triangleIndex)
    , m_scene(scene)
    , m_emission(emission)
    , m_intensity(intensity)
    , m_strategy(SamplingStrategy::Uniform)
{
    cacheGeometry();
    initializeSampler();
}

void TriangleAreaLight::cacheGeometry() {
    const Triangle* tri = m_scene->getTriangle(m_triangleIndex);
    if (!tri) {
        m_v0 = m_v1 = m_v2 = glm::vec3(0.0f);
        m_normal = glm::vec3(0.0f, 1.0f, 0.0f);
        m_area = 0.0f;
        return;
    }

    const auto& vertices = tri->getOriginalVertices();
    if (vertices.size() < 3) {
        m_v0 = m_v1 = m_v2 = glm::vec3(0.0f);
        m_normal = glm::vec3(0.0f, 1.0f, 0.0f);
        m_area = 0.0f;
        return;
    }

    m_v0 = glm::vec3(vertices[0].position.x, vertices[0].position.y, vertices[0].position.z);
    m_v1 = glm::vec3(vertices[1].position.x, vertices[1].position.y, vertices[1].position.z);
    m_v2 = glm::vec3(vertices[2].position.x, vertices[2].position.y, vertices[2].position.z);

    const glm::vec3 edge1 = m_v1 - m_v0;
    const glm::vec3 edge2 = m_v2 - m_v0;
    m_normal = glm::normalize(glm::cross(edge1, edge2));
    m_area = 0.5f * glm::length(glm::cross(edge1, edge2));
}

void TriangleAreaLight::initializeSampler() {
    switch (m_strategy) {
        case SamplingStrategy::Uniform:
            m_sampler = std::make_unique<UniformTriangleSampler>(
                m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
            break;

        case SamplingStrategy::AreaImportance:
            m_sampler = std::make_unique<AreaImportanceTriangleSampler>(
                m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
            break;

        case SamplingStrategy::HierarchicalFlux:
        case SamplingStrategy::VisibilityAwareHierarchical:
            // TODO: Implement other strategies
            m_sampler = std::make_unique<UniformTriangleSampler>(
                m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
            break;
    }
}

AreaLightSample TriangleAreaLight::sample(const glm::vec3& shadingPoint,
                                         float u1, float u2) const {
    return m_sampler->sample(shadingPoint, u1, u2);
}

float TriangleAreaLight::pdf(const glm::vec3& shadingPoint,
                            const glm::vec3& lightPoint) const {
    return m_sampler->pdf(shadingPoint, lightPoint);
}

float TriangleAreaLight::getPower() const {
    return m_sampler->getTotalFlux();
}

float TriangleAreaLight::calculateSolidAngle(const glm::vec3& p) const {
    if (m_strategy == SamplingStrategy::AreaImportance) {
        if (const auto* importanceSampler = dynamic_cast<AreaImportanceTriangleSampler*>(m_sampler.get())) {
            return importanceSampler->calculateSolidAngle(p);
        }
    }
    // For uniform sampling, we still need solid angle for mesh area importance
    auto tempSampler = AreaImportanceTriangleSampler(
        m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
    return tempSampler.calculateSolidAngle(p);
}

void TriangleAreaLight::setSamplingStrategy(SamplingStrategy strategy) {
    if (m_strategy != strategy) {
        m_strategy = strategy;
        initializeSampler();
    }
}