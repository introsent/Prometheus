//
// Created by ivans on 20/12/2025.
//

#include "area_light.h"

#include "triangle.h"
#include "render/scene_manager.h"

/// S1: Uniform Triangle Sampler implementation
UniformTriangleSampler::UniformTriangleSampler(const glm::vec3& v0,
                                               const glm::vec3& v1,
                                               const glm::vec3& v2,
                                               const glm::vec3& normal,
                                               float area,
                                               const glm::vec3& emission,
                                               float intensity)
    : m_v0(v0), m_v1(v1), m_v2(v2)
    , m_normal(normal)
    , m_area(area)
    , m_emission(emission)
    , m_intensity(intensity)
{}

AreaLightSample UniformTriangleSampler::sample(const glm::vec3& shadingPoint,
                                               float u1, float u2) const {
    // uniform sampling using SQUARE ROOT parameterization
    const float su0 = std::sqrt(u1);
    const float b0 = 1.0f - su0;
    const float b1 = u2 * su0;
    const float b2 = 1.0f - b0 - b1;

    // sample point on triangle using barycentric coordinates
    const glm::vec3 position = b0 * m_v0 + b1 * m_v1 + b2 * m_v2;

    // uniform PDF (1/Area)
    const float pdf = 1.0f / m_area;

    // radiance (emission * intensity)
    const glm::vec3 radiance = m_emission * m_intensity;

    return AreaLightSample{
        position,
        m_normal,
        pdf,
        radiance,
        m_area
    };
}

float UniformTriangleSampler::pdf(const glm::vec3& shadingPoint,
                                  const glm::vec3& lightPoint) const {
    // Uniform sampling: constant PDF = 1/Area
    return 1.0f / m_area;
}

float UniformTriangleSampler::getTotalFlux() const {
    // Power = Emission * Intensity * Area * π (for Lambertian)
    // For simplicity: Power = avg(emission) * intensity * area
    const float avgEmission = (m_emission.r + m_emission.g + m_emission.b) / 3.0f;
    return avgEmission * m_intensity * m_area;
}

/// Triangle Area light implementation
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
    Triangle* tri = m_scene->getTriangle(m_triangleIndex);
    if (!tri) {
        m_v0 = m_v1 = m_v2 = glm::vec3(0.0f);
        m_normal = glm::vec3(0.0f, 1.0f, 0.0f);
        m_area = 0.0f;
        return;
    }

    // get vertices from triangle
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

    // calculate normal
    const glm::vec3 edge1 = m_v1 - m_v0;
    const glm::vec3 edge2 = m_v2 - m_v0;
    m_normal = glm::normalize(glm::cross(edge1, edge2));

    // calculate area: 0.5 * |cross(edge1, edge2)|
    m_area = 0.5f * glm::length(glm::cross(edge1, edge2));
}

void TriangleAreaLight::initializeSampler() {
    switch (m_strategy) {
        case SamplingStrategy::Uniform:
            m_sampler = std::make_unique<UniformTriangleSampler>(
                m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
            break;

        case SamplingStrategy::AreaImportance:
            // TODO: Implement area importance sampler
            // for now, fall back to uniform
            m_sampler = std::make_unique<UniformTriangleSampler>(
                m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
            break;

        case SamplingStrategy::HierarchicalFlux:
            // TODO: Implement hierarchical flux sampler
            m_sampler = std::make_unique<UniformTriangleSampler>(
                m_v0, m_v1, m_v2, m_normal, m_area, m_emission, m_intensity);
            break;

        case SamplingStrategy::VisibilityAwareHierarchical:
            // TODO: Implement visibility aware sampler
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

glm::vec3 TriangleAreaLight::getNormal() const {
    return m_normal;
}

float TriangleAreaLight::getPower() const {
    return m_sampler->getTotalFlux();
}

void TriangleAreaLight::setSamplingStrategy(SamplingStrategy strategy) {
    if (m_strategy != strategy) {
        m_strategy = strategy;
        initializeSampler();
    }
}


/// Mesh Area light implementation
MeshAreaLight::MeshAreaLight(unsigned int meshIndex,
                            const glm::vec3& emission,
                            float intensity,
                            SceneManager* scene)
    : m_meshIndex(meshIndex)
    , m_scene(scene)
    , m_emission(emission)
    , m_intensity(intensity)
    , m_totalArea(0.0f)
{
    // get mesh and extract triangles
    Mesh* mesh = scene->getMesh(meshIndex);
    if (!mesh) return;

    const auto& vertices = mesh->getOriginalVertices();
    const auto& indices = mesh->getIndices();

    // create triangle lights for each face WITHOUT adding to scene
    for (size_t i = 0; i < indices.size(); i += 3) {
        const uint32_t idx0 = indices[i];
        const uint32_t idx1 = indices[i + 1];
        const uint32_t idx2 = indices[i + 2];

        // extract triangle vertices
        const glm::vec3 v0(vertices[idx0].position.x,
                          vertices[idx0].position.y,
                          vertices[idx0].position.z);
        const glm::vec3 v1(vertices[idx1].position.x,
                          vertices[idx1].position.y,
                          vertices[idx1].position.z);
        const glm::vec3 v2(vertices[idx2].position.x,
                          vertices[idx2].position.y,
                          vertices[idx2].position.z);

        // calculate normal and area directly
        const glm::vec3 edge1 = v1 - v0;
        const glm::vec3 edge2 = v2 - v0;
        const glm::vec3 crossProd = glm::cross(edge1, edge2);
        const float area = 0.5f * glm::length(crossProd);
        const glm::vec3 normal = (area > 0.0f) ? glm::normalize(crossProd) : glm::vec3(0, 1, 0);

        // create sampler directly (skip TriangleAreaLight wrapper for now)
        auto sampler = std::make_unique<UniformTriangleSampler>(
            v0, v1, v2, normal, area, emission, intensity
        );

        m_totalArea += area;
        m_triangleSamplers.push_back(std::move(sampler));
    }

    // build CDF for triangle selection (area-weighted)
    m_triangleCDF.resize(m_triangleSamplers.size());
    float cumulativeArea = 0.0f;
    for (size_t i = 0; i < m_triangleSamplers.size(); ++i) {
        // use the area from the sampler
        const float triArea = m_totalArea / m_triangleSamplers.size(); // approximate
        // better: store areas separately or get from sampler
        cumulativeArea += triArea;
        m_triangleCDF[i] = cumulativeArea;
    }

    // normalize CDF
    if (cumulativeArea > 0.0f) {
        for (float& val : m_triangleCDF) {
            val /= cumulativeArea;
        }
    }
}

AreaLightSample MeshAreaLight::sample(const glm::vec3& shadingPoint,
                                     float u1, float u2, float u3) const {
    if (m_triangleSamplers.empty()) {
        return AreaLightSample{
            glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f),
            0.0f, glm::vec3(0.0f), 0.0f
        };
    }

    // select triangle based on area-weighted CDF
    size_t triIndex = 0;
    for (size_t i = 0; i < m_triangleCDF.size(); ++i) {
        if (u1 <= m_triangleCDF[i]) {
            triIndex = i;
            break;
        }
    }

    // sample the selected triangle
    AreaLightSample sample = m_triangleSamplers[triIndex]->sample(shadingPoint, u2, u3);

    // adjust PDF for triangle selection probability
    const float triSelectionProb = 1.0f / m_triangleSamplers.size(); // uniform for now
    sample.pdf *= triSelectionProb;
    sample.area = m_totalArea;

    return sample;
}

float MeshAreaLight::getTotalPower() const {
    float power = 0.0f;
    for (const auto& triLight : m_triangleLights) {
        power += triLight->getPower();
    }
    return power;
}

void MeshAreaLight::setSamplingStrategy(SamplingStrategy strategy) const
{
    for (auto& triLight : m_triangleLights) {
        triLight->setSamplingStrategy(strategy);
    }
}
