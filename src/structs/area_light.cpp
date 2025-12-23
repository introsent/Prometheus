//
// Created by ivans on 20/12/2025.
//

#include "area_light.h"

#include "hit_result.h"
#include "ray.h"
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
    // uniform sampling: constant PDF = 1/Area
    return 1.0f / m_area;
}

float UniformTriangleSampler::getTotalFlux() const {
    // power = emission * intensity * area * pi (for Lambertian)
    // for simplicity: power = avg(emission) * intensity * area
    const float avgEmission = (m_emission.r + m_emission.g + m_emission.b) / 3.0f;
    return avgEmission * m_intensity * m_area;
}

/// S2: Area Importance Triangle Sampler implementation
AreaImportanceTriangleSampler::AreaImportanceTriangleSampler(const glm::vec3& v0,
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

AreaLightSample AreaImportanceTriangleSampler::sample(const glm::vec3& shadingPoint,
                                               float u1, float u2) const {
    // sample direction using solid angle (spherical triangle projection)
    float pdfSolidAngle = 0.0f;
    glm::vec3 direction = sampleSphericalTriangle(shadingPoint, u1, u2, &pdfSolidAngle);

    if (pdfSolidAngle <= 0.0f) {
        // degenerate: fallback to centroid

        // compute centroid as average of vertices for fallback position
        glm::vec3 centroid = (m_v0 + m_v1 + m_v2) / 3.0f;
        glm::vec3 toCentroid = centroid - shadingPoint;
        float r = glm::length(toCentroid);
        if (r <= 0.0f) return {shadingPoint, m_normal, 0.0f, glm::vec3(0), m_area};

        direction = toCentroid / r;  // normalize direction to unit vector
        pdfSolidAngle = 1.0f;  // arbitrary for fallback
    }

    // compute barycentric coordinates from direction (stable projection back to plane)
    glm::vec3 bary = computeBarycentricsFromDirection(shadingPoint, direction);
    float sum = bary.x + bary.y + bary.z;  // sum of barycentric coordinates (should ideally be 1)
    if (sum > 0.0f) {
        bary /= sum;  // normalize by dividing each coordinate by the sum to ensure they add to 1
    } else {
        // rare degenerate case: fallback to centroid
        bary = glm::vec3(1.0f / 3.0f);  // uniform barycentric coordinates for centroid
    }

    // compute position on triangle
    glm::vec3 position = bary.x * m_v0 + bary.y * m_v1 + bary.z * m_v2;  // weighted sum of vertices using barycentric coordinates

    glm::vec3 fromLightToShading = shadingPoint - position;
    float r = glm::length(fromLightToShading);
    if (r <= 0.0f) return {position, m_normal, 0.0f, glm::vec3(0), m_area};  // degenerate

    glm::vec3 dirToShading = fromLightToShading / r;
    // cosine of angle between light normal and direction to shading point
    float cos_theta_light = glm::dot(m_normal, dirToShading);

    if (cos_theta_light <= 0.0f) return {position, m_normal, 0.0f, glm::vec3(0), m_area};

    // convert solid-angle pdf to area pdf using jacobian: dA = dw * (r^2 / cos_theta_light)
    float pdfArea = pdfSolidAngle * (cos_theta_light / (r * r));

    return {
        position,
        m_normal,
        pdfArea,
        m_emission * m_intensity,
        m_area
    };
}

float AreaImportanceTriangleSampler::pdf(const glm::vec3& shadingPoint,
                                  const glm::vec3& lightPoint) const {
    // compute solid angle subtended by triangle from shadingPoint
    const float solidAngle = calculateSolidAngle(shadingPoint);
    if (solidAngle <= 0.0f) return 0.0f;

    // vector from light point to shading point
    glm::vec3 fromLightToShading = shadingPoint - lightPoint;
    float r = glm::length(fromLightToShading);
    if (r <= 0.0f) return 0.0f;

    // normalize to unit vector
    glm::vec3 dirToShading = fromLightToShading / r;

    // cosine of angle at light point
    float cos_theta_light = glm::dot(m_normal, dirToShading);
    if (cos_theta_light <= 0.0f) return 0.0f;

    // uniform pdf over solid angle: 1 / omega
    float pdfSolidAngle = 1.0f / solidAngle;
    // transform to area pdf: pdf_area = pdf_solid * (cos_theta_light / r^2)
    return pdfSolidAngle * (cos_theta_light / (r * r));
}

float AreaImportanceTriangleSampler::getTotalFlux() const {
    // power = emission * intensity * area * pi (for Lambertian)
    // for simplicity: power = avg(emission) * intensity * area
    const float avgEmission = (m_emission.r + m_emission.g + m_emission.b) / 3.0f;  // average rgb components of emission for scalar approximation
    return avgEmission * m_intensity * m_area;  // total flux as product of average emission, intensity, and area (simplified, ignoring pi for lambertian)
}

float AreaImportanceTriangleSampler::calculateSolidAngle(const glm::vec3 &p) const {
    // van oosterom & strackee formula for solid angle of a triangle
    // ref: "The Solid Angle of a Plane Triangle" (1983)

    // transform triangle vertices to vectors from point p
    const glm::vec3 a = m_v0 - p;
    const glm::vec3 b = m_v1 - p;
    const glm::vec3 c = m_v2 - p;

    // compute lengths (distances from p to each vertex)
    const float la = glm::length(a);
    const float lb = glm::length(b);
    const float lc = glm::length(c);

    // handle degenerate cases: if point is on or very close to triangle
    if (la == 0.0f || lb == 0.0f || lc == 0.0f) {
        return 0.0f;  // point is on a vertex - triangle subtends no solid angle
    }

    // normalize vectors to get unit vectors
    const glm::vec3 a_unit = a / la;
    const glm::vec3 b_unit = b / lb;
    const glm::vec3 c_unit = c / lc;

    // compute scalar triple product: a dot (b cross c)
    // in result, the signed volume of the parallelepiped spanned by a, b, c
    const float tripleProduct = glm::dot(a_unit, glm::cross(b_unit, c_unit));

    // compute dot products between unit vectors
    const float dot_ab = glm::dot(a_unit, b_unit);
    const float dot_bc = glm::dot(b_unit, c_unit);
    const float dot_ca = glm::dot(c_unit, a_unit);

    // denominator for the formula
    const float denominator = 1.0f + dot_ab + dot_bc + dot_ca;

    // if denominator is zero or negative, triangle is back-facing or degenerate
    if (denominator <= 0.0f || std::abs(tripleProduct) < 1e-10f) {
        // check triangle orientation using signed triple product
        if (const float signedVolume = glm::dot(a, glm::cross(b, c)); signedVolume <= 0.0f) {
            // triangle is back-facing from point p
            return 0.0f;
        }
        // for degenerate cases with positive volume but zero triple_product,
        // returned a small epsilon to avoid division by zero
        return 1e-10f;
    }
    // van oosterom & strackee formula:
    // solidAngle =  2 * arctan( |a·(b×c)| / (1 + a·b + b·c + c·a) )
    // where a, b, c are unit vectors from p to triangle vertices

    // use atan2(y, x) for numerical stability (handles cases where x limit is 0 ( x-> 0) )
    const float solidAngle = 2.0f * std::atan2(std::abs(tripleProduct), denominator);

    // clamp to valid range [0, 2 pi]
    return std::max(0.0f, std::min(solidAngle, 2.0f * glm::pi<float>()));
}

glm::vec3 AreaImportanceTriangleSampler::sampleSphericalTriangle(const glm::vec3 &p, float u1, float u2,
    float *outPdf) const {
        glm::vec3 a = m_v0 - p;
        glm::vec3 b = m_v1 - p;
        glm::vec3 c = m_v2 - p;

        float la = glm::length(a), lb = glm::length(b), lc = glm::length(c);  // lengths of vectors to vertices
        if (la <= 0.0f || lb <= 0.0f || lc <= 0.0f) {
            if (outPdf) *outPdf = 0.0f;
            return glm::normalize(m_v0 - p);  // fallback to normalized direction to vertex 0
        }
        a /= la; b /= lb; c /= lc;  // normalize vectors to unit sphere

        // normal to great circle arcs
        glm::vec3 n_ab = glm::cross(a, b);
        glm::vec3 n_bc = glm::cross(b, c);
        glm::vec3 n_ca = glm::cross(c, a);

        constexpr float eps = 1e-10f;
        if (glm::dot(n_ab, n_ab) <= eps || glm::dot(n_bc, n_bc) <= eps || glm::dot(n_ca, n_ca) <= eps) {  // check for degenerate arcs (near-zero length)
            if (outPdf) *outPdf = 0.0f;
            return glm::normalize(a + b + c);  // fallback to normalized average direction
        }

        // unit normal for arc
        n_ab = glm::normalize(n_ab);
        n_bc = glm::normalize(n_bc);
        n_ca = glm::normalize(n_ca);

        // spherical angle at vertices (arvo's method: angle between planes)
        float alpha = safeAcos(glm::dot(n_ab, -n_ca));
        float beta = safeAcos(glm::dot(n_bc, -n_ab));
        float gamma = safeAcos(glm::dot(n_ca, -n_bc));

        // sum of spherical angles
        float A_pi = alpha + beta + gamma;
        // spherical excess: solid angle omega = alpha + beta + gamma - pi
        float Omega = A_pi - glm::pi<float>();
        if (Omega <= 0.0f) {
            if (outPdf) *outPdf = 0.0f;
            return glm::normalize(a + b + c);  // degenerate solid angle
        }
        // uniform pdf over solid angle
        if (outPdf) *outPdf = 1.0f / Omega;

        // map u1 to angle in [pi, A_pi] for marginal sampling
        float Ap_pi = glm::pi<float>() + u1 * (A_pi - glm::pi<float>());

        // precompute trig for alpha
        float cosAlpha = std::cos(alpha);
        float sinAlpha = std::sin(alpha);

        // trig for adjusted angle Ap_pi
        float sinAp = std::sin(Ap_pi);
        float cosAp = std::cos(Ap_pi);

        // sin(phi) using angle subtraction formula
        float sinPhi = sinAp * cosAlpha - cosAp * sinAlpha;
        // cos(phi) using angle addition
        float cosPhi = cosAp * cosAlpha + sinAp * sinAlpha;

        // cosine of spherical side c (angle between a and b)
        float cosc = glm::dot(a, b);

        float k1 = cosPhi + cosAlpha;
        float k2 = sinPhi - sinAlpha * cosc;
        // denominator for cos(Bp) in spherical law of cosines
        float denominator = (k2 * sinPhi + k1 * cosPhi) * sinAlpha;
        // numerator for cos(Bp)
        float numerator = k2 + (k2 * cosPhi - k1 * sinPhi) * cosAlpha;

        float cosBp = (std::abs(denominator) < 1e-10f) ? clampf(cosc * (1.0f - u1) + glm::dot(c, b) * u1, -1.0f, 1.0f)
                                                 : clampf(numerator / denominator, -1.0f, 1.0f);
        // cos of angle Bp, with fallback linear interp

        float sinBp = safeSqrt(1.0f - cosBp * cosBp);  // sin(Bp) from trig identity

        // orthogonalize to get basis vector perpendicular to a in plane of a and c
        glm::vec3 axis_ac = gramSchmidtNormalize(c, a);
        // point cp on arc from a to projection
        glm::vec3 cp = cosBp * a + sinBp * axis_ac;
        cp = glm::normalize(cp);

        // cosine of angle between cp and b
        float dot_cp_b = glm::dot(cp, b);
        // map u2 to cos(theta) along arc
        float cosTheta = clampf(1.0f - u2 * (1.0f - dot_cp_b), -1.0f, 1.0f);
        float sinTheta = safeSqrt(1.0f - cosTheta * cosTheta);

        // orthogonal basis for rotation around b
        glm::vec3 axis_bp = gramSchmidtNormalize(cp, b);
        // rotate around b by theta
        glm::vec3 w = cosTheta * b + sinTheta * axis_bp;
        return glm::normalize(w);
}

glm::vec3 AreaImportanceTriangleSampler::computeBarycentricsFromDirection(const glm::vec3 &origin,
const glm::vec3 &dir) const {
    // moller-trumbore barycentrics (without t, since intersection is guaranteed by sampling)
    const glm::vec3 e1 = m_v1 - m_v0;  // edge vector from v0 to v1
    const glm::vec3 e2 = m_v2 - m_v0;  // edge vector from v0 to v2

    const glm::vec3 s1 = glm::cross(dir, e2);  // cross product for perpendicular vector (part of determinant computation)
    const float divisor = glm::dot(s1, e1);  // determinant of matrix [dir, e2, e1] for barycentric solve
    if (std::abs(divisor) < 1e-10f) {
        // degenerate: fallback to uniform bary (centroid)
        return glm::vec3(1.0f / 3.0f);
    }
    const float invDivisor = 1.0f / divisor;  // inverse determinant for scaling

    const glm::vec3 s = origin - m_v0;  // vector from v0 to ray origin
    float b1 = glm::dot(s, s1) * invDivisor;  // barycentric u = det([s, dir, e2]) / det
    float b2 = glm::dot(dir, glm::cross(s, e1)) * invDivisor;  // barycentric v = det([s, e1, dir]) / det (note cross order for sign)

    b1 = clampf(b1, 0.0f, 1.0f);  // clamp to triangle bounds
    b2 = clampf(b2, 0.0f, 1.0f);  // clamp to triangle bounds
    if (b1 + b2 > 1.0f) {
        const float sum = b1 + b2;
        b1 /= sum;  // renormalize if outside (project to edge)
        b2 /= sum;
    }

    return {1.0f - b1 - b2, b1, b2};  // b0 = 1 - b1 - b2
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
    , m_strategy(SamplingStrategy::AreaImportance)
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
            m_sampler = std::make_unique<AreaImportanceTriangleSampler>(
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
