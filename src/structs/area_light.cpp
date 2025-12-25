//
// Created by ivans on 20/12/2025.
//

#include "area_light.h"

#include <iostream>

#include "hit_result.h"
#include "ray.h"
#include "triangle.h"
#include "render/scene_manager.h"
#include "math_helpers.h"

// static helpers
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
    const glm::vec3 direction = sampleSphericalTriangle(shadingPoint, u1, u2, &pdfSolidAngle);

    if (pdfSolidAngle <= 0.0f) {
        return {shadingPoint, m_normal, 0.0f, glm::vec3(0), m_area};
    }

    // intersect ray with triangle plane to get exact position
    auto [t, position] = getRayPlaneIntersection(shadingPoint, direction);
    if (position == glm::vec3(0.f)) {
        return {shadingPoint, m_normal, 0.0f, glm::vec3(0), m_area};
    }
    //recomputePositionUsingBarycentricCoordinates(t, position);

    // use the ray parameter t directly for distance (more numerically stable)
    float r = t;
    if (r <= 1e-6f) return {position, m_normal, 0.0f, glm::vec3(0), m_area};

    // cosine at light (direction is already normalized)
    const float cos_theta_light = glm::dot(m_normal, -direction);

    if (cos_theta_light <= 0.0f) return {position, m_normal, 0.0f, glm::vec3(0), m_area};

    // convert solid-angle pdf to area pdf using jacobian
    float pdfArea = pdfSolidAngle * (cos_theta_light / (r * r));
    if (!std::isfinite(pdfArea) || pdfArea <= 0.0f) {
        return {position, m_normal, 0.0f, glm::vec3(0), m_area};
    }

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
        float alpha = MathHelpers::safeAcos(glm::dot(n_ab, -n_ca));
        float beta = MathHelpers::safeAcos(glm::dot(n_bc, -n_ab));
        float gamma = MathHelpers::safeAcos(glm::dot(n_ca, -n_bc));

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

        float cosBp = (std::abs(denominator) < 1e-10f) ? std::clamp(cosc * (1.0f - u1) + glm::dot(c, b) * u1, -1.0f, 1.0f)
                                                 : std::clamp(numerator / denominator, -1.0f, 1.0f);
        // cos of angle Bp, with fallback linear interp

        float sinBp = MathHelpers::safeSqrt(1.0f - cosBp * cosBp);  // sin(Bp) from trig identity

        // orthogonalize to get basis vector perpendicular to a in plane of a and c
        glm::vec3 axis_ac = MathHelpers::gramSchmidtNormalize(c, a);
        // point cp on arc from a to projection
        glm::vec3 cp = cosBp * a + sinBp * axis_ac;
        cp = glm::normalize(cp);

        // cosine of angle between cp and b
        float dot_cp_b = glm::dot(cp, b);
        // map u2 to cos(theta) along arc
        float cosTheta = std::clamp(1.0f - u2 * (1.0f - dot_cp_b), -1.0f, 1.0f);
        float sinTheta = MathHelpers::safeSqrt(1.0f - cosTheta * cosTheta);

        // orthogonal basis for rotation around b
        glm::vec3 axis_bp = MathHelpers::gramSchmidtNormalize(cp, b);
        // rotate around b by theta
        glm::vec3 w = cosTheta * b + sinTheta * axis_bp;
        return glm::normalize(w);
}

std::pair<float, glm::vec3> AreaImportanceTriangleSampler::getRayPlaneIntersection(const glm::vec3& shadingPoint, const glm::vec3& direction) const
{
    float denominator = glm::dot(m_normal, direction);
    if (std::abs(denominator) < 1e-6f) {
        // ray parallel to plane - degenerate
        return {0.f, {0.f, 0.f, 0.f}};
    }

    float t = glm::dot(m_v0 - shadingPoint, m_normal) / denominator;
    if (t <= 0.0f) {
        // intersection behind origin - should not happen but safety check
        return {0.f, {0.f, 0.f, 0.f}};
    }

    // compute exact intersection point
    return {t, shadingPoint + t * direction};
}

void AreaImportanceTriangleSampler::recomputePositionUsingBarycentricCoordinates(float t, glm::vec3 &position) const {
    // compute barycentric coordinates of position
    glm::vec3 v0v1 = m_v1 - m_v0;
    glm::vec3 v0v2 = m_v2 - m_v0;
    glm::vec3 v0p = position - m_v0;

    float d00 = glm::dot(v0v1, v0v1);
    float d01 = glm::dot(v0v1, v0v2);
    float d11 = glm::dot(v0v2, v0v2);
    float d20 = glm::dot(v0p, v0v1);
    float d21 = glm::dot(v0p, v0v2);

    float invDenom = 1.0f / (d00 * d11 - d01 * d01);
    float b1 = (d11 * d20 - d01 * d21) * invDenom;
    float b2 = (d00 * d21 - d01 * d20) * invDenom;
    float b0 = 1.0f - b1 - b2;

    // clamp to triangle (for numerical safety)
    b0 = std::max(0.0f, std::min(1.0f, b0));
    b1 = std::max(0.0f, std::min(1.0f, b1));
    b2 = std::max(0.0f, std::min(1.0f, b2));

    float sum = b0 + b1 + b2;
    if (sum > 0.0f) {
        b0 /= sum;
        b1 /= sum;
        b2 /= sum;
    }

    // recompute position with clamped barycentrics for consistency
    position = b0 * m_v0 + b1 * m_v1 + b2 * m_v2;
}

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

glm::vec3 TriangleAreaLight::getNormal() const {
    return m_normal;
}

float TriangleAreaLight::getPower() const {
    return m_sampler->getTotalFlux();
}

float TriangleAreaLight::calculateSolidAngle(const glm::vec3& p) const {
    if (m_strategy == SamplingStrategy::AreaImportance) {
        if (auto* importanceSampler = dynamic_cast<AreaImportanceTriangleSampler*>(m_sampler.get())) {
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

/// Mesh Area light implementation
MeshAreaLight::MeshAreaLight(unsigned int meshIndex,
                            const glm::vec3& emission,
                            float intensity,
                            SceneManager* scene)
    : m_meshIndex(meshIndex)
    , m_scene(scene)
    , m_emission(emission)
    , m_intensity(intensity)
    , m_strategy(SamplingStrategy::HierarchicalFlux)
    , m_totalArea(0.0f),
    m_rootNodeIndex(-1)
{
    Mesh* mesh = scene->getMesh(meshIndex);
    if (!mesh) return;

    const auto& vertices = mesh->getOriginalVertices();
    const auto& indices = mesh->getIndices();

    const size_t numTriangles = indices.size() / 3;
    m_triangles.reserve(numTriangles);
    m_areaCDF.reserve(numTriangles);

    for (size_t i = 0; i < indices.size(); i += 3) {
        const uint32_t idx0 = indices[i];
        const uint32_t idx1 = indices[i + 1];
        const uint32_t idx2 = indices[i + 2];

        const glm::vec3 v0(vertices[idx0].position.x, vertices[idx0].position.y, vertices[idx0].position.z);
        const glm::vec3 v1(vertices[idx1].position.x, vertices[idx1].position.y, vertices[idx1].position.z);
        const glm::vec3 v2(vertices[idx2].position.x, vertices[idx2].position.y, vertices[idx2].position.z);

        const glm::vec3 edge1 = v1 - v0;
        const glm::vec3 edge2 = v2 - v0;
        const glm::vec3 crossProd = glm::cross(edge1, edge2);
        const float area = 0.5f * glm::length(crossProd);

        if (area <= 1e-8f) continue;

        const glm::vec3 normal = glm::normalize(crossProd);

        m_triangles.emplace_back(v0, v1, v2, normal, area, emission, intensity);
        m_totalArea += area;
    }

    // build area CDF for uniform sampling
    if (!m_triangles.empty() && m_totalArea > 0.0f) {
        float cumulativeArea = 0.0f;
        for (const auto& tri : m_triangles) {
            cumulativeArea += tri.area;
            m_areaCDF.push_back(cumulativeArea / m_totalArea);
        }
        m_areaCDF.back() = 1.0f; // Ensure exact 1.0
    }

    // build BVH for hierarchical flux sampling
    buildBVH();

    // optional: Print BVH statistics
    printBVHStats();
}

void MeshAreaLight::buildBVH() {
    if (m_triangles.empty()) return;

    // prepare triangles info for BVH construction
    std::vector<TriangleInfo> triangleInfos;
    triangleInfos.reserve(m_triangles.size() * 2);

    for (size_t i = 0; i < m_triangles.size(); ++i) {
        const auto& triangle = m_triangles[i];

        // compute bounding box
        glm::vec3 bboxMin {};
        glm::vec3 bboxMax {};
        MathHelpers::computeTriangleBoundingBox(triangle.v0, triangle.v1, triangle.v2, bboxMin, bboxMax);

        // compute centroid
        glm::vec3 centroid = MathHelpers::computeTriangleCentroid(triangle.v0, triangle.v1, triangle.v2);

        // compute flux
        float flux = MathHelpers::computeTriangleFlux(m_emission, m_intensity, triangle.area);

        triangleInfos.emplace_back(centroid, bboxMin, bboxMax, flux, triangle.area, i);
    }

    m_rootNodeIndex = buildBVHNode(triangleInfos, 0, static_cast<int>(triangleInfos.size()));

    calculateNodeFlux(m_rootNodeIndex);
}

// Print BVH values to check tree correctness
void MeshAreaLight::printBVHStats() const {
    if (m_bvhNodes.empty()) return;

    int leafCount = 0;
    int maxDepth = 0;
    const float totalFlux = m_bvhNodes[m_rootNodeIndex].totalFlux;

    // simple traversal to collect stats
    std::function<void(int, int)> traverse = [&](int nodeIndex, int depth) {
        maxDepth = std::max(maxDepth, depth);

        if (const BVHNode& node = m_bvhNodes[nodeIndex]; node.isLeaf) {
            leafCount++;
        } else {
            traverse(node.leftChild, depth + 1);
            traverse(node.rightChild, depth + 1);
        }
    };

    traverse(m_rootNodeIndex, 0);

    std::cout << "BVH Statistics:\n";
    std::cout << "  Total nodes: " << m_bvhNodes.size() << "\n";
    std::cout << "  Leaf nodes: " << leafCount << "\n";
    std::cout << "  Max depth: " << maxDepth << "\n";
    std::cout << "  Total flux: " << totalFlux << "\n";
}

int MeshAreaLight::buildBVHNode(std::vector<TriangleInfo> &triInfos, int start, int end, int depth) {
    constexpr int MAX_TRIANGLES_PER_LEAF = 4;
    constexpr int MAX_DEPTH = 20;

    // create new node
    int nodeIndex = static_cast<int>(m_bvhNodes.size());
    m_bvhNodes.emplace_back();

    // compute bounding box for all triangles in this node
    auto bboxMin = glm::vec3(FLT_MAX);
    auto bboxMax = glm::vec3(-FLT_MAX);

    for (int i = start; i < end; ++i) {
        MathHelpers::expandBoundingBox(bboxMin, bboxMax, triInfos[i].bboxMin);
        MathHelpers::expandBoundingBox(bboxMin, bboxMax, triInfos[i].bboxMax);
    }

    // Update node using index (safe from reallocation)
    m_bvhNodes[nodeIndex].bboxMin = bboxMin;
    m_bvhNodes[nodeIndex].bboxMax = bboxMax;

    // check termination criteria
    int numTriangles = end - start;
    if (numTriangles <= MAX_TRIANGLES_PER_LEAF || depth >= MAX_DEPTH) {
        // create leaf node
        m_bvhNodes[nodeIndex].startTri = start;
        m_bvhNodes[nodeIndex].endTri = end;
        m_bvhNodes[nodeIndex].isLeaf = true;
        m_bvhNodes[nodeIndex].totalArea = 0.0f;

        // store the original triangle indices
        m_bvhNodes[nodeIndex].triangleIndices.reserve(numTriangles);
        for (int i = start; i < end; ++i) {
            m_bvhNodes[nodeIndex].triangleIndices.push_back(triInfos[i].originalIndex);
            m_bvhNodes[nodeIndex].totalArea += triInfos[i].area;
        }
        return nodeIndex;
    }

    // choose split axis (longest axis)
    glm::vec3 diagonal = bboxMax - bboxMin;
    int splitAxis = 0;
    if (diagonal.y > diagonal.x) splitAxis = 1;
    if (diagonal.z > diagonal[splitAxis]) splitAxis = 2;

    // sort triangles by centroid along split axis
    std::sort(triInfos.begin() + start, triInfos.begin() + end,
              [splitAxis](const TriangleInfo& a, const TriangleInfo& b) {
                  return a.centroid[splitAxis] < b.centroid[splitAxis];
              });

    // split at median
    int mid = start + numTriangles / 2;

    // recursively build children
    int leftChild = buildBVHNode(triInfos, start, mid, depth + 1);
    int rightChild = buildBVHNode(triInfos, mid, end, depth + 1);

    m_bvhNodes[nodeIndex].leftChild = leftChild;
    m_bvhNodes[nodeIndex].rightChild = rightChild;
    m_bvhNodes[nodeIndex].isLeaf = false;

    return nodeIndex;
}

void MeshAreaLight::calculateNodeFlux(int nodeIndex) {
    BVHNode& node = m_bvhNodes[nodeIndex];

    if (node.isLeaf) {
        node.totalFlux = 0.f;
        for (int originalIdx : node.triangleIndices) {
            const auto& triangle = m_triangles[originalIdx];
            node.totalFlux += MathHelpers::computeTriangleFlux(m_emission, m_intensity, triangle.area);
        }
    } else {
        // internal node: recursively calculate children, then sum
        calculateNodeFlux(node.leftChild);
        calculateNodeFlux(node.rightChild);

        // only sum if both children are valid
        if (node.leftChild >= 0 && node.rightChild >= 0 &&
            node.leftChild < static_cast<int>(m_bvhNodes.size()) &&
            node.rightChild < static_cast<int>(m_bvhNodes.size())) {

            const BVHNode& left = m_bvhNodes[node.leftChild];
            const BVHNode& right = m_bvhNodes[node.rightChild];

            node.totalFlux = left.totalFlux + right.totalFlux;
            node.totalArea = left.totalArea + right.totalArea;
            } else {
                // Fallback: initialize to 0
                node.totalFlux = 0.0f;
                node.totalArea = 0.0f;
                std::cerr << "ERROR: Cannot sum invalid children for node " << nodeIndex << std::endl;
            }
    }
}


int MeshAreaLight::selectBVHNode(int nodeIndex, float u, const glm::vec3& shadingPoint, float& outPdf) const {
    const BVHNode& node = m_bvhNodes[nodeIndex];

    // Base case: Leaf node
    if (node.isLeaf) {
        return nodeIndex;
    }

    const BVHNode& left = m_bvhNodes[node.leftChild];
    const BVHNode& right = m_bvhNodes[node.rightChild];

    // 1. Calculate Distance Squared to centroids
    // Using a small epsilon (1e-4f) to prevent division by zero if shadingPoint is inside the node
    float distSqLeft = glm::distance(shadingPoint, left.getCentroid());
    distSqLeft = std::max(distSqLeft * distSqLeft, 1e-4f);

    float distSqRight = glm::distance(shadingPoint, right.getCentroid());
    distSqRight = std::max(distSqRight * distSqRight, 1e-4f);

    // 2. Calculate Importance Weights (Flux / Distance^2)
    // This estimates how much light reaches the shading point from this node
    float weightLeft = left.totalFlux / distSqLeft;
    float weightRight = right.totalFlux / distSqRight;
    float totalWeight = weightLeft + weightRight;

    if (totalWeight <= 0.0f) {
        return -1; // Should not happen if flux > 0
    }

    // 3. Calculate Probability of going Left
    float probLeft = weightLeft / totalWeight;

    // 4. Traverse and Update PDF
    if (u < probLeft) {
        // Update the accumulated PDF for the path
        outPdf *= probLeft;

        // Rescale u for the next level
        float newU = (probLeft > 0.0f) ? u / probLeft : 0.0f;
        return selectBVHNode(node.leftChild, newU, shadingPoint, outPdf);
    } else {
        // Probability of going Right is (1.0 - probLeft)
        float probRight = 1.0f - probLeft;
        outPdf *= probRight;

        // Rescale u
        float newU = (probRight > 0.0f) ? (u - probLeft) / probRight : 0.0f;
        return selectBVHNode(node.rightChild, newU, shadingPoint, outPdf);
    }
}

bool MeshAreaLight::isPointInTriangle(const glm::vec3& p,
                                      const glm::vec3& v0,
                                      const glm::vec3& v1,
                                      const glm::vec3& v2,
                                      float& u, float& v, float& w) const {
    // Compute vectors
    glm::vec3 v0v1 = v1 - v0;
    glm::vec3 v0v2 = v2 - v0;
    glm::vec3 v0p = p - v0;

    // Compute dot products
    float d00 = glm::dot(v0v1, v0v1);
    float d01 = glm::dot(v0v1, v0v2);
    float d11 = glm::dot(v0v2, v0v2);
    float d20 = glm::dot(v0p, v0v1);
    float d21 = glm::dot(v0p, v0v2);

    float denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < 1e-8f) return false;

    // Compute barycentric coordinates
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;

    // Check if point is in triangle
    constexpr float eps = 1e-4f;
    return (u >= -eps) && (v >= -eps) && (w >= -eps);
}

AreaLightSample MeshAreaLight::sampleUniform(const glm::vec3& shadingPoint,
                                             float u1, float u2, float u3) const {
    if (m_triangles.empty()) {
        return AreaLightSample{};
    }

    // 1. select triangle by area
    const auto it = std::lower_bound(m_areaCDF.begin(), m_areaCDF.end(), u1);
    size_t triIndex = std::distance(m_areaCDF.begin(), it);
    triIndex = std::min(triIndex, m_triangles.size() - 1);

    // 2. sample uniformly on selected triangle
    const TriangleData& tri = m_triangles[triIndex];
    AreaLightSample sample = tri.uniformSampler->sample(shadingPoint, u2, u3);

    // 3. adjust PDF for triangle selection
    float selectionProb = tri.area / m_totalArea;
    sample.pdf *= selectionProb;
    sample.area = m_totalArea;

    return sample;
}

AreaLightSample MeshAreaLight::sampleAreaImportance(const glm::vec3& shadingPoint,
                                                   float u1, float u2, float u3) const {
    if (m_triangles.empty()) {
        return AreaLightSample{};
    }

    // 1. compute solid angles for all triangles
    std::vector<float> solidAngles(m_triangles.size());
    float totalSolidAngle = 0.0f;

    for (size_t i = 0; i < m_triangles.size(); ++i) {
        solidAngles[i] = m_triangles[i].areaImportanceSampler->calculateSolidAngle(shadingPoint);
        if (solidAngles[i] > 0.0f) {
            totalSolidAngle += solidAngles[i];
        }
    }

    if (totalSolidAngle <= 0.0f) {
        // fall back to uniform sampling if no triangle is visible
        return sampleUniform(shadingPoint, u1, u2, u3);
    }

    // 2. build solid angle CDF
    std::vector<float> solidAngleCDF(m_triangles.size());
    float cumulative = 0.0f;
    for (size_t i = 0; i < m_triangles.size(); ++i) {
        cumulative += solidAngles[i] / totalSolidAngle;
        solidAngleCDF[i] = cumulative;
    }
    solidAngleCDF.back() = 1.0f;

    // 3. select triangle by solid angle
    const auto it = std::lower_bound(solidAngleCDF.begin(), solidAngleCDF.end(), u1);
    size_t triIndex = std::distance(solidAngleCDF.begin(), it);
    triIndex = std::min(triIndex, m_triangles.size() - 1);

    // 4. sample using area importance on selected triangle
    const TriangleData& tri = m_triangles[triIndex];
    AreaLightSample sample = tri.areaImportanceSampler->sample(shadingPoint, u2, u3);

    // 5. adjust PDF for triangle selection
    float selectionProb = solidAngles[triIndex] / totalSolidAngle;
    sample.pdf *= selectionProb;
    sample.area = m_totalArea;

    return sample;
}

AreaLightSample MeshAreaLight::sampleHierarchicalFlux(const glm::vec3 &shadingPoint, float u1, float u2, float u3) const
{
    if (m_bvhNodes.empty() || m_rootNodeIndex == -1)
    {
        return sampleUniform(shadingPoint, u1, u2, u3); // fallback to uniform
    }


    // 1. select BVH node using flux based importance sampling
    // initialize path probability
    float pathPdf = 1.0f;
    int selectedNodeIndex = selectBVHNode(m_rootNodeIndex, u1, shadingPoint, pathPdf);
    if (selectedNodeIndex == -1) {
        return AreaLightSample{};
    }

    const BVHNode& node = m_bvhNodes[selectedNodeIndex];
    if (!node.isLeaf) {
        return sampleUniform(shadingPoint, u1, u2, u3); // fallback to uniform
    }

    // 2. within the leaf node, select a triangle
    int numTrianglesInLeaf = node.endTri - node.startTri;

    float triangleSelectU = u2;
    int triangleIndexInLeaf = static_cast<int>(triangleSelectU * static_cast<float>(numTrianglesInLeaf));
    triangleIndexInLeaf = std::min(triangleIndexInLeaf, numTrianglesInLeaf - 1);

    int actualTriangleIndex = node.startTri + triangleIndexInLeaf;
    const TriangleData& triangle = m_triangles[actualTriangleIndex];

    // 3. sample the selected triangle (uniform sampling for simplicity)
    AreaLightSample sample = triangle.uniformSampler->sample(shadingPoint, u3, std::fmod(u1 + u2, 1.0f));
    if (sample.pdf <= 0.0f) {
        return AreaLightSample{};
    }

    // 4: compute PDF
    // P(node) * P(triangle | node) * P(point | triangle)

    // probability of selecting this node = node.flux / root.flux
   //float rootFlux = m_bvhNodes[m_rootNodeIndex].totalFlux;
   //float nodeProbability = node.totalFlux / rootFlux;

    float triangleProbability = 1.0f / static_cast<float>(numTrianglesInLeaf);
    sample.pdf = pathPdf * sample.pdf * triangleProbability;

    return sample;
}

AreaLightSample MeshAreaLight::sample(const glm::vec3& shadingPoint,
                                      float u1, float u2, float u3) const {
    switch (m_strategy) {
        case SamplingStrategy::Uniform:
            return sampleUniform(shadingPoint, u1, u2, u3);
        case SamplingStrategy::AreaImportance:
            return sampleAreaImportance(shadingPoint, u1, u2, u3);
        case SamplingStrategy::HierarchicalFlux:
            return sampleHierarchicalFlux(shadingPoint, u1, u2, u3);
        default:
            return sampleUniform(shadingPoint, u1, u2, u3);
    }
}

float MeshAreaLight::pdfUniform(const glm::vec3& shadingPoint,
                               const glm::vec3& lightPoint) const {
    // find which triangle contains the point
    for (size_t i = 0; i < m_triangles.size(); ++i) {
        const TriangleData& tri = m_triangles[i];

        float u, v, w;
        if (isPointInTriangle(lightPoint, tri.v0, tri.v1, tri.v2, u, v, w)) {
            // point is on this triangle
            float trianglePdf = tri.uniformSampler->pdf(shadingPoint, lightPoint);
            float selectionProb = tri.area / m_totalArea;
            return trianglePdf * selectionProb;
        }
    }

    // point not on any triangle
    return 0.0f;
}

float MeshAreaLight::pdfAreaImportance(const glm::vec3& shadingPoint,
                                      const glm::vec3& lightPoint) const {
    // 1. compute solid angles for all triangles
    std::vector<float> solidAngles(m_triangles.size());
    float totalSolidAngle = 0.0f;

    for (size_t i = 0; i < m_triangles.size(); ++i) {
        solidAngles[i] = m_triangles[i].areaImportanceSampler->calculateSolidAngle(shadingPoint);
        if (solidAngles[i] > 0.0f) {
            totalSolidAngle += solidAngles[i];
        }
    }

    if (totalSolidAngle <= 0.0f) {
        return pdfUniform(shadingPoint, lightPoint);
    }

    // 2. find which triangle contains the point
    for (size_t i = 0; i < m_triangles.size(); ++i) {
        const TriangleData& tri = m_triangles[i];

        float u, v, w;
        if (isPointInTriangle(lightPoint, tri.v0, tri.v1, tri.v2, u, v, w)) {
            // Point is on this triangle
            float trianglePdf = tri.areaImportanceSampler->pdf(shadingPoint, lightPoint);
            float selectionProb = solidAngles[i] / totalSolidAngle;
            return trianglePdf * selectionProb;
        }
    }

    // point not on any triangle
    return 0.0f;
}

float MeshAreaLight::pdfHierarchicalFlux(const glm::vec3 &shadingPoint, const glm::vec3 &lightPoint) const {
    if (m_bvhNodes.empty()) {
        return pdfUniform(shadingPoint, lightPoint);
    }

    // find which triangle contains the point
    for (size_t i = 0; i < m_triangles.size(); ++i) {
        const TriangleData& triangle = m_triangles[i];

        float u, v, w;
        if (isPointInTriangle(lightPoint, triangle.v0, triangle.v1, triangle.v2, u, v, w)) {
            // found the triangle
            // now need to find which leaf node contains this triangle

            // simple linear search through leaf nodes
            // (potentially can be optimized by storing triangle-to-node mapping)
            for (const auto& node : m_bvhNodes) {
                if (!node.isLeaf) continue;

                // check if triangle i is in this leaf
                bool found = false;
                for (int j = node.startTri; j < node.endTri; ++j) {
                    if (static_cast<size_t>(j) == i) {
                        found = true;
                        break;
                    }
                }

                if (found) {
                    // compute PDF
                    const float rootFlux = m_bvhNodes[m_rootNodeIndex].totalFlux;
                    const float nodeProbability = node.totalFlux / rootFlux;

                    const int numTrianglesInLeaf = node.endTri - node.startTri;
                    const float triProbability = 1.0f / static_cast<float>(numTrianglesInLeaf);

                    const float trianglePdf = triangle.uniformSampler->pdf(shadingPoint, lightPoint);

                    return trianglePdf * nodeProbability * triProbability;
                }
            }

            // triangle found but not in BVH, fallback to uniform
            return pdfUniform(shadingPoint, lightPoint);
        }
    }

    // point not on any triangle
    return 0.0f;
}

float MeshAreaLight::pdf(const glm::vec3& shadingPoint,
                         const glm::vec3& lightPoint) const {
    switch (m_strategy) {
        case SamplingStrategy::Uniform:
            return pdfUniform(shadingPoint, lightPoint);
        case SamplingStrategy::AreaImportance:
            return pdfAreaImportance(shadingPoint, lightPoint);
        case SamplingStrategy::HierarchicalFlux:
            return pdfHierarchicalFlux(shadingPoint, lightPoint);
        default:
            return pdfUniform(shadingPoint, lightPoint);
    }
}

float MeshAreaLight::getTotalPower() const {
    float power = 0.0f;
    for (const auto& tri : m_triangles) {
        power += tri.uniformSampler->getTotalFlux();
    }
    return power;
}

void MeshAreaLight::setSamplingStrategy(SamplingStrategy strategy) {
    m_strategy = strategy;

    if ((strategy == SamplingStrategy::HierarchicalFlux ||
         strategy == SamplingStrategy::VisibilityAwareHierarchical) &&
        m_bvhNodes.empty()) {
        buildBVH();
    }
}


BVHNode::BVHNode() : isLeaf(false) {
}

BVHNode::BVHNode(int left, int right, float flux, float area, bool leaf) {
    isLeaf = leaf;
    if (isLeaf) {
        initializeLeafNode(left, right, flux, area);
    } else {
        initializeInternalNode(left, right, flux, area);
    }
}

void BVHNode::initializeInternalNode(int left, int right, float flux, float a) {
    leftChild = left;
    rightChild = right;
    startTri = 0;
    endTri = 0;
    totalFlux = flux;
    totalArea = a;
}

void BVHNode::initializeLeafNode(int start, int end, float flux, float a) {
    leftChild = -1;
    rightChild = -1;
    startTri = start;
    endTri = end;
    totalFlux = flux;
    totalArea = a;
}



