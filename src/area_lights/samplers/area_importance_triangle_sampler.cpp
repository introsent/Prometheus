//
// Created by ivans on 28/12/2025.
//

#include "area_importance_triangle_sampler.h"

#include <glm/ext/scalar_constants.hpp>

AreaImportanceTriangleSampler::AreaImportanceTriangleSampler(
    const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2,
    const glm::vec3& normal, float area,
    const glm::vec3& emission, float intensity)
    : m_v0(v0), m_v1(v1), m_v2(v2)
    , m_normal(normal)
    , m_area(area)
    , m_emission(emission)
    , m_intensity(intensity)
{}

AreaLightSample AreaImportanceTriangleSampler::sample(
    const glm::vec3& shadingPoint,
    float u1, float u2) const {

    // sample direction using spherical triangle projection
    float pdfSolidAngle = 0.0f;
    const glm::vec3 direction = sampleSphericalTriangle(
        shadingPoint, u1, u2, &pdfSolidAngle);

    if (pdfSolidAngle <= 0.0f) {
        return {shadingPoint, m_normal, 0.0f, 1.0f, glm::vec3(0), m_area};
    }

    // intersect ray with triangle plane to find exact sample position
    auto [t, position] = intersectTrianglePlane(shadingPoint, direction);

    if (position == glm::vec3(0.0f) || t <= 1e-6f) {
        return {shadingPoint, m_normal, 0.0f, 1.0f, glm::vec3(0), m_area};
    }

    // convert solid angle pdf to area pdf using jacobian
    // pdf_area = pdf_solidAngle * |cos(theta)| / r²
    //
    // where:
    //   theta = angle between light normal and direction to shading point
    //   r = distance from shading point to light sample
    const float cosTheta = glm::dot(m_normal, -direction);

    if (cosTheta <= 0.0f) {
        return {position, m_normal, 0.0f, 1.0f, glm::vec3(0), m_area};
    }

    const float pdfArea = pdfSolidAngle * (cosTheta / (t * t));

    if (!std::isfinite(pdfArea) || pdfArea <= 0.0f) {
        return {position, m_normal, 0.0f, 1.0f, glm::vec3(0), m_area};
    }

    return {
        position,
        m_normal,
        pdfArea,
        1.0f,
        m_emission * m_intensity,
        m_area
    };
}

float AreaImportanceTriangleSampler::pdf(
    const glm::vec3& shadingPoint,
    const glm::vec3& lightPoint) const {

    // compute solid angle subtended by triangle
    const float solidAngle = calculateSolidAngle(shadingPoint);
    if (solidAngle <= 0.0f) return 0.0f;

    // direction from light sample to shading point
    const glm::vec3 toShading = shadingPoint - lightPoint;
    const float distSq = glm::dot(toShading, toShading);
    if (distSq <= 0.0f) return 0.0f;

    const float dist = std::sqrt(distSq);
    const glm::vec3 direction = toShading / dist;

    // cosine term at light surface
    const float cosTheta = glm::dot(m_normal, direction);
    if (cosTheta <= 0.0f) return 0.0f;

    // uniform pdf over solid angle: 1/omega
    const float pdfSolidAngle = 1.0f / solidAngle;

    // convert to area pdf
    return pdfSolidAngle * (cosTheta / distSq);
}

void AreaImportanceTriangleSampler::setIntensity(float intensity) {
    m_intensity = intensity;
}

float AreaImportanceTriangleSampler::getTotalFlux() const {
    return MathHelpers::computeTriangleFlux(m_emission, m_intensity, m_area);
}

float AreaImportanceTriangleSampler::calculateSolidAngle(const glm::vec3& p) const {
    // van Oosterom & Strackee formula for solid angle of triangle
    // ref: "The Solid Angle of a Plane Triangle" (1983)
    // omega = 2 * arctan( |a dot (b cross c)| / (1 + a dot b + b dot c + c dot a) )
    // where a, b, c are unit vectors from point p to triangle vertices

    // vectors from point to vertices
    const glm::vec3 a = m_v0 - p;
    const glm::vec3 b = m_v1 - p;
    const glm::vec3 c = m_v2 - p;

    // lengths
    const float la = glm::length(a);
    const float lb = glm::length(b);
    const float lc = glm::length(c);

    // degenerate case: point on vertex
    if (la == 0.0f || lb == 0.0f || lc == 0.0f) {
        return 0.0f;
    }

    // normalize to unit vectors
    const glm::vec3 aUnit = a / la;
    const glm::vec3 bUnit = b / lb;
    const glm::vec3 cUnit = c / lc;

    // scalar triple product: a dot (b cross c)
    const float tripleProduct = glm::dot(aUnit, glm::cross(bUnit, cUnit));

    // dot products between unit vectors
    const float dotAB = glm::dot(aUnit, bUnit);
    const float dotBC = glm::dot(bUnit, cUnit);
    const float dotCA = glm::dot(cUnit, aUnit);

    // denominator
    const float denominator = 1.0f + dotAB + dotBC + dotCA;

    // check for back-facing or degenerate triangle
    if (denominator <= 0.0f || std::abs(tripleProduct) < 1e-10f) {
        const float signedVolume = glm::dot(a, glm::cross(b, c));
        if (signedVolume <= 0.0f) {
            return 0.0f;  // back-facing
        }
        return 1e-10f;  // degenerate but front-facing
    }

    // compute solid angle using atan2 for numerical stability
    const float solidAngle = 2.0f * std::atan2(
        std::abs(tripleProduct),
        denominator
    );

    // clamp to valid range [0, 2pi]
    return std::max(0.0f, std::min(solidAngle, 2.0f * glm::pi<float>()));
}

glm::vec3 AreaImportanceTriangleSampler::sampleSphericalTriangle(
    const glm::vec3& p,
    float u1, float u2,
    float* outPdfSolidAngle) const {

    // Arvo's method for sampling spherical triangles
    // ref: "Stratified Sampling of Spherical Triangles" (1995)

    // project triangle vertices onto unit sphere centered at p
    glm::vec3 a = m_v0 - p;
    glm::vec3 b = m_v1 - p;
    glm::vec3 c = m_v2 - p;

    const float la = glm::length(a);
    const float lb = glm::length(b);
    const float lc = glm::length(c);

    if (la <= 0.0f || lb <= 0.0f || lc <= 0.0f) {
        if (outPdfSolidAngle) *outPdfSolidAngle = 0.0f;
        return glm::normalize(m_v0 - p);
    }

    a /= la; b /= lb; c /= lc;

    // normals to great circle arcs (edges of spherical triangle)
    glm::vec3 nAB = glm::cross(a, b);
    glm::vec3 nBC = glm::cross(b, c);
    glm::vec3 nCA = glm::cross(c, a);

    constexpr float eps = 1e-10f;
    if (glm::dot(nAB, nAB) <= eps ||
        glm::dot(nBC, nBC) <= eps ||
        glm::dot(nCA, nCA) <= eps) {
        if (outPdfSolidAngle) *outPdfSolidAngle = 0.0f;
        return glm::normalize(a + b + c);
    }

    nAB = glm::normalize(nAB);
    nBC = glm::normalize(nBC);
    nCA = glm::normalize(nCA);

    // spherical angles at vertices (dihedral angles between planes)
    const float alpha = MathHelpers::safeAcos(glm::dot(nAB, -nCA));
    const float beta = MathHelpers::safeAcos(glm::dot(nBC, -nAB));
    const float gamma = MathHelpers::safeAcos(glm::dot(nCA, -nBC));

    // spherical excess gives solid angle
    // omega = alpha + beta + gamma - pi (Girard's theorem)
    const float sumAngles = alpha + beta + gamma;
    const float solidAngle = sumAngles - glm::pi<float>();

    if (solidAngle <= 0.0f) {
        if (outPdfSolidAngle) *outPdfSolidAngle = 0.0f;
        return glm::normalize(a + b + c);
    }

    if (outPdfSolidAngle) {
        *outPdfSolidAngle = 1.0f / solidAngle;
    }

    // sample using Arvo's parameterization
    // map u1 to adjusted angle in [pi, alpha + beta + gamma]
    const float Ap = glm::pi<float>() + u1 * (sumAngles - glm::pi<float>());

    // trigonometric calculations for spherical law of cosines
    const float cosAlpha = std::cos(alpha);
    const float sinAlpha = std::sin(alpha);
    const float sinAp = std::sin(Ap);
    const float cosAp = std::cos(Ap);

    const float sinPhi = sinAp * cosAlpha - cosAp * sinAlpha;
    const float cosPhi = cosAp * cosAlpha + sinAp * sinAlpha;

    const float cosC = glm::dot(a, b);

    const float k1 = cosPhi + cosAlpha;
    const float k2 = sinPhi - sinAlpha * cosC;
    const float denominator = (k2 * sinPhi + k1 * cosPhi) * sinAlpha;
    const float numerator = k2 + (k2 * cosPhi - k1 * sinPhi) * cosAlpha;

    // compute cos(B') with fallback for degenerate cases
    float cosBp;
    if (std::abs(denominator) < 1e-10f) {
        cosBp = std::clamp(
            cosC * (1.0f - u1) + glm::dot(c, b) * u1,
            -1.0f, 1.0f
        );
    } else {
        cosBp = std::clamp(numerator / denominator, -1.0f, 1.0f);
    }

    const float sinBp = MathHelpers::safeSqrt(1.0f - cosBp * cosBp);

    // construct point on arc AC
    const glm::vec3 axisAC = MathHelpers::gramSchmidtNormalize(c, a);
    glm::vec3 cp = cosBp * a + sinBp * axisAC;
    cp = glm::normalize(cp);

    // rotate around B to final position
    const float dotCpB = glm::dot(cp, b);
    const float cosTheta = std::clamp(
        1.0f - u2 * (1.0f - dotCpB),
        -1.0f, 1.0f
    );
    const float sinTheta = MathHelpers::safeSqrt(1.0f - cosTheta * cosTheta);

    const glm::vec3 axisBp = MathHelpers::gramSchmidtNormalize(cp, b);
    glm::vec3 w = cosTheta * b + sinTheta * axisBp;

    return glm::normalize(w);
}

std::pair<float, glm::vec3> AreaImportanceTriangleSampler::intersectTrianglePlane(
    const glm::vec3& rayOrigin,
    const glm::vec3& rayDirection) const {

    // ray-plane intersection: solve for t in P = O + t*D
    // where plane is defined by point v0 and normal n
    const float denominator = glm::dot(m_normal, rayDirection);

    if (std::abs(denominator) < 1e-6f) {
        return {0.0f, glm::vec3(0.0f)};  // parallel
    }

    const float t = glm::dot(m_v0 - rayOrigin, m_normal) / denominator;

    if (t <= 0.0f) {
        return {0.0f, glm::vec3(0.0f)};  // behind ray origin
    }

    const glm::vec3 position = rayOrigin + t * rayDirection;
    return {t, position};
}