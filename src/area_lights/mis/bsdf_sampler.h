//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_BSDF_SAMPLER_H
#define PROMETHEUS_BSDF_SAMPLER_H

#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

/// BSDF sample result
struct BSDFSample {
    glm::vec3 direction;     // sampled direction (in world space)
    float pdf;               // probability density function
    glm::vec3 reflectance;   // BRDF value (f_r / pdf for importance sampling)
};

/// BSDF sampler
// generates importance-sampled directions for surface reflection
class BSDFSampler {
public:
    // sample diffuse BRDF using cosine-weighted hemisphere sampling
    // pdf = cos(theta) / pi
    // BRDF = ro / pi  (where ro is albedo, assumed 1.0 here)

    // ref: "Monte Carlo Techniques for Direct Lighting Calculations" (Shirley et al., 1996)
    static BSDFSample sampleDiffuse(
        const glm::vec3& normal,
        float u1, float u2) {
        u1 = std::max(1e-6f, u1);
        u1 = std::min(u1, 1.0f - 1e-6f);


        BSDFSample sample{};

        // cosine-weighted hemisphere sampling
        // maps uniform [0,1]^2 to directions with density proportional to cos(θ)
        const glm::vec3 localDir = cosineSampleHemisphere(u1, u2);

        // transform to world space aligned with surface normal
        sample.direction = alignWithNormal(localDir, normal);

        // pdf for cosine-weighted sampling
        const float cosTheta = std::max(0.0f, glm::dot(sample.direction, normal));
        sample.pdf = cosTheta / glm::pi<float>();

        // diffuse BRDF (albedo = 1.0)
        sample.reflectance = glm::vec3(1.0f / glm::pi<float>());

        return sample;
    }

    // evaluate pdf for diffuse BRDF
    static float pdfDiffuse(const glm::vec3& normal, const glm::vec3& direction) {
        float cosTheta = std::max(0.0f, glm::dot(direction, normal));
        return cosTheta / glm::pi<float>();
    }

    // Sample GGX distribution for specular materials
    static BSDFSample sampleGGX(
        const glm::vec3& normal,
        const glm::vec3& viewDir,
        float roughness,
        float u1, float u2) {

        BSDFSample sample{};

        float alpha = roughness * roughness;

        // sample half-vector from GGX distribution
        float phi = 2.0f * glm::pi<float>() * u1;
        float cosTheta = std::sqrt((1.0f - u2) / (1.0f + (alpha * alpha - 1.0f) * u2));
        float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);

        // half-vector in local space
        glm::vec3 hLocal(
            sinTheta * std::cos(phi),
            sinTheta * std::sin(phi),
            cosTheta
        );

        // transform to world space
        glm::vec3 h = alignWithNormal(hLocal, normal);

        // reflect view direction around half-vector to get light direction
        sample.direction = glm::reflect(-viewDir, h);

        // check if direction is valid (above surface)
        float NdotL = glm::dot(normal, sample.direction);
        if (NdotL <= 0.0f) {
            sample.pdf = 0.0f;
            return sample;
        }

        float NdotH = std::max(glm::dot(normal, h), 0.0f);
        float VdotH = std::max(glm::dot(viewDir, h), 0.0f);

        // GGX PDF: D(h) * NdotH / (4 * VdotH)
        float D = normalDistributionGGX(normal, h, roughness);
        sample.pdf = (D * NdotH) / (4.0f * VdotH + 1e-6f);

        return sample;
    }

    static float pdfGGX(
        const glm::vec3& normal,
        const glm::vec3& viewDir,
        const glm::vec3& lightDir,
        float roughness) {

        glm::vec3 h = glm::normalize(viewDir + lightDir);

        float NdotH = std::max(glm::dot(normal, h), 0.0f);
        float VdotH = std::max(glm::dot(viewDir, h), 0.0f);

        float D = normalDistributionGGX(normal, h, roughness);
        return (D * NdotH) / (4.0f * VdotH + 1e-6f);
    }
private:
    static float normalDistributionGGX(const glm::vec3& n, const glm::vec3& h, float roughness) {
        const float alpha = roughness * roughness;
        const float alphaSq = alpha * alpha;
        const float nDotH = glm::dot(n, h);
        const float denom = nDotH * nDotH * (alphaSq - 1.0f) + 1.0f;
        return alphaSq / (glm::pi<float>() * denom * denom);
    }

    // cosine-weighted hemisphere sampling
    // returns direction in local space (z+ is up)
    static glm::vec3 cosineSampleHemisphere(float u1, float u2) {
        float z = std::sqrt(u1);
        const float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
        const float phi = 2.0f * glm::pi<float>() * u2;
        return {r * std::cos(phi), r * std::sin(phi), z};
    }

    // align local direction with surface normal
    static glm::vec3 alignWithNormal(
        const glm::vec3& localDir,
        const glm::vec3& normal) {

        // construct orthonormal basis around normal
        const glm::vec3 up = std::abs(normal.y) < 0.999f ?
            glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);

        const glm::vec3 tangent = glm::normalize(glm::cross(up, normal));
        const glm::vec3 bitangent = glm::cross(normal, tangent);

        // transform local direction to world space
        return tangent * localDir.x +
               bitangent * localDir.y +
               normal * localDir.z;
    }
};

#endif //PROMETHEUS_BSDF_SAMPLER_H