//
// Created by minaj on 10/27/2025.
//

#ifndef MATERIAL_H
#define MATERIAL_H

#include "glm/vec3.hpp"
#include <algorithm>
#include <cmath>

#include "glm/detail/func_geometric.inl"

#ifndef PI
#define PI 3.14159265358979323846f
#endif

enum class MaterialType {
    SolidColor,
    Lambert,
    LambertPhong,
    CookTorrence
};

// BRDF Helper Functions
namespace BRDF {
    inline glm::vec3 Lambert(float kd, const glm::vec3& cd) {
        return kd * cd / PI;
    }

    inline glm::vec3 Lambert(const glm::vec3& kd, const glm::vec3& cd) {
        return kd * cd / PI;
    }

    inline glm::vec3 Phong(float ks, float exp, const glm::vec3& l, const glm::vec3& v, const glm::vec3& n) {
        const glm::vec3 reflect = glm::reflect(l, n);
        const float cosAngle = std::abs(glm::dot(reflect, v));
        const float specular = ks * std::pow(cosAngle, exp);
        return glm::vec3(specular);
    }

    inline glm::vec3 FresnelFunction_Schlick(const glm::vec3& h, const glm::vec3& v, const glm::vec3& f0) {
        return f0 + (glm::vec3(1.f) - f0) * std::pow(1.f - glm::dot(h, v), 5.f);
    }

    inline float NormalDistribution_GGX(const glm::vec3& n, const glm::vec3& h, float roughness) {
        const float alpha = roughness * roughness;
        const float alphaSq = alpha * alpha;
        const float nDotH = glm::dot(n, h);
        const float denom = nDotH * nDotH * (alphaSq - 1.f) + 1.f;
        return alphaSq / (PI * denom * denom);
    }

    inline float GeometryFunction_SchlickGGX(const glm::vec3& n, const glm::vec3& v, float roughness) {
        float k = ((roughness + 1.f) * (roughness + 1.f)) / 8.f;
        float nDotV = glm::dot(n, v);
        return nDotV / (nDotV * (1.f - k) + k);
    }

    inline float GeometryFunction_Smith(const glm::vec3& n, const glm::vec3& v, const glm::vec3& l, float roughness) {
        return GeometryFunction_SchlickGGX(n, v, roughness) * GeometryFunction_SchlickGGX(n, l, roughness);
    }
}


struct Material {
    virtual ~Material() = default;
    [[nodiscard]] virtual MaterialType getType() const = 0;
    [[nodiscard]] virtual glm::vec3 shade(const glm::vec3& hitPoint, const glm::vec3& normal,
                                          const glm::vec3& viewDir, const glm::vec3& lightDir) const = 0;
};

struct Material_SolidColor final : public Material {
    glm::vec3 color;

    explicit Material_SolidColor(const glm::vec3& color) : color(color) {}

    [[nodiscard]] MaterialType getType() const override { return MaterialType::SolidColor; }

    [[nodiscard]] glm::vec3 shade(const glm::vec3& hitPoint, const glm::vec3& normal,
                                  const glm::vec3& viewDir, const glm::vec3& lightDir) const override {
        return color;
    }
};

struct Material_Lambert final : public Material {
    glm::vec3 diffuseColor;
    float diffuseReflectance;

    Material_Lambert(const glm::vec3& diffuseColor, float kd)
        : diffuseColor(diffuseColor), diffuseReflectance(kd) {}

    [[nodiscard]] MaterialType getType() const override { return MaterialType::Lambert; }

    [[nodiscard]] glm::vec3 shade(const glm::vec3& hitPoint, const glm::vec3& normal,
                                  const glm::vec3& viewDir, const glm::vec3& lightDir) const override {
        return BRDF::Lambert(diffuseReflectance, diffuseColor);
    }
};

struct Material_LambertPhong final : public Material {
    glm::vec3 diffuseColor;
    float diffuseReflectance;  // kd
    float specularReflectance; // ks
    float phongExponent;

    Material_LambertPhong(const glm::vec3& diffuseColor, const float kd, const float ks, const float exp)
        : diffuseColor(diffuseColor), diffuseReflectance(kd),
          specularReflectance(ks), phongExponent(exp) {}

        [[nodiscard]] MaterialType getType() const override { return MaterialType::LambertPhong; }

        [[nodiscard]] glm::vec3 shade(const glm::vec3& hitPoint, const glm::vec3& normal,
                                  const glm::vec3& viewDir, const glm::vec3& lightDir) const override {
            const glm::vec3 diffuse = BRDF::Lambert(diffuseReflectance, diffuseColor);
            const glm::vec3 specular = BRDF::Phong(specularReflectance, phongExponent, lightDir, -viewDir, normal);
            return diffuse + specular;
        }
};

struct Material_CookTorrence final : public Material {
    glm::vec3 albedo;
    float metalness;
    float roughness;

    Material_CookTorrence(const glm::vec3& albedo, const float metalness, const float roughness)
        : albedo(albedo), metalness(metalness), roughness(roughness) {}

        [[nodiscard]] MaterialType getType() const override { return MaterialType::CookTorrence; }

        [[nodiscard]] glm::vec3 shade(const glm::vec3& hitPoint, const glm::vec3& normal,
                                  const glm::vec3& viewDir, const glm::vec3& lightDir) const override {
        const glm::vec3 f0 = (metalness == 0.f) ? glm::vec3(0.04f) : albedo;

        const glm::vec3 h = glm::normalize(viewDir + lightDir);

        const glm::vec3 F = BRDF::FresnelFunction_Schlick(h, viewDir, f0);
        const float D = BRDF::NormalDistribution_GGX(normal, h, roughness);
        const float G = BRDF::GeometryFunction_Smith(normal, viewDir, lightDir, roughness);

        const float nDotV = std::max(glm::dot(viewDir, normal), 0.001f);
        const float nDotL = std::max(glm::dot(lightDir, normal), 0.001f);

        const glm::vec3 specular = (D * F * G) / (4.f * nDotV * nDotL);

        const glm::vec3 kd = (metalness == 1.f) ? glm::vec3(0.f) : (glm::vec3(1.f) - F);

        return BRDF::Lambert(kd, albedo) + specular;
    }
};

namespace colors {
    constexpr glm::vec3 red{1.f, 0.f, 0.f};
    constexpr glm::vec3 blue{0.f, 0.f, 1.f};
    constexpr glm::vec3 yellow{1.f, 1.f, 0.f};
    constexpr glm::vec3 green{0.f, 1.f, 0.f};
    constexpr glm::vec3 magenta{1.f, 0.f, 1.f};
    constexpr glm::vec3 white{1.f, 1.f, 1.f};
    constexpr glm::vec3 black{0.f, 0.f, 0.f};
}


#endif //MATERIAL_H
