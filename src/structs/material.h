//
// Created by minaj on 10/27/2025.
//

#ifndef MATERIAL_H
#define MATERIAL_H

#include "glm/vec3.hpp"

enum class MaterialType {
    SolidColor
};

struct Material {
    virtual ~Material() = default;
    [[nodiscard]] virtual MaterialType getType() const = 0;
    [[nodiscard]] virtual glm::vec3 getColor() const = 0;
};

struct Material_SolidColor final : public Material {
    glm::vec3 color;

    explicit Material_SolidColor(const glm::vec3& color) : color(color) {}

    [[nodiscard]] MaterialType getType() const override { return MaterialType::SolidColor; }
    [[nodiscard]] glm::vec3 getColor() const override { return color; }
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
