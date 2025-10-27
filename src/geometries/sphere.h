//
// Created by minaj on 10/27/2025.
//

#ifndef SPHERE_H
#define SPHERE_H
#include "geometry.h"
#include "glm/vec3.hpp"


class EmbreeDevice;

class Sphere : public Geometry {
public:
    Sphere(const glm::vec3& center, float radius, const EmbreeDevice* devicePtr);

    void setMaterialId(unsigned char materialId) { m_materialId = materialId; }
    [[nodiscard]] unsigned char getMaterialId() const { return m_materialId; }

private:
    unsigned char m_materialId = 0;
};



#endif //SPHERE_H
