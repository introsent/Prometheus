//
// Created by minaj on 10/27/2025.
//

#ifndef PLANE_H
#define PLANE_H

#include "../init/embree_device.h"
#include "geometry.h"
#include "glm/vec3.hpp"

class Plane : public Geometry {
public:
    Plane(const glm::vec3& origin, const glm::vec3& normal, const EmbreeDevice* devicePtr);

    void setMaterialId(unsigned char materialId) { m_materialId = materialId; }
    [[nodiscard]] unsigned char getMaterialId() const { return m_materialId; }

private:
    unsigned char m_materialId = 0;
};

#endif //PLANE_H
