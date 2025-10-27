//
// Created by minaj on 10/4/2025.
//

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

#include "init/embree_scene.h"
#include "embree4/rtcore.h"

struct Vertex;

class Geometry {
public:
    explicit Geometry(RTCGeometryType geometryType, RTCDevice device);

    [[nodiscard]] RTCGeometryType getType() const;
    [[nodiscard]] RTCGeometry getGeometry() const;

    ~Geometry() = default;

    void commit() const;
    unsigned attach(const EmbreeScene* embreeScenePtr) const;
    void release() const;
protected:
    RTCGeometryType m_type;
    RTCGeometry m_geometry = nullptr;
};



#endif //GEOMETRY_H
