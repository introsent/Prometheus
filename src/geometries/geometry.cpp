//
// Created by minaj on 10/4/2025.
//
#include "geometry.h"

Geometry::Geometry(RTCGeometryType geometryType, RTCDevice device) : m_type(geometryType) {
    m_geometry = rtcNewGeometry(device, m_type);
}

RTCGeometryType Geometry::getType() const {
    return m_type;
}

RTCGeometry Geometry::getGeometry() const {
    return m_geometry;
}

void Geometry::commit() const {
    rtcCommitGeometry(m_geometry);
}

unsigned Geometry::attach(const EmbreeScene* embreeScenePtr) const {
    return rtcAttachGeometry(embreeScenePtr->handle(), m_geometry);
}

void Geometry::release() const {
    rtcReleaseGeometry(m_geometry);
}

