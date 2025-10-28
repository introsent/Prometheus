//
// Created by minaj on 10/27/2025.
//
#ifndef RAY_TRACER_H
#define RAY_TRACER_H

#include "ray.h"
#include "hit_result.h"
#include <embree4/rtcore.h>
#include <vector>
#include <unordered_map>
#include "render/scene_manager.h"

class RayTracer {
public:
    explicit RayTracer(const SceneManager* sceneManager);

    HitResult intersect(const Ray& ray) const;
    void intersectPacket16(const std::vector<Ray>& rays, std::vector<HitResult>& results) const;
    bool isOccluded(const Ray& ray) const;
    glm::vec3 computeNormal(const HitResult& hit, const Ray& ray) const;

private:
    RTCScene m_scene;
    const SceneManager* m_sceneManager;
    RTCRayQueryContext m_rqctx{};

    // NEW: Cached geometry data using unordered_map for O(1) lookup
    struct GeometryCache {
        const float* vertices;
        const unsigned* indices;
        RTCGeometryType type;
    };
    mutable std::unordered_map<unsigned, GeometryCache> m_geometryCache;

    // Cache geometry on first access
    const GeometryCache& getCachedGeometry(unsigned geomID) const;
};


#endif //RAY_TRACER_H
