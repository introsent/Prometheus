//
// Created by minaj on 10/27/2025.
//
#ifndef RAY_TRACER_H
#define RAY_TRACER_H
#include <embree4/rtcore.h>
#include <glm/glm.hpp>
#include <vector>

#include "hit_result.h"
#include "ray.h"
#include "render/scene_manager.h"

class RayTracer {
public:
    explicit RayTracer(const SceneManager* sceneManager);

    // Single ray intersection
    HitResult intersect(const Ray& ray) const;

    // Packet intersection (up to 16 rays)
    void intersectPacket16(const std::vector<Ray>& rays, std::vector<HitResult>& results) const;

    // Occlusion test (faster than full intersection)
    bool isOccluded(const Ray& ray) const;
    void intersectionFilter(const struct RTCFilterFunctionNArguments* args) const;

private:
    glm::vec3 computeNormal(const HitResult& hit, const Ray& ray) const;

    RTCScene m_scene;
    const SceneManager* m_sceneManager;

    // Cached ray query context for reuse
    mutable RTCRayQueryContext m_rqctx{};
};



#endif //RAY_TRACER_H
