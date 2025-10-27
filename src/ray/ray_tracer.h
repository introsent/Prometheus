//
// Created by minaj on 10/27/2025.
//

#ifndef RAY_TRACER_H
#define RAY_TRACER_H
#include <vector>

#include "hit_result.h"
#include "ray.h"
#include "embree4/rtcore_device.h"
#include "embree4/rtcore_ray.h"


class RayTracer {
public:
    explicit RayTracer(RTCScene scene);
    [[nodiscard]] HitResult intersect(const Ray& ray) const;

    void intersectPacket16(const std::vector<Ray>& rays, std::vector<HitResult>& results) const;

    [[nodiscard]] bool isOccluded(const Ray& ray) const;

private:
    RTCScene m_scene;

    static void setupRay(RTCRayHit& rayHit, const Ray& ray);
};



#endif //RAY_TRACER_H
