//
// Created by minaj on 10/27/2025.
//

#include "ray_tracer.h"

#include <embree4/rtcore.h>

RayTracer::RayTracer(RTCScene scene) : m_scene(scene)
{
}

HitResult RayTracer::intersect(const Ray &ray) const {
    HitResult result;

    RTCRayHit rayHit = {};
    setupRay(rayHit, ray);

    RTCRayQueryContext rqctx{};
    rtcInitRayQueryContext(&rqctx);

    RTCIntersectArguments args{};
    rtcInitIntersectArguments(&args);
    args.context = &rqctx;

    rtcIntersect1(m_scene, &rayHit, &args);

    if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
        result.didHit = true;
        result.distance = rayHit.ray.tfar;
        result.geomID = rayHit.hit.geomID;
        result.primID = rayHit.hit.primID;
        result.u = rayHit.hit.u;
        result.v = rayHit.hit.v;
    }

    return result;
}

void RayTracer::setupRay(RTCRayHit &rayHit, const Ray &ray) {
    rayHit.ray.org_x = ray.origin.x;
    rayHit.ray.org_y = ray.origin.y;
    rayHit.ray.org_z = ray.origin.z;
    rayHit.ray.dir_x = ray.direction.x;
    rayHit.ray.dir_y = ray.direction.y;
    rayHit.ray.dir_z = ray.direction.z;
    rayHit.ray.tnear = ray.tNear;
    rayHit.ray.tfar = ray.tFar;
    rayHit.ray.flags = 0;
    rayHit.ray.time = 0.0f;
    rayHit.ray.mask = ray.mask;
    rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
}
