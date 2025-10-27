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
        result.origin = ray.origin + ray.direction * result.distance;
    }

    return result;
}

void RayTracer::intersectPacket16(const std::vector<Ray> &rays, std::vector<HitResult> &results) const {
    const size_t N = rays.size();
    assert(N <= 16); // this helper expects at most 16 rays per call

    // Prepare results
    results.clear();
    results.resize(N);

    // Create and zero the SOA rayhit structure with required alignment
    alignas(64) RTCRayHit16 rayhit = {};

    // valid mask: -1 for active lane, 0 for inactive
    int valid[16];
    for (int & i : valid) i = 0;

    // Fill per-lane ray data (SOA layout)
    for (size_t i = 0; i < N; ++i) {
        const Ray& r = rays[i];

        rayhit.ray.org_x[i] = r.origin.x;
        rayhit.ray.org_y[i] = r.origin.y;
        rayhit.ray.org_z[i] = r.origin.z;

        rayhit.ray.dir_x[i] = r.direction.x;
        rayhit.ray.dir_y[i] = r.direction.y;
        rayhit.ray.dir_z[i] = r.direction.z;

        rayhit.ray.tnear[i] = r.tNear;
        rayhit.ray.tfar[i] = r.tFar;

        rayhit.ray.time[i] = 0.0f; // motion blur time if you use it
        rayhit.ray.flags[i] = 0;
        rayhit.ray.mask[i] = r.mask;

        // initialize hits to invalid
        rayhit.hit.geomID[i] = RTC_INVALID_GEOMETRY_ID;

        valid[i] = -1; // mark active
    }

    RTCRayQueryContext rqctx{};
    rtcInitRayQueryContext(&rqctx);

    RTCIntersectArguments args{};
    rtcInitIntersectArguments(&args);
    args.context = &rqctx;

    // Call Embree packet intersect
    rtcIntersect16(valid, m_scene, &rayhit, &args);

    // Gather results per lane
    for (size_t i = 0; i < N; ++i) {
        HitResult hr;
        if (rayhit.hit.geomID[i] != RTC_INVALID_GEOMETRY_ID) {
            hr.didHit = true;
            hr.distance = rayhit.ray.tfar[i];
            hr.geomID = rayhit.hit.geomID[i];
            hr.primID = rayhit.hit.primID[i];
            hr.u = rayhit.hit.u[i];
            hr.v = rayhit.hit.v[i];
            hr.origin = rays[i].origin + rays[i].direction * hr.distance;
        } else {
            hr.didHit = false;
        }
        results[i] = hr;
    }
}

bool RayTracer::isOccluded(const Ray& ray) const {
    RTCRay rtcRay = {};
    rtcRay.org_x = ray.origin.x;
    rtcRay.org_y = ray.origin.y;
    rtcRay.org_z = ray.origin.z;
    rtcRay.dir_x = ray.direction.x;
    rtcRay.dir_y = ray.direction.y;
    rtcRay.dir_z = ray.direction.z;
    rtcRay.tnear = ray.tNear;
    rtcRay.tfar = ray.tFar;
    rtcRay.mask = ray.mask;
    rtcRay.flags = 0;

    RTCRayQueryContext rqctx{};
    rtcInitRayQueryContext(&rqctx);

    RTCOccludedArguments args{};
    rtcInitOccludedArguments(&args);
    args.context = &rqctx;

    rtcOccluded1(m_scene, &rtcRay, &args);

    return rtcRay.tfar < 0.0f;
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
