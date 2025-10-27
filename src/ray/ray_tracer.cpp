//
// Created by minaj on 10/27/2025.
// Optimized for maximum performance
//

#include "ray_tracer.h"
#include <embree4/rtcore.h>

RayTracer::RayTracer(const SceneManager* sceneManager)
    : m_scene(sceneManager->getScene())
    , m_sceneManager(sceneManager)
{
    // Pre-initialize reusable context
    rtcInitRayQueryContext(&m_rqctx);
}

HitResult RayTracer::intersect(const Ray &ray) const {
    HitResult result;

    RTCRayHit rayHit{
        .ray = {
            .org_x = ray.origin.x,
            .org_y = ray.origin.y,
            .org_z = ray.origin.z,
            .tnear = ray.tNear,
            .dir_x = ray.direction.x,
            .dir_y = ray.direction.y,
            .dir_z = ray.direction.z,
            .time = 0.0f,
            .tfar = ray.tFar,
            .mask = ray.mask,
            .id = 0,
            .flags = 0
        },
        .hit = {
            .geomID = RTC_INVALID_GEOMETRY_ID
        }
    };

    RTCIntersectArguments args{};
    rtcInitIntersectArguments(&args);
    args.context = const_cast<RTCRayQueryContext*>(&m_rqctx);

    rtcIntersect1(m_scene, &rayHit, &args);

    if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
        return result;
    }

    // Check culling BEFORE computing full normal
    CullingMode culling = m_sceneManager->getCullingMode(rayHit.hit.geomID);
    if (culling != CullingMode::NONE) {
        // Fast culling check, compute unnormalized normal
        RTCGeometry geom = rtcGetGeometry(m_scene, rayHit.hit.geomID);
        const RTCGeometryType geomType = m_sceneManager->getGeometryType(rayHit.hit.geomID);

        if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
            const auto* vertices = static_cast<const float*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
            const auto* indices = static_cast<const unsigned*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

            const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
            const unsigned i0 = indices[rayHit.hit.primID * stride + 0];
            const unsigned i1 = indices[rayHit.hit.primID * stride + 1];
            const unsigned i2 = indices[rayHit.hit.primID * stride + 2];

            const float* v0 = vertices + i0 * 3;
            const float* v1 = vertices + i1 * 3;
            const float* v2 = vertices + i2 * 3;

            // Compute cross product (unnormalized is fine for sign check)
            const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
            const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
            const glm::vec3 normal = glm::cross(edge1, edge2);

            const float dot = glm::dot(normal, ray.direction);

            // Cull based on face orientation
            if ((culling == CullingMode::BACK_FACE && dot > 0.0f) ||
                (culling == CullingMode::FRONT_FACE && dot < 0.0f)) {
                return result; // Return no hit (culled)
            }
        }
    }

    // Hit passed culling test - fill result
    result.didHit = true;
    result.distance = rayHit.ray.tfar;
    result.geomID = rayHit.hit.geomID;
    result.primID = rayHit.hit.primID;
    result.u = rayHit.hit.u;
    result.v = rayHit.hit.v;
    result.origin = ray.origin + ray.direction * result.distance;
    result.normal = computeNormal(result, ray);

    return result;
}

void RayTracer::intersectPacket16(const std::vector<Ray> &rays, std::vector<HitResult> &results) const {
    const size_t N = rays.size();
    assert(N <= 16);

    results.resize(N);

    // Aligned SOA structure
    alignas(64) RTCRayHit16 rayhit{};

    // Stack-allocate valid mask
    alignas(64) int valid[16];

    // Unrolled initialization and data setup
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

        rayhit.ray.time[i] = 0.0f;
        rayhit.ray.flags[i] = 0;
        rayhit.ray.mask[i] = r.mask;

        rayhit.hit.geomID[i] = RTC_INVALID_GEOMETRY_ID;

        valid[i] = -1;
    }

    // Zero out inactive lanes
    for (size_t i = N; i < 16; ++i) {
        valid[i] = 0;
    }

    RTCIntersectArguments args{};
    rtcInitIntersectArguments(&args);
    args.context = &m_rqctx;

    rtcIntersect16(valid, m_scene, &rayhit, &args);

    // Gather results with culling check
    for (size_t i = 0; i < N; ++i) {
        HitResult& hr = results[i];
        const bool hit = (rayhit.hit.geomID[i] != RTC_INVALID_GEOMETRY_ID);

        hr.didHit = false;

        if (hit) {
            // Check culling BEFORE computing full results
            CullingMode culling = m_sceneManager->getCullingMode(rayhit.hit.geomID[i]);
            bool culled = false;

            if (culling != CullingMode::NONE) {
                RTCGeometry geom = rtcGetGeometry(m_scene, rayhit.hit.geomID[i]);
                const RTCGeometryType geomType = m_sceneManager->getGeometryType(rayhit.hit.geomID[i]);

                if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
                    const auto* vertices = static_cast<const float*>(
                        rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
                    const auto* indices = static_cast<const unsigned*>(
                        rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

                    const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
                    const unsigned primID = rayhit.hit.primID[i];
                    const unsigned i0 = indices[primID * stride + 0];
                    const unsigned i1 = indices[primID * stride + 1];
                    const unsigned i2 = indices[primID * stride + 2];

                    const float* v0 = vertices + i0 * 3;
                    const float* v1 = vertices + i1 * 3;
                    const float* v2 = vertices + i2 * 3;

                    // Compute unnormalized normal
                    const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
                    const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
                    const glm::vec3 normal = glm::cross(edge1, edge2);

                    const float dot = glm::dot(normal, rays[i].direction);

                    // Check if face should be culled
                    if ((culling == CullingMode::BACK_FACE && dot > 0.0f) ||
                        (culling == CullingMode::FRONT_FACE && dot < 0.0f)) {
                        culled = true;
                    }
                }
            }

            // Only fill result if not culled
            if (!culled) {
                hr.didHit = true;
                hr.distance = rayhit.ray.tfar[i];
                hr.geomID = rayhit.hit.geomID[i];
                hr.primID = rayhit.hit.primID[i];
                hr.u = rayhit.hit.u[i];
                hr.v = rayhit.hit.v[i];
                hr.origin = rays[i].origin + rays[i].direction * hr.distance;
                hr.normal = computeNormal(hr, rays[i]);
            }
        }
    }
}

bool RayTracer::isOccluded(const Ray& ray) const {
    RTCRay rtcRay{
        .org_x = ray.origin.x,
        .org_y = ray.origin.y,
        .org_z = ray.origin.z,
        .tnear = ray.tNear,
        .dir_x = ray.direction.x,
        .dir_y = ray.direction.y,
        .dir_z = ray.direction.z,
        .time = 0.0f,
        .tfar = ray.tFar,
        .mask = ray.mask,
        .id = 0,
        .flags = 0
    };

    RTCOccludedArguments args{};
    rtcInitOccludedArguments(&args);
    args.context = &m_rqctx;

    rtcOccluded1(m_scene, &rtcRay, &args);

    return rtcRay.tfar < 0.0f;
}


glm::vec3 RayTracer::computeNormal(const HitResult& hit, const Ray& ray) const {
    const RTCGeometryType geomType = m_sceneManager->getGeometryType(hit.geomID);
    RTCGeometry geom = rtcGetGeometry(m_scene, hit.geomID);

    if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
        // Unified triangle/quad handling
        const auto* vertices = static_cast<const float*>(rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
        const auto* indices = static_cast<const unsigned*>(rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

        const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
        const unsigned i0 = indices[hit.primID * stride + 0];
        const unsigned i1 = indices[hit.primID * stride + 1];
        const unsigned i2 = indices[hit.primID * stride + 2];

        // Direct pointer arithmetic for vertex access
        const float* v0 = vertices + i0 * 3;
        const float* v1 = vertices + i1 * 3;
        const float* v2 = vertices + i2 * 3;

        // Compute edges and cross product inline
        const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
        const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
        glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));

        // Flip normal if facing away
        if (glm::dot(normal, ray.direction) > 0.0f) {
            normal = -normal;
        }

        return normal;
    }

    if (geomType == RTC_GEOMETRY_TYPE_SPHERE_POINT) {
        struct Sphere4 {
            float x, y, z, r;
        };

        const auto* spheres = static_cast<const Sphere4*>(rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
        const glm::vec3 center(spheres[hit.primID].x, spheres[hit.primID].y, spheres[hit.primID].z);

        return glm::normalize(hit.origin - center);
    }

    // Fallback
    return -glm::normalize(ray.direction);
}