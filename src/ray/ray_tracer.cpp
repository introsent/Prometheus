//
// Created by minaj on 10/27/2025.
// Optimized for maximum performance with Embree 4
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

// Non-static member function for the filter
void RayTracer::intersectionFilter(const struct RTCFilterFunctionNArguments* args) const {
    int* valid = args->valid;
    RTCHitN* hit = args->hit;
    RTCRayN* ray = args->ray;

    const unsigned int N = args->N;

    // Handle each lane in the packet
    for (unsigned int i = 0; i < N; i++) {
        if (!valid[i]) continue;

        // Extract geometry ID and primitive ID for this lane
        const unsigned int geomID = RTCHitN_geomID(hit, N, i);
        const unsigned int primID = RTCHitN_primID(hit, N, i);

        // Check culling for this geometry
        CullingMode culling = m_sceneManager->getCullingMode(geomID);
        if (culling == CullingMode::NONE) {
            continue; // Keep this hit
        }

        RTCGeometry geom = rtcGetGeometry(m_scene, geomID);
        const RTCGeometryType geomType = m_sceneManager->getGeometryType(geomID);

        if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
            const auto* vertices = static_cast<const float*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
            const auto* indices = static_cast<const unsigned*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

            const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
            const unsigned i0 = indices[primID * stride + 0];
            const unsigned i1 = indices[primID * stride + 1];
            const unsigned i2 = indices[primID * stride + 2];

            const float* v0 = vertices + i0 * 3;
            const float* v1 = vertices + i1 * 3;
            const float* v2 = vertices + i2 * 3;

            // Compute face normal
            const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
            const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
            const glm::vec3 normal = glm::cross(edge1, edge2);

            // Get ray direction for this lane
            const glm::vec3 rayDir(
                RTCRayN_dir_x(ray, N, i),
                RTCRayN_dir_y(ray, N, i),
                RTCRayN_dir_z(ray, N, i)
            );

            const float dot = glm::dot(normal, rayDir);

            // Cull based on face orientation
            if ((culling == CullingMode::BACK_FACE && dot > 0.0f) ||
                (culling == CullingMode::FRONT_FACE && dot < 0.0f)) {
                valid[i] = 0; // Reject this hit - ray will continue through
            }
        }
    }
}
HitResult RayTracer::intersect(const Ray &ray) const {
    HitResult result;
    float currentT = ray.tNear;
    constexpr int maxHits = 8; // Prevent infinite loops
    int hitCount = 0;

    while (currentT < ray.tFar && hitCount < maxHits) {
        RTCRayHit rayHit{
            .ray = {
                .org_x = ray.origin.x,
                .org_y = ray.origin.y,
                .org_z = ray.origin.z,
                .tnear = currentT,
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
            return result; // No more hits
        }

        // Check if this hit should be culled
        bool culled = false;
        CullingMode culling = m_sceneManager->getCullingMode(rayHit.hit.geomID);
        if (culling != CullingMode::NONE) {
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

                const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
                const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
                const glm::vec3 normal = glm::cross(edge1, edge2);
                const float dot = glm::dot(normal, ray.direction);

                if ((culling == CullingMode::BACK_FACE && dot > 0.0f) ||
                    (culling == CullingMode::FRONT_FACE && dot < 0.0f)) {
                    culled = true;
                }
            }
        }

        if (!culled) {
            // Valid non-culled hit
            result.didHit = true;
            result.distance = rayHit.ray.tfar;
            result.geomID = rayHit.hit.geomID;
            result.primID = rayHit.hit.primID;
            result.u = rayHit.hit.u;
            result.v = rayHit.hit.v;
            result.origin = ray.origin + ray.direction * result.distance;
            result.normal = computeNormal(result, ray);
            return result;
        } else {
            // Continue from just beyond this culled hit
            currentT = rayHit.ray.tfar + 0.0001f;
            hitCount++;
        }
    }

    return result; // No non-culled hits found
}

void RayTracer::intersectPacket16(const std::vector<Ray> &rays, std::vector<HitResult> &results) const {
    const size_t N = rays.size();
    assert(N <= 16);

    results.resize(N);

    // For packet rays, we'll do a simpler approach - single intersection per ray
    // A proper continuation approach for packets is complex
    for (size_t i = 0; i < N; ++i) {
        results[i] = intersect(rays[i]);
    }
}

bool RayTracer::isOccluded(const Ray& ray) const {
    // For occlusion, we need to check if ANY non-culled geometry occludes
    float currentT = ray.tNear;
    const int maxHits = 8;
    int hitCount = 0;

    while (currentT < ray.tFar && hitCount < maxHits) {
        RTCRay rtcRay{
            .org_x = ray.origin.x,
            .org_y = ray.origin.y,
            .org_z = ray.origin.z,
            .tnear = currentT,
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
        args.context = const_cast<RTCRayQueryContext*>(&m_rqctx);

        rtcOccluded1(m_scene, &rtcRay, &args);

        if (rtcRay.tfar >= 0.0f) {
            return false; // No occlusion
        }

        // We got an occlusion, but we need to check if it was culled
        // Since rtcOccluded1 doesn't give us hit info, we need to do a separate intersect
        // to check if the occluding geometry should be culled

        RTCRayHit rayHit{
            .ray = {
                .org_x = ray.origin.x,
                .org_y = ray.origin.y,
                .org_z = ray.origin.z,
                .tnear = currentT,
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

        RTCIntersectArguments intersectArgs{};
        rtcInitIntersectArguments(&intersectArgs);
        args.context = const_cast<RTCRayQueryContext*>(&m_rqctx);

        rtcIntersect1(m_scene, &rayHit, &intersectArgs);

        if (rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
            return false; // Shouldn't happen if we got occlusion
        }

        // Check if this occluding hit should be culled
        bool culled = false;
        CullingMode culling = m_sceneManager->getCullingMode(rayHit.hit.geomID);
        if (culling != CullingMode::NONE) {
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

                const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
                const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
                const glm::vec3 normal = glm::cross(edge1, edge2);
                const float dot = glm::dot(normal, ray.direction);

                if ((culling == CullingMode::BACK_FACE && dot > 0.0f) ||
                    (culling == CullingMode::FRONT_FACE && dot < 0.0f)) {
                    culled = true;
                }
            }
        }

        if (!culled) {
            return true; // Valid non-culled occlusion
        } else {
            // Continue from just beyond this culled hit
            currentT = rayHit.ray.tfar + 0.0001f;
            hitCount++;
        }
    }

    return false; // No non-culled occlusions found
}

// computeNormal remains the same as your original
glm::vec3 RayTracer::computeNormal(const HitResult& hit, const Ray& ray) const {
    const RTCGeometryType geomType = m_sceneManager->getGeometryType(hit.geomID);
    RTCGeometry geom = rtcGetGeometry(m_scene, hit.geomID);

    if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
        const auto* vertices = static_cast<const float*>(rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
        const auto* indices = static_cast<const unsigned*>(rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

        const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
        const unsigned i0 = indices[hit.primID * stride + 0];
        const unsigned i1 = indices[hit.primID * stride + 1];
        const unsigned i2 = indices[hit.primID * stride + 2];

        const float* v0 = vertices + i0 * 3;
        const float* v1 = vertices + i1 * 3;
        const float* v2 = vertices + i2 * 3;

        const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
        const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
        glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));

        // Always face the normal toward the camera
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