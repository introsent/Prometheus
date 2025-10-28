//
// Fixed: Proper geometry ID handling for Embree 4
//

#include "ray_tracer.h"
#include <embree4/rtcore.h>

RayTracer::RayTracer(const SceneManager* sceneManager)
    : m_scene(sceneManager->getScene())
    , m_sceneManager(sceneManager)
{
    rtcInitRayQueryContext(&m_rqctx);
}

// Inline helper for fast culling check
inline bool shouldCullHit(
    const glm::vec3& rayDir,
    const float* v0, const float* v1, const float* v2,
    CullingMode culling
) {
    if (culling == CullingMode::NONE) return false;

    const glm::vec3 edge1(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
    const glm::vec3 edge2(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
    const glm::vec3 normal = glm::cross(edge1, edge2);
    const float dot = glm::dot(normal, rayDir);

    return (culling == CullingMode::BACK_FACE && dot > 0.0f) ||
           (culling == CullingMode::FRONT_FACE && dot < 0.0f);
}

HitResult RayTracer::intersect(const Ray &ray) const {
    HitResult result;
    float currentT = ray.tNear;
    constexpr int maxIterations = 4;

    for (int iteration = 0; iteration < maxIterations && currentT < ray.tFar; ++iteration) {
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
            return result;
        }

        const unsigned geomID = rayHit.hit.geomID;
        const CullingMode culling = m_sceneManager->getCullingMode(geomID);

        // Fast path: no culling
        if (culling == CullingMode::NONE) {
            result.didHit = true;
            result.distance = rayHit.ray.tfar;
            result.geomID = geomID;
            result.primID = rayHit.hit.primID;
            result.u = rayHit.hit.u;
            result.v = rayHit.hit.v;
            result.origin = ray.origin + ray.direction * result.distance;
            result.normal = computeNormal(result, ray);
            return result;
        }

        // Check culling - get geometry directly from Embree
        const RTCGeometryType geomType = m_sceneManager->getGeometryType(geomID);
        bool culled = false;

        if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
            RTCGeometry geom = rtcGetGeometry(m_scene, geomID);
            const float* vertices = static_cast<const float*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
            const unsigned* indices = static_cast<const unsigned*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

            if (vertices && indices) {
                const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
                const unsigned i0 = indices[rayHit.hit.primID * stride + 0];
                const unsigned i1 = indices[rayHit.hit.primID * stride + 1];
                const unsigned i2 = indices[rayHit.hit.primID * stride + 2];

                culled = shouldCullHit(
                    ray.direction,
                    vertices + i0 * 3,
                    vertices + i1 * 3,
                    vertices + i2 * 3,
                    culling
                );
            }
        }

        if (!culled) {
            result.didHit = true;
            result.distance = rayHit.ray.tfar;
            result.geomID = geomID;
            result.primID = rayHit.hit.primID;
            result.u = rayHit.hit.u;
            result.v = rayHit.hit.v;
            result.origin = ray.origin + ray.direction * result.distance;
            result.normal = computeNormal(result, ray);
            return result;
        }

        currentT = rayHit.ray.tfar + 0.0001f;
    }

    return result;
}

void RayTracer::intersectPacket16(const std::vector<Ray> &rays, std::vector<HitResult> &results) const {
    const size_t N = rays.size();
    results.resize(N);

    // Process rays in scalar mode (proper packet requires aligned structures)
    for (size_t i = 0; i < N; ++i) {
        results[i] = intersect(rays[i]);
    }
}

bool RayTracer::isOccluded(const Ray& ray) const {
    float currentT = ray.tNear;
    constexpr int maxIterations = 4;

    for (int iteration = 0; iteration < maxIterations && currentT < ray.tFar; ++iteration) {
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
            return false;
        }

        const unsigned geomID = rayHit.hit.geomID;
        const CullingMode culling = m_sceneManager->getCullingMode(geomID);

        // Fast path
        if (culling == CullingMode::NONE) {
            return true;
        }

        // Check culling
        const RTCGeometryType geomType = m_sceneManager->getGeometryType(geomID);
        bool culled = false;

        if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
            RTCGeometry geom = rtcGetGeometry(m_scene, geomID);
            const float* vertices = static_cast<const float*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
            const unsigned* indices = static_cast<const unsigned*>(
                rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

            if (vertices && indices) {
                const unsigned stride = (geomType == RTC_GEOMETRY_TYPE_QUAD) ? 4 : 3;
                const unsigned i0 = indices[rayHit.hit.primID * stride + 0];
                const unsigned i1 = indices[rayHit.hit.primID * stride + 1];
                const unsigned i2 = indices[rayHit.hit.primID * stride + 2];

                culled = shouldCullHit(
                    ray.direction,
                    vertices + i0 * 3,
                    vertices + i1 * 3,
                    vertices + i2 * 3,
                    culling
                );
            }
        }

        if (!culled) {
            return true;
        }

        currentT = rayHit.ray.tfar + 0.0001f;
    }

    return false;
}

glm::vec3 RayTracer::computeNormal(const HitResult& hit, const Ray& ray) const {
    const RTCGeometryType geomType = m_sceneManager->getGeometryType(hit.geomID);
    RTCGeometry geom = rtcGetGeometry(m_scene, hit.geomID);

    if (geomType == RTC_GEOMETRY_TYPE_TRIANGLE || geomType == RTC_GEOMETRY_TYPE_QUAD) {
        const float* vertices = static_cast<const float*>(
            rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
        const unsigned* indices = static_cast<const unsigned*>(
            rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_INDEX, 0));

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

        if (glm::dot(normal, ray.direction) > 0.0f) {
            normal = -normal;
        }

        return normal;
    }

    if (geomType == RTC_GEOMETRY_TYPE_SPHERE_POINT) {
        struct Sphere4 { float x, y, z, r; };
        const auto* spheres = static_cast<const Sphere4*>(
            rtcGetGeometryBufferData(geom, RTC_BUFFER_TYPE_VERTEX, 0));
        const glm::vec3 center(spheres[hit.primID].x, spheres[hit.primID].y, spheres[hit.primID].z);
        return glm::normalize(hit.origin - center);
    }

    return -glm::normalize(ray.direction);
}