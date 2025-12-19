#include "ray_tracer.h"
#include <embree4/rtcore.h>

RayTracer::RayTracer(const SceneManager* sceneManager)
    : m_scene(sceneManager->getScene())
    , m_sceneManager(sceneManager)
{
    rtcInitRayQueryContext(&m_rqctx);
}

HitResult RayTracer::intersect(const Ray &ray) const {
    HitResult result;

    RTCRayHit rayHit{
        .ray = {
            .org_x = ray.origin.x,
            .org_y = ray.origin.y,
            .org_z = ray.origin.z,
            .tnear = 0.0001f,
            .dir_x = ray.direction.x,
            .dir_y = ray.direction.y,
            .dir_z = ray.direction.z,
            .time = 0.0f,
            .tfar = ray.tFar,
            .mask = ray.mask,
            .id = 0,
            .flags = RTC_RAY_QUERY_FLAG_COHERENT
        },
        .hit = {
            .geomID = RTC_INVALID_GEOMETRY_ID
        }
    };
    RTCIntersectArguments args{};
    rtcInitIntersectArguments(&args);
    args.context = const_cast<RTCRayQueryContext*>(&m_rqctx);
    rtcIntersect1(m_scene, &rayHit, &args);

    if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
    {
        result.didHit = true;
        result.distance = rayHit.ray.tfar;
        result.geomID = rayHit.hit.geomID;
        result.primID = rayHit.hit.primID;
        result.u = rayHit.hit.u;
        result.v = rayHit.hit.v;
        result.origin = ray.origin + ray.direction * result.distance;

        glm::vec3 geometricNormal(
            rayHit.hit.Ng_x,
            rayHit.hit.Ng_y,
            rayHit.hit.Ng_z
        );

        result.normal = glm::normalize(geometricNormal);

        // ensure facing the ray
        if (glm::dot(result.normal, ray.direction) > 0.0f)
            result.normal = -result.normal;
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
        RTCRay shadow{
                .org_x = ray.origin.x,
                .org_y = ray.origin.y,
                .org_z = ray.origin.z,
                .tnear = 0.001f,
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

        rtcOccluded1(m_scene, &shadow, &args);

        return shadow.tfar < 0.0f;
}