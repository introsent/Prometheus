//
// Created by minaj on 4/4/2025.
//

#include <embree4/rtcore.h>
#include <limits>
#include <iostream>
#include <memory>

#include "geometries/triangle.h"
#include "init/embree_device.h"
#include "init/embree_scene.h"
#include "structs/vertex.h"


int main()
{
    auto devicePtr = std::make_unique<EmbreeDevice>();
    auto scenePtr = std::make_unique<EmbreeScene>(devicePtr.get());

    std::vector <Vertex> vertices = {Vertex{ 0.f, 0.f, 0.f}, Vertex{1.f, 0.f, 0.f}, Vertex{0.f, 1.f, 0.f}};
    auto trianglePtr = std::make_unique<Triangle>(vertices, devicePtr.get());

    trianglePtr->commit();
    trianglePtr->attach(scenePtr.get());
    trianglePtr->release();

    scenePtr->commit();


    RTCRayHit rayhit = {};
    rayhit.ray.org_x  = 0.f; rayhit.ray.org_y = 0.f; rayhit.ray.org_z = -1.f;
    rayhit.ray.dir_x  = 0.f; rayhit.ray.dir_y = 0.f; rayhit.ray.dir_z =  1.f;
    rayhit.ray.tnear  = 0.f;
    rayhit.ray.tfar   = std::numeric_limits<float>::infinity();
    rayhit.ray.flags  = 0;
    rayhit.ray.time   = 0.0f;
    rayhit.ray.mask   = 0xFFFFFFFF;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

    RTCRayQueryContext rqctx;
    rtcInitRayQueryContext(&rqctx);

    RTCIntersectArguments args;
    rtcInitIntersectArguments(&args);
    args.context = &rqctx;

    rtcIntersect1(scenePtr->getScene(), &rayhit, &args);

    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
        std::cout << "Intersection at t = " << rayhit.ray.tfar << std::endl;
    } else {
        std::cout << "No Intersection" << std::endl;
    }

    scenePtr->release();
    devicePtr->release();
}
