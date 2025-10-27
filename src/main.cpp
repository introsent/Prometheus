//
// Created by minaj on 4/4/2025.
//
#include <iostream>
#include <memory>

#include "geometries/triangle.h"
#include "init/embree_device.h"
#include "init/embree_scene.h"
#include "ray/ray_tracer.h"
#include "structs/vertex.h"


int main()
{
    const auto devicePtr = std::make_unique<EmbreeDevice>();
    const auto scenePtr = std::make_unique<EmbreeScene>(devicePtr.get());

    std::vector <Vertex> vertices = {Vertex{ 0.f, 0.f, 0.f}, Vertex{1.f, 0.f, 0.f}, Vertex{0.f, 1.f, 0.f}};
    const auto trianglePtr = std::make_unique<Triangle>(vertices, devicePtr.get());

    trianglePtr->commit();
    trianglePtr->attach(scenePtr.get());
    trianglePtr->release();

    scenePtr->commit();

    const auto rayTracerPtr = std::make_unique<RayTracer>(scenePtr->handle());
    Ray ray{glm::vec3(0, 0, -1), glm::vec3(0, 0, 1)};

    if (const HitResult hit = rayTracerPtr->intersect(ray); hit.didHit) {
        std::cout << "Intersection at t = " << hit.distance
                  << ", geometry ID: " << hit.geomID << std::endl;
    }

    scenePtr->release();
    devicePtr->release();
}
