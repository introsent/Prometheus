//
// Created by minaj on 4/4/2025.
//
#include <iostream>
#include <memory>

// Project includes
#include "render/timer.h"
#include "render/scene_manager.h"
#include "camera/camera.h"
#include "render/renderer.h"
#include "structs/material.h"

void createBaseColorScene(SceneManager* pScene)
{
    // Setup materials
    constexpr unsigned char matId_Solid_Red = 0;  // Default red material
    const unsigned char matId_Solid_Blue = pScene->addMaterial(new Material_SolidColor{colors::blue});
    const unsigned char matId_Solid_Yellow = pScene->addMaterial(new Material_SolidColor{colors::yellow});
    const unsigned char matId_Solid_Green = pScene->addMaterial(new Material_SolidColor{colors::green});
    const unsigned char matId_Solid_Magenta = pScene->addMaterial(new Material_SolidColor{colors::magenta});

    // Add spheres
    pScene->addSphere({-25.f, 0.f, 100.f}, 50.f, matId_Solid_Red);
    pScene->addSphere({25.f, 0.f, 100.f}, 50.f, matId_Solid_Blue);

    // Add planes
    pScene->addPlane({-75.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matId_Solid_Green);
    pScene->addPlane({75.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matId_Solid_Green);
    pScene->addPlane({0.f, -75.f, 0.f}, {0.f, 1.f, 0.f}, matId_Solid_Yellow);
    pScene->addPlane({0.f, 75.f, 0.f}, {0.f, -1.f, 0.f}, matId_Solid_Yellow);
    pScene->addPlane({0.f, 0.f, 125.f}, {0.f, 0.f, -1.f}, matId_Solid_Magenta);
}

void createBaseColorWithLightScene(SceneManager* pScene)
{
    // Setup materials
    constexpr unsigned char matId_Solid_Red = 0;  // Default red material
    const unsigned char matId_Solid_Blue = pScene->addMaterial(new Material_SolidColor{colors::blue});
    const unsigned char matId_Solid_Yellow = pScene->addMaterial(new Material_SolidColor{colors::yellow});
    const unsigned char matId_Solid_Green = pScene->addMaterial(new Material_SolidColor{colors::green});
    const unsigned char matId_Solid_Magenta = pScene->addMaterial(new Material_SolidColor{colors::magenta});

    //Plane
    pScene->addPlane({ -5.f, 0.f, 0.f }, { 1.f, 0.f, 0.f }, matId_Solid_Green );
    pScene->addPlane({ 5.f, 0.f, 0.f }, { -1.f, 0.f, 0.f }, matId_Solid_Green);
    pScene->addPlane({ 0.f, 0.f, 0.f }, { 0.f,  1.f, 0.f }, matId_Solid_Yellow);
    pScene->addPlane({ 0.f, 10.f, 0.f }, { 0.f, -1.f, 0.f}, matId_Solid_Yellow);
    pScene->addPlane({ 0.f,  0.f, 10.f }, { 0.f,  0.f, -1.f }, matId_Solid_Magenta);

    //Spheres
    pScene->addSphere({ -1.75f, 1.f, 0.f }, 0.75f, matId_Solid_Red);
    pScene->addSphere({ 0.f, 1.f, 0.f }, 0.75f, matId_Solid_Blue);
    pScene->addSphere({ 1.75f, 1.f, 0.f }, 0.75f, matId_Solid_Red);
    pScene->addSphere({ -1.75f, 3.f, 0.f }, 0.75f, matId_Solid_Blue);
    pScene->addSphere({ 0.f, 3.f, 0.f }, 0.75f, matId_Solid_Red);
    pScene->addSphere({ 1.75f, 3.f, 0.f }, 0.75f, matId_Solid_Blue);

    //Lights
    pScene->addLight(new Light({ 0.f, 5.f, -5.f }, glm::vec3{0.f, 0.f, 0.f},
        colors::white, 70.f, LightType::Point));
}

int main(int argc, char* args[])
{
    // Unreferenced parameters
    (void)argc;
    (void)args;

    constexpr uint32_t WIDTH = 640;
    constexpr uint32_t HEIGHT = 480;

    // Initialize framework
    auto pTimer = std::make_unique<Timer>();
    auto pRenderer = std::make_unique<Renderer>(WIDTH, HEIGHT);

    if (!pRenderer->initialize()) {
        std::cerr << "Failed to initialize renderer" << std::endl;
        return -1;
    }

    // Create scene
    auto pScene = std::make_unique<SceneManager>();
    auto pCamera = std::make_unique<Camera>(glm::vec3{0.f, 3.f, -9.f},
                                            45.f,
                                            static_cast<float>(WIDTH) / static_cast<float>(HEIGHT));
    createBaseColorWithLightScene(pScene.get());
    pScene->commit();


    // Start loop
    pTimer->reset();
    pTimer->start();

    // Uncomment to start benchmark
    // pTimer->startBenchmark(1000); // 1000 frames benchmark

    float printTimer = 0.f;

    while (!pRenderer->shouldQuit())
    {
        //--------- Update ---------
        // If your scene has update functionality, uncomment:
        // pScene->Update(pTimer.get());

        //--------- Render ---------
        pRenderer->render(*pCamera, *pScene);
        pRenderer->present();

        //--------- Timer ---------
        pTimer->update();
        printTimer += pTimer->getElapsed();

        if (printTimer >= 1.f)
        {
            printTimer = 0.f;
            std::cout << "dFPS: " << pTimer->getdFPS() << std::endl;
        }
    }

    pTimer->stop();

    std::cout << "Application shutdown complete" << std::endl;
    return 0;
}