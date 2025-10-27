//
// Created by minaj on 4/4/2025.
//
#include <iostream>
#include <memory>

// Project includes
#include "vertex.h"
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

void createReferenceScene(SceneManager* pScene)
{
    // Setup Cook-Torrence materials
    const unsigned char matCT_GrayRoughMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 1.f));
    const unsigned char matCT_GrayMediumMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 0.6f));
    const unsigned char matCT_GraySmoothMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 0.1f));

    const unsigned char matCT_GrayRoughPlastic = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.75f, 0.75f, 0.75f), 0.f, 1.f));
    const unsigned char matCT_GrayMediumPlastic = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.75f, 0.75f, 0.75f), 0.f, 0.6f));
    const unsigned char matCT_GraySmoothPlastic = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.75f, 0.75f, 0.75f), 0.f, 0.1f));

    // Lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));
    const unsigned char matLambert_White = pScene->addMaterial(
        new Material_Lambert(colors::white, 1.f));

    // Room planes
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    // Bottom row spheres (metals with varying roughness)
    pScene->addSphere({-1.75f, 1.f, 0.f}, 0.75f, matCT_GrayRoughMetal);
    pScene->addSphere({0.f, 1.f, 0.f}, 0.75f, matCT_GrayMediumMetal);
    pScene->addSphere({1.75f, 1.f, 0.f}, 0.75f, matCT_GraySmoothMetal);

    // Top row spheres (plastics with varying roughness)
    pScene->addSphere({-1.75f, 3.f, 0.f}, 0.75f, matCT_GrayRoughPlastic);
    pScene->addSphere({0.f, 3.f, 0.f}, 0.75f, matCT_GrayMediumPlastic);
    pScene->addSphere({1.75f, 3.f, 0.f}, 0.75f, matCT_GraySmoothPlastic);

    // Add triangles for culling demonstration
    // Triangle 1: Back-face culling
    std::vector<Vertex> triangle1 = {
        {-2.5f, 6.f, 0.f},   // top
        {-1.f, 4.5f, 0.f},   // bottom-right
        {-2.5f, 4.5f, 0.f}   // bottom-left
    };
    pScene->addTriangle(triangle1, matLambert_White, CullingMode::BACK_FACE);

    // Triangle 2: Front-face culling
    std::vector<Vertex> triangle2 = {
        {0.f, 6.f, 0.f},     // top
        {-0.75f, 4.5f, 0.f}, // bottom-left
        {0.75f, 4.5f, 0.f}   // bottom-right
    };
    pScene->addTriangle(triangle2, matLambert_White, CullingMode::FRONT_FACE);

    // Triangle 3: No culling
    std::vector<Vertex> triangle3 = {
        {2.5f, 6.f, 0.f},    // top
        {1.75f, 4.5f, 0.f},  // bottom-right
        {2.5f, 4.5f, 0.f}    // bottom-left
    };
    pScene->addTriangle(triangle3, matLambert_White, CullingMode::NONE);

    // Lights
    pScene->addLight(new Light(
        glm::vec3(0.f, 5.f, 5.f),
        glm::vec3(0.f, 0.f, 0.f),
        glm::vec3(1.f, 0.61f, 0.45f),  // Warm color
        50.f,
        LightType::Point
    ));

    pScene->addLight(new Light(
        glm::vec3(-2.5f, 5.f, -5.f),
        glm::vec3(0.f, 0.f, 0.f),
        glm::vec3(1.f, 0.8f, 0.45f),   // Warm-neutral color
        70.f,
        LightType::Point
    ));

    pScene->addLight(new Light(
        glm::vec3(2.5f, 2.5f, -5.f),
        glm::vec3(0.f, 0.f, 0.f),
        glm::vec3(0.34f, 0.47f, 0.68f), // Cool color
        50.f,
        LightType::Point
    ));
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
    createReferenceScene(pScene.get());
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