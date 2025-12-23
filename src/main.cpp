//
// Created by minaj on 4/4/2025.
//
#include <iostream>
#include <memory>

// Project includes
#include "mesh.h"
#include "vertex.h"
#include "render/timer.h"
#include "render/scene_manager.h"
#include "camera/camera.h"
#include "parser/obj_parser.h"
#include "render/renderer.h"
#include "structs/material.h"
#include "structs/global_indices.h"


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

    // Add triangle
    std::vector<Vertex> baseTriangle = {
        {{-0.75f, 1.5f, 0.f}},   // top
        {{-0.75f, 0.f, 0.f}},    // bottom-left (swapped)
        {{0.75f, 0.f, 0.f}},     // bottom-right (swapped)
    };

    // Triangle
    const unsigned int triId = pScene->addTriangle(baseTriangle, matLambert_White);
    if (Triangle* tri = pScene->getTriangle(triId)) {
        tri->translate({0.f, 4.5f, 0.f});

        tri->setUpdateFunc([tri](float, float totalTime)
        {
            const float yaw =
                (std::cos(totalTime) + 1.f) * glm::pi<float>();

            tri->rotateY(yaw);
            tri->updateAABB();
            tri->updateTransforms();
        });
    }


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

void createBunnyScene(SceneManager* pScene)
{
    // Lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));
    const unsigned char matLambert_White = pScene->addMaterial(
        new Material_Lambert(colors::white, 1.f));

    // Room planes (same as reference scene)
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    // Load bunny mesh
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {
        const unsigned int bunnyId = pScene->addMesh(bunnyVertices, bunnyIndices, matLambert_White);

        if (Mesh* pBunny = pScene->getMesh(bunnyId)) {
            // Apply transformations
            pBunny->rotateY(glm::radians(180.0f)); // PI radians = 180 degrees
            pBunny->scale(glm::vec3(2.f, 2.f, 2.f));
            pBunny->updateAABB();
            pBunny->updateTransforms();

            //pBunny->setUpdateFunc([pBunny](float, float totalTime)
            //{
            //    const float yaw =
            //        (std::cos(totalTime) + 1.f) * glm::pi<float>();
            //
            //    pBunny->rotateY(yaw);
            //    pBunny->updateAABB();
            //    pBunny->updateTransforms();
            //});
        }
    } else {
        std::cerr << "Failed to load lowpoly_bunny.obj" << std::endl;
    }

    // Lights (same as reference scene)
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

void createSceneA(SceneManager* pScene)
{
    // Lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    // Room planes (same as reference scene)
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    std::vector<Vertex> ceilingLight_tri1 = {
        {{-2.0f, 9.99f, 3.0f}, {0.0f, -1.0f, 0.0f}},  // Front-left
        {{ 2.0f, 9.99f, 3.0f}, {0.0f, -1.0f, 0.0f}},  // Front-right
        {{ 2.0f, 9.99f, 7.0f}, {0.0f, -1.0f, 0.0f}}   // Back-right
    };

    std::vector<Vertex> ceilingLight_tri2 = {
        {{-2.0f, 9.99f, 3.0f}, {0.0f, -1.0f, 0.0f}},  // Front-left
        {{ 2.0f, 9.99f, 7.0f}, {0.0f, -1.0f, 0.0f}},  // Back-right
        {{-2.0f, 9.99f, 7.0f}, {0.0f, -1.0f, 0.0f}}   // Back-left
    };

    glm::vec3 emission(1.0f, 0.0f, 0.0f);
    float intensity = 20.0f;  // higher intensity for smaller area

    pScene->addTriangleAreaLight(ceilingLight_tri1, emission, intensity);
    pScene->addTriangleAreaLight(ceilingLight_tri2, emission, intensity);
}

void createSceneB(SceneManager* pScene)
{
     // Lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    // Room planes (same as reference scene)
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    // Load bunny mesh AS AN AREA LIGHT
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {

        // apply transformations to vertices BEFORE adding as area light
        glm::mat4 transform = glm::mat4(1.0f);
        transform = glm::rotate(transform, glm::radians(180.0f), glm::vec3(0.f, 1.f, 0.f));
        transform = glm::scale(transform, glm::vec3(2.f, 2.f, 2.f));

        // transform all vertices
        for (auto& vertex : bunnyVertices) {
            glm::vec4 pos(vertex.position.x, vertex.position.y, vertex.position.z, 1.0f);
            glm::vec4 transformedPos = transform * pos;
            vertex.position.x = transformedPos.x;
            vertex.position.y = transformedPos.y;
            vertex.position.z = transformedPos.z;

            // transform normals (use inverse transpose for normals)
            glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(transform)));
            glm::vec3 normal(vertex.normal.x, vertex.normal.y, vertex.normal.z);
            glm::vec3 transformedNormal = glm::normalize(normalMatrix * normal);
            vertex.normal.x = transformedNormal.x;
            vertex.normal.y = transformedNormal.y;
            vertex.normal.z = transformedNormal.z;
        }

        // define emission properties
        glm::vec3 emission(1.0f, 1.0f, 1.0f);
        float intensity = 10.0f;  // high intensity since bunny is the only light source

        // add bunny as a MESH AREA LIGHT
        pScene->addMeshAreaLight(
            bunnyVertices,
            bunnyIndices,
            emission,
            intensity
        );

        std::cout << "Emissive bunny added as mesh area light with "
                  << (bunnyIndices.size() / 3) << " triangle lights" << std::endl;

    } else {
        std::cerr << "Failed to load lowpoly_bunny.obj" << std::endl;
    }
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
    createSceneA(pScene.get());
    pScene->commit();


    // Start loop
    pTimer->reset();
    pTimer->start();

    // Uncomment to start benchmark
    //pTimer->startBenchmark(10); // 1000 frames benchmark

    float printTimer = 0.f;

    while (!pRenderer->shouldQuit())
    {
        //--------- Update ---------
        pScene->update(pTimer.get());

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