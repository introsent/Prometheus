//
// Created by minaj on 4/4/2025.
//
#include <iostream>
#include <memory>

// Project includes
#include <fstream>

#include "mesh.h"
#include "vertex.h"
#include "render/timer.h"
#include "render/scene_manager.h"
#include "camera/camera.h"
#include "lights/mesh_area_light.h"
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
    // setup Cook-Torrence materials
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

    // lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));
    const unsigned char matLambert_White = pScene->addMaterial(
        new Material_Lambert(colors::white, 1.f));

    // room planes
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    // bottom row spheres (metals with varying roughness)
    pScene->addSphere({-1.75f, 1.f, 0.f}, 0.75f, matCT_GrayRoughMetal);
    pScene->addSphere({0.f, 1.f, 0.f}, 0.75f, matCT_GrayMediumMetal);
    pScene->addSphere({1.75f, 1.f, 0.f}, 0.75f, matCT_GraySmoothMetal);

    // top row spheres (plastics with varying roughness)
    pScene->addSphere({-1.75f, 3.f, 0.f}, 0.75f, matCT_GrayRoughPlastic);
    pScene->addSphere({0.f, 3.f, 0.f}, 0.75f, matCT_GrayMediumPlastic);
    pScene->addSphere({1.75f, 3.f, 0.f}, 0.75f, matCT_GraySmoothPlastic);

    // add triangle
    std::vector<Vertex> baseTriangle = {
        {{-0.75f, 1.5f, 0.f}},   // top
        {{-0.75f, 0.f, 0.f}},    // bottom-left (swapped)
        {{0.75f, 0.f, 0.f}},     // bottom-right (swapped)
    };

    // triangle
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
    // lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));
    const unsigned char matLambert_White = pScene->addMaterial(
        new Material_Lambert(colors::white, 1.f));

    // room planes (same as reference scene)
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    // load bunny mesh
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {
        const unsigned int bunnyId = pScene->addMesh(bunnyVertices, bunnyIndices, matLambert_White);

        if (Mesh* pBunny = pScene->getMesh(bunnyId)) {
            // apply transformations
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

    // lights (same as reference scene)
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
    // lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    // room planes (same as reference scene)
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

    // room planes (same as reference scene)
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT

    // load bunny mesh AS AN AREA LIGHT
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {

        // apply transformations to vertices BEFORE adding as area light
        auto transform = glm::mat4(1.0f);
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
    {0.9f, 0.9f, 0.9f},
            emission,
            intensity
        );

        std::cout << "Emissive bunny added as mesh area light with "
                  << (bunnyIndices.size() / 3) << " triangle lights" << std::endl;

    } else {
        std::cerr << "Failed to load lowpoly_bunny.obj" << std::endl;
    }
}

void createSceneC(SceneManager* pScene)
{
    // lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    // setup Cook-Torrence materials
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

    // room planes (same as reference scene)
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT



    const unsigned char matLambert_White = pScene->addMaterial(
        new Material_Lambert(colors::white, 1.f));

    // bottom row spheres (metals with varying roughness)
    pScene->addSphere({-1.75f, 1.f, 5.f}, 0.75f, matCT_GrayRoughMetal);
    pScene->addSphere({0.f, 1.f, 5.f}, 0.75f, matCT_GrayMediumMetal);
    pScene->addSphere({1.75f, 1.f, 5.f}, 0.75f, matCT_GraySmoothMetal);

    // top row spheres (plastics with varying roughness)
    pScene->addSphere({-1.75f, 3.f, 5.f}, 0.75f, matCT_GrayRoughPlastic);
    pScene->addSphere({0.f, 3.f, 5.f}, 0.75f, matCT_GrayMediumPlastic);
    pScene->addSphere({1.75f, 3.f, 5.f}, 0.75f, matCT_GraySmoothPlastic);

    // load bunny mesh AS AN AREA LIGHT
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {

        // apply transformations to vertices BEFORE adding as area light
        auto transform = glm::mat4(1.0f);
        transform = glm::rotate(transform, glm::radians(180.0f), glm::vec3(0.f, 2.f, 0.f));
        transform = glm::scale(transform, glm::vec3(2.f, 2.f, 2.f));
        transform = glm::translate(transform, glm::vec3(0.f, 0.0f, 0.f));

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
            {0.9f, 0.9f, 0.9f},
            emission,
            intensity
        );

        std::cout << "Emissive bunny added as mesh area light with "
                  << (bunnyIndices.size() / 3) << " triangle lights" << std::endl;

    } else {
        std::cerr << "Failed to load lowpoly_bunny.obj" << std::endl;
    }


    std::vector<Vertex> ceilingLight_tri1 = {
        {{-2.0f, 9.99f, 3.0f}, {0.0f, -1.0f, 0.0f}},  // Front-left
        {{ 2.0f, 9.99f, 3.0f}, {0.0f, -1.0f, 0.0f}},  // Front-right
        {{ 2.0f, 9.99f, 7.0f}, {0.0f, -1.0f, 0.0f}}   // Back-right
    };
    std::vector<uint32_t>  ceilingLight_indices1 = {0, 1, 2};

    std::vector<Vertex> ceilingLight_tri2 = {
        {{-2.0f, 9.99f, 3.0f}, {0.0f, -1.0f, 0.0f}},  // Front-left
        {{ 2.0f, 9.99f, 7.0f}, {0.0f, -1.0f, 0.0f}},  // Back-right
        {{-2.0f, 9.99f, 7.0f}, {0.0f, -1.0f, 0.0f}}   // Back-left
    };
    std::vector<uint32_t>  ceilingLight_indices2 = {0, 1, 2};

    glm::vec3 trianglesEmission(1.0f, 0.0f, 0.0f);
    float trianglesIntensity = 20.0f;  // higher intensity for smaller area

    pScene->addMeshAreaLight(ceilingLight_tri1, ceilingLight_indices1, {0.9f, 0.f, 0.f}, trianglesEmission, trianglesIntensity);
    pScene->addMeshAreaLight(ceilingLight_tri2, ceilingLight_indices2, {0.9f, 0.f, 0.f}, trianglesEmission, trianglesIntensity);
}

void createSceneC_BunnyOnly(SceneManager* pScene)
{
    // Lambert materials for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    const unsigned char matCT_GrayRoughMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 1.f));
    const unsigned char matCT_GraySmoothMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 0.1f));

    // Room planes
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);

    // Just three spheres for testing
    pScene->addSphere({-1.75f, 5.f, 5.f}, 0.75f, matCT_GrayRoughMetal);
    pScene->addSphere({0.f, 5.f, 5.f}, 0.75f, matCT_GraySmoothMetal);
    pScene->addSphere({1.75f, 5.f, 5.f}, 0.75f, matCT_GrayRoughMetal);

    // Load bunny mesh AS ONLY LIGHT SOURCE (NO ceiling lights!)
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {
        auto transform = glm::mat4(1.0f);
        transform = glm::rotate(transform, glm::radians(180.0f), glm::vec3(0.f, 1.f, 0.f));
        transform = glm::scale(transform, glm::vec3(2.f, 2.f, 2.f));

        for (auto& vertex : bunnyVertices) {
            glm::vec4 pos(vertex.position.x, vertex.position.y, vertex.position.z, 1.0f);
            glm::vec4 transformedPos = transform * pos;
            vertex.position.x = transformedPos.x;
            vertex.position.y = transformedPos.y;
            vertex.position.z = transformedPos.z;

            glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(transform)));
            glm::vec3 normal(vertex.normal.x, vertex.normal.y, vertex.normal.z);
            glm::vec3 transformedNormal = glm::normalize(normalMatrix * normal);
            vertex.normal.x = transformedNormal.x;
            vertex.normal.y = transformedNormal.y;
            vertex.normal.z = transformedNormal.z;
        }

        // INCREASED intensity - bunny is ONLY light
        glm::vec3 emission(1.0f, 1.0f, 1.0f);
        float intensity = 50.0f;  // INCREASED from 10 to 50!

        pScene->addMeshAreaLight(
            bunnyVertices,
            bunnyIndices,
            {0.9f, 0.9f, 0.9f},
            emission,
            intensity
        );

        std::cout << "Scene: Bunny as ONLY light (" << (bunnyIndices.size() / 3) << " triangles)" << std::endl;
    }
}

void createSceneD_FluxStressTest(SceneManager* pScene)
{
    // Camera is at (0, 3, -9) looking forward
    // We need visible objects and visible lights

    // --- Materials ---
    const unsigned char matLambert_Gray = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.7f, 0.7f, 0.7f), 1.f));

    const unsigned char matCT_Shiny = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.9f, 0.9f, 0.9f), 0.0f, 0.2f));

    const unsigned char matCT_Rough = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.7f, 0.7f, 0.7f), 0.0f, 0.8f));

    // --- Room (same as SceneC) ---
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_Gray);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_Gray);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_Gray);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_Gray);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_Gray);    // LEFT

    // --- Test objects (visible from camera) ---
    std::vector<Vertex> lightVerts;
    std::vector<uint32_t> lightIndices;

    // Create 100 triangles in a grid
    int gridSize = 10;
    float cellSize = 0.8f;
    float yPos = 9.8f;
    float zStart = 4.0f;

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            float x0 = -4.0f + i * cellSize;
            float z0 = zStart + j * cellSize;
            float x1 = x0 + cellSize;
            float z1 = z0 + cellSize;

            // Two triangles per cell
            uint32_t base = static_cast<uint32_t>(lightVerts.size());

            lightVerts.push_back({{x0, yPos, z0}, {0, -1, 0}});
            lightVerts.push_back({{x1, yPos, z0}, {0, -1, 0}});
            lightVerts.push_back({{x1, yPos, z1}, {0, -1, 0}});
            lightVerts.push_back({{x0, yPos, z1}, {0, -1, 0}});

            lightIndices.push_back(base);
            lightIndices.push_back(base + 1);
            lightIndices.push_back(base + 2);

            lightIndices.push_back(base);
            lightIndices.push_back(base + 2);
            lightIndices.push_back(base + 3);
        }
    }

    // Add as mesh light
    unsigned int lightIndex = pScene->addMeshAreaLight(
        lightVerts,
        lightIndices,
        {1.0f, 1.0f, 1.0f},
        glm::vec3(1.0f, 1.0f, 1.0f),
        1.0f  // Base intensity
    );

    MeshAreaLight* light = pScene->getMeshAreaLight(lightIndex);
    int totalTris = lightIndices.size() / 3;

    // Make EXTREME flux variations:
    for (int i = 0; i < totalTris; ++i) {
        if (i == 5) {
            light->setTriangleIntensity(i, 10000.0f);  // ONE super bright triangle
        } else if (i < 20) {
            light->setTriangleIntensity(i, 100.0f);    // Some medium ones
        } else {
            light->setTriangleIntensity(i, 0.01f);     // Most very dim
        }
    }

    light->updateBVH();

    std::cout << "\n=== FLUX STRESS TEST ===\n";
    std::cout << "Total triangles: " << totalTris << "\n";
    std::cout << "Triangle 5 intensity: 10000.0 (super bright)\n";
    std::cout << "Triangles 0-19 intensity: 100.0 (medium)\n";
    std::cout << "Other triangles intensity: 0.01 (very dim)\n";
    std::cout << "=========================\n";

    pScene->addSphere({-2.0f, 3.0f, 6.0f}, 0.8f, matCT_Shiny);
    pScene->addSphere({0.0f, 3.0f, 6.0f}, 0.8f, matCT_Shiny);
    pScene->addSphere({2.0f, 3.0f, 6.0f}, 0.8f, matCT_Shiny);
}


void createSimpleTest(SceneManager* pScene) {
    // just 3 triangles with DIFFERENT fluxes
    std::vector<Vertex> verts = {
        {{-1, 9, 4}, {0,-1,0}},
        {{ 1, 9, 4}, {0,-1,0}},
        {{ 0, 9, 6}, {0,-1,0}},

        {{ 2, 9, 4}, {0,-1,0}},
        {{ 4, 9, 4}, {0,-1,0}},
        {{ 3, 9, 6}, {0,-1,0}}
    };

    std::vector<uint32_t> indices = {0,1,2, 3,4,5};

    unsigned int lightIndex = pScene->addMeshAreaLight(
        verts,
        indices,
        {1,1,1},
        glm::vec3(1,1,1),
        1.0f
    );

    MeshAreaLight* light = pScene->getMeshAreaLight(lightIndex);
    // make triangle 0 BRIGHT, triangle 1 dim
    light->setTriangleIntensity(0, 1000.0f);
    light->setTriangleIntensity(1, 0.1f);

    std::cout << "Simple test: 2 triangles, one bright (1000), one dim (0.1)" << std::endl;
}

void setupScene(SceneManager& scene) {
    // calculate scene bounds
    glm::vec3 sceneMin(FLT_MAX);
    glm::vec3 sceneMax(-FLT_MAX);

    // compute from all meshes
    for (int i = 0; i < scene.getNumMeshes(); ++i) {
        Mesh* mesh = scene.getMesh(i);
        if (!mesh) continue;

        const auto& vertices = mesh->getOriginalVertices();
        for (const auto& v : vertices) {
            glm::vec3 pos(v.position.x, v.position.y, v.position.z);
            sceneMin = glm::min(sceneMin, pos);
            sceneMax = glm::max(sceneMax, pos);
        }
    }

    // check all lights
    const auto& lights = scene.getLights();
    for (const auto& light : lights) {
        if (light->type == LightType::MeshArea && light->meshAreaLight) {
            // mesh area lights are already included in meshes
            continue;
        }
        sceneMin = glm::min(sceneMin, light->origin);
        sceneMax = glm::max(sceneMax, light->origin);
    }

    // add padding (10% of scene size)
    glm::vec3 extent = sceneMax - sceneMin;
    glm::vec3 padding = extent * 0.1f;
    sceneMin -= padding;
    sceneMax += padding;

    std::cout << "Scene bounds: " << sceneMin.x << "," << sceneMin.y << "," << sceneMin.z
              << " to " << sceneMax.x << "," << sceneMax.y << "," << sceneMax.z << std::endl;

    // init visibility cache
    MeshAreaLight::initializeVisibilityCache(sceneMin, sceneMax, 16);

    std::cout << "Visibility cache initialized for scene" << std::endl;
}

void printVisibilityStats() {
    int totalCells, activeCells, totalSamples;
    MeshAreaLight::getVisibilityStats(totalCells, activeCells, totalSamples);

    std::cout << "Visibility Cache Stats:" << std::endl;
    std::cout << "  Total cells: " << totalCells << std::endl;
    std::cout << "  Active cells: " << activeCells << " ("
              << (100.0f * activeCells / totalCells) << "%)" << std::endl;
    std::cout << "  Total samples: " << totalSamples << std::endl;
}


int main(int argc, char* argv[])
{
    constexpr uint32_t WIDTH = 640;
    constexpr uint32_t HEIGHT = 480;

    bool testMode = false;
    std::string strategyStr = "vis";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--test-samples" || arg == "-t") {
            testMode = true;
        } else if (arg == "--strategy" || arg == "-s") {
            if (i + 1 < argc) {
                strategyStr = argv[++i];
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -t, --test-samples           Run sampling test" << std::endl;
            std::cout << "  -s, --strategy <name>        uniform, areaimportance, hierarchical" << std::endl;
            std::cout << "  -h, --help                   Show this help" << std::endl;
            return 0;
        }
    }

    SamplingStrategy strategy;
    if (strategyStr == "uniform") {
        strategy = SamplingStrategy::Uniform;
    } else if (strategyStr == "areaimportance" || strategyStr == "area") {
        strategy = SamplingStrategy::AreaImportance;
    } else if (strategyStr == "hierarchical" || strategyStr == "hier") {
        strategy = SamplingStrategy::HierarchicalFlux;
    } else if (strategyStr == "visibilityaware" || strategyStr == "vis") {
        strategy = SamplingStrategy::VisibilityAwareHierarchical;
    } else {
        std::cerr << "Unknown strategy: " << strategyStr << std::endl;
        return -1;
    }

    auto pRenderer = std::make_unique<Renderer>(WIDTH, HEIGHT);
    if (!pRenderer->initialize()) {
        return -1;
    }

    auto pScene = std::make_unique<SceneManager>();
    auto pCamera = std::make_unique<Camera>(
        glm::vec3{0.f, 3.f, -9.f}, 45.f,
        static_cast<float>(WIDTH) / static_cast<float>(HEIGHT)
    );

    createSceneC(pScene.get());
    pScene->commit();
    setupScene(*pScene);

    // CRITICAL: Set strategy AFTER scene creation
    pRenderer->setAreaLightStrategy(*pScene, strategy);

    if (testMode) {
        std::cout << "\n==================================" << std::endl;
        std::cout << "Testing: " << strategyStr << std::endl;
        std::cout << "==================================" << std::endl;
        pRenderer->setTestMode(true);
    }

    auto pTimer = std::make_unique<Timer>();
    pTimer->reset();
    pTimer->start();

    while (!pRenderer->shouldQuit()) {
        pScene->update(pTimer.get());
        pRenderer->render(*pCamera, *pScene);
        pRenderer->present();
        pTimer->update();
    }

    return 0;
}