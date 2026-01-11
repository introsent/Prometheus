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
#include "mis/mis_validation.h"
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


// ============================================================================
// SCENE 1: QuadScene - Simple baseline for all sampling methods
// ============================================================================
void createQuadScene(SceneManager* pScene)
{
    // Lambert material for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    // Cornell box walls
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);    // LEFT


    // Create BOTH front and back faces
    std::vector<Vertex> quadVerts = {
        // Front face (facing camera at -Z)
        {{-1.5f, 4.0f, 3.0f}, {0.0f, 0.0f, -1.0f}},  // 0: Bottom-left
        {{ 1.5f, 4.0f, 3.0f}, {0.0f, 0.0f, -1.0f}},  // 1: Bottom-right
        {{ 1.5f, 6.0f, 3.0f}, {0.0f, 0.0f, -1.0f}},  // 2: Top-right
        {{-1.5f, 6.0f, 3.0f}, {0.0f, 0.0f, -1.0f}},  // 3: Top-left

        // Back face (facing spheres at +Z) - DUPLICATE vertices with opposite normal
        {{-1.5f, 4.0f, 3.0f}, {0.0f, 0.0f, 1.0f}},   // 4: Bottom-left (back)
        {{ 1.5f, 4.0f, 3.0f}, {0.0f, 0.0f, 1.0f}},   // 5: Bottom-right (back)
        {{ 1.5f, 6.0f, 3.0f}, {0.0f, 0.0f, 1.0f}},   // 6: Top-right (back)
        {{-1.5f, 6.0f, 3.0f}, {0.0f, 0.0f, 1.0f}}    // 7: Top-left (back)
    };

    std::vector<uint32_t> quadIndices = {
        // Front face (visible from camera)
        0, 1, 2,
        0, 2, 3,

        // Back face (illuminates spheres) - REVERSED winding order
        4, 6, 5,
        4, 7, 6
    };

    // White emission with moderate intensity
    glm::vec3 emission(1.0f, 1.0f, 1.0f);
    float intensity = 50.0f;

    pScene->addMeshAreaLight(
        quadVerts,
        quadIndices,
        {1.0f, 1.0f, 1.0f},  // White albedo
        emission,
        intensity
    );

    std::cout << "\n=== QUAD SCENE ===\n";
    std::cout << "Purpose: Baseline test for sampling method comparison\n";
    std::cout << "Light: Single DOUBLE-SIDED quad (4 triangles total)\n";
    std::cout << "  - Front face visible to camera\n";
    std::cout << "  - Back face illuminates spheres\n";
    std::cout << "Test objects: 3 spheres behind the light\n";
    std::cout << "Expected: All methods should perform similarly (simple case)\n";
    std::cout << "==================\n";
}

// ============================================================================
// SCENE 2: BunnyScene - Complex mesh light test
// ============================================================================
void createBunnyScene(SceneManager* pScene)
{
    // Lambert material for room
    const unsigned char matLambert_GrayBlue = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.49f, 0.57f, 0.57f), 1.f));

    // Cook-Torrance materials for test spheres
    const unsigned char matCT_RoughMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 0.8f));
    const unsigned char matCT_SmoothMetal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.972f, 0.960f, 0.915f), 1.f, 0.2f));

    // Cornell box walls
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_GrayBlue);
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_GrayBlue);
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_GrayBlue);
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_GrayBlue);
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_GrayBlue);

    // Three test spheres to show bunny lighting effect
    pScene->addSphere({-2.0f, 4.5f, 5.5f}, 0.75f, matCT_RoughMetal);
    pScene->addSphere({0.0f, 4.5f, 5.5f}, 0.75f, matCT_SmoothMetal);
    pScene->addSphere({2.0f, 4.5f, 5.5f}, 0.75f, matCT_RoughMetal);

    // Load bunny mesh as ONLY light source
    std::vector<Vertex> bunnyVertices;
    std::vector<uint32_t> bunnyIndices;

    if (ParseOBJ("resources/lowpoly_bunny.obj", bunnyVertices, bunnyIndices)) {
        // Transform bunny to appropriate position and orientation
        auto transform = glm::mat4(1.0f);// Above spheres
        transform = glm::rotate(transform, glm::radians(180.0f), glm::vec3(0.f, 1.f, 0.f));
        transform = glm::scale(transform, glm::vec3(1.5f, 1.5f, 1.5f));

        // Apply transformations
        for (auto& vertex : bunnyVertices) {
            glm::vec4 pos(vertex.position, 1.0f);
            glm::vec4 transformedPos = transform * pos;
            vertex.position = glm::vec3(transformedPos);

            glm::mat3 normalMatrix = glm::transpose(glm::inverse(glm::mat3(transform)));
            vertex.normal = glm::normalize(normalMatrix * vertex.normal);
        }

        // White emission with high intensity (bunny is sole light source)
        glm::vec3 emission(1.0f, 1.0f, 1.0f);
        float intensity = 40.0f;  // Bright enough to illuminate scene

        pScene->addMeshAreaLight(
            bunnyVertices,
            bunnyIndices,
            {1.0f, 1.0f, 1.0f},  // White albedo
            emission,
            intensity
        );

        std::cout << "\n=== BUNNY SCENE ===\n";
        std::cout << "Purpose: Test hierarchical sampling on complex mesh light\n";
        std::cout << "Light: Stanford Bunny mesh (" << (bunnyIndices.size() / 3) << " triangles)\n";
        std::cout << "Expected: Hierarchical methods should outperform uniform sampling\n";
        std::cout << "===================\n";
    } else {
        std::cerr << "ERROR: Failed to load lowpoly_bunny.obj\n";
    }
}

// ============================================================================
// SCENE 3: FluxStressScene - Tests importance sampling with extreme flux variation
// ============================================================================
void createFluxStressScene(SceneManager* pScene)
{
    // Materials
    const unsigned char matLambert_Gray = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.7f, 0.7f, 0.7f), 1.f));

    const unsigned char matCT_Shiny = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.9f, 0.9f, 0.9f), 0.0f, 0.2f));

    // Cornell box
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_Gray);
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_Gray);
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_Gray);
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_Gray);
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_Gray);

    // Test objects (visible from camera)
    pScene->addSphere({-2.0f, 3.0f, 6.0f}, 0.8f, matCT_Shiny);
    pScene->addSphere({0.0f, 3.0f, 6.0f}, 0.8f, matCT_Shiny);
    pScene->addSphere({2.0f, 3.0f, 6.0f}, 0.8f, matCT_Shiny);

    // Create grid of triangles on ceiling with extreme flux variation
    std::vector<Vertex> lightVerts;
    std::vector<uint32_t> lightIndices;

    int gridSize = 10;  // 10x10 = 200 triangles total
    float cellSize = 0.8f;
    float yPos = 9.8f;
    float zStart = 4.0f;

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            float x0 = -4.0f + i * cellSize;
            float z0 = zStart + j * cellSize;
            float x1 = x0 + cellSize;
            float z1 = z0 + cellSize;

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
        1.0f  // Base intensity (will be overridden per triangle)
    );

    MeshAreaLight* light = pScene->getMeshAreaLight(lightIndex);
    int totalTris = lightIndices.size() / 3;

    // Create EXTREME flux distribution:
    // - 1 super bright triangle (10000x base)
    // - 20 medium triangles (100x base)
    // - Rest very dim (0.01x base)
    for (int i = 0; i < totalTris; ++i) {
        if (i == 5) {
            light->setTriangleIntensity(i, 10000.0f);  // ONE hotspot
        } else if (i < 20) {
            light->setTriangleIntensity(i, 100.0f);    // Some medium
        } else {
            light->setTriangleIntensity(i, 0.01f);     // Most very dim
        }
    }

    light->updateBVH();

    std::cout << "\n=== FLUX STRESS SCENE ===\n";
    std::cout << "Purpose: Test flux-based importance sampling\n";
    std::cout << "Light: 200 triangles with extreme flux variation\n";
    std::cout << "  - Triangle 5: 10000x intensity (0.5% of triangles)\n";
    std::cout << "  - Triangles 0-19: 100x intensity (10% of triangles)\n";
    std::cout << "  - Rest: 0.01x intensity (89.5% of triangles)\n";
    std::cout << "Expected: Flux-based BVH should heavily sample triangle 5\n";
    std::cout << "==========================\n";
}

// ============================================================================
// SCENE 4: OcclusionStressScene - Tests visibility-aware sampling
// ============================================================================
void createOcclusionStressScene(SceneManager* pScene)
{
    // Materials
    const unsigned char matLambert_Gray = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.7f, 0.7f, 0.7f), 1.f));

    const unsigned char matLambert_Red = pScene->addMaterial(
        new Material_Lambert(glm::vec3(0.8f, 0.2f, 0.2f), 1.f));

    const unsigned char matCT_Metal = pScene->addMaterial(
        new Material_CookTorrence(glm::vec3(0.9f, 0.9f, 0.9f), 1.0f, 0.3f));

    // Cornell box
    pScene->addPlane({0.f, 0.f, 10.f}, {0.f, 0.f, -1.f}, matLambert_Gray);   // BACK
    pScene->addPlane({0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, matLambert_Gray);     // BOTTOM
    pScene->addPlane({0.f, 10.f, 0.f}, {0.f, -1.f, 0.f}, matLambert_Gray);   // TOP
    pScene->addPlane({5.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, matLambert_Gray);    // RIGHT
    pScene->addPlane({-5.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, matLambert_Gray);    // LEFT

    // ========================================
    // Test objects FIRST (so they're visible)
    // ========================================
    // Objects below occluder (y < 7.5) - should be mostly dark
    pScene->addSphere({-2.5f, 2.5f, 5.0f}, 0.8f, matCT_Metal);
    pScene->addSphere({0.0f, 2.5f, 5.0f}, 0.8f, matCT_Metal);
    pScene->addSphere({2.5f, 2.5f, 5.0f}, 0.8f, matCT_Metal);

    // One sphere ABOVE occluder (y > 7.5) - should be well-lit
    pScene->addSphere({0.0f, 8.5f, 5.0f}, 0.6f, matCT_Metal);

    // ========================================
    // Create OCCLUDER: Large horizontal plane at y=6.0 (lower position)
    // ========================================
    // This creates a "ceiling" at y=6 that blocks light from reaching lower spheres
    // Using a mesh instead of infinite plane so we can see it
    std::vector<Vertex> occluderVerts = {
        {{-4.5f, 6.0f, 2.0f}, {0.0f, 1.0f, 0.0f}},   // Front-left (normal pointing UP)
        {{ 4.5f, 6.0f, 2.0f}, {0.0f, 1.0f, 0.0f}},   // Front-right
        {{ 4.5f, 6.0f, 9.5f}, {0.0f, 1.0f, 0.0f}},   // Back-right
        {{-4.5f, 6.0f, 9.5f}, {0.0f, 1.0f, 0.0f}}    // Back-left
    };

    std::vector<uint32_t> occluderIndices = {
        0, 1, 2,  // First triangle
        0, 2, 3   // Second triangle
    };

    // Add occluder as regular mesh (not a light)
    pScene->addMesh(
        occluderVerts,
        occluderIndices,
        matLambert_Red  // Red color so it's visible
    );

    // ========================================
    // Create LARGE mesh light on ceiling (200 triangles) ABOVE the occluder
    // ========================================
    std::vector<Vertex> ceilingLightVerts;
    std::vector<uint32_t> ceilingLightIndices;

    int gridSize = 10;  // 10x10 grid = 200 triangles
    float cellSize = 0.8f;
    float yPos = 9.8f;   // Near ceiling
    float zStart = 2.0f;

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            float x0 = -4.0f + i * cellSize;
            float z0 = zStart + j * cellSize;
            float x1 = x0 + cellSize;
            float z1 = z0 + cellSize;

            uint32_t base = static_cast<uint32_t>(ceilingLightVerts.size());

            ceilingLightVerts.push_back({{x0, yPos, z0}, {0, -1, 0}});
            ceilingLightVerts.push_back({{x1, yPos, z0}, {0, -1, 0}});
            ceilingLightVerts.push_back({{x1, yPos, z1}, {0, -1, 0}});
            ceilingLightVerts.push_back({{x0, yPos, z1}, {0, -1, 0}});

            ceilingLightIndices.push_back(base);
            ceilingLightIndices.push_back(base + 1);
            ceilingLightIndices.push_back(base + 2);

            ceilingLightIndices.push_back(base);
            ceilingLightIndices.push_back(base + 2);
            ceilingLightIndices.push_back(base + 3);
        }
    }

    pScene->addMeshAreaLight(
        ceilingLightVerts,
        ceilingLightIndices,
        {1.0f, 1.0f, 1.0f},
        glm::vec3(1.0f, 1.0f, 1.0f),
        30.0f  // Higher intensity to compensate for occlusion
    );

    std::cout << "\n=== OCCLUSION STRESS SCENE ===\n";
    std::cout << "Purpose: Test visibility-aware sampling with heavy occlusion\n";
    std::cout << "Setup:\n";
    std::cout << "  - Large ceiling light (200 triangles at y=9.8)\n";
    std::cout << "  - Red occluder plane at y=6.0 blocking direct light\n";
    std::cout << "  - 3 test spheres below occluder at y=2.5 (should be dark)\n";
    std::cout << "  - 1 test sphere above occluder at y=8.5 (should be bright)\n";
    std::cout << "Expected:\n";
    std::cout << "  - Uniform sampling: wastes samples on occluded triangles\n";
    std::cout << "  - Area/Flux BVH: still samples many occluded triangles\n";
    std::cout << "  - Visibility-aware MIS: learns to avoid occluded regions\n";
    std::cout << "  - Upper sphere should be much brighter than lower spheres\n";
    std::cout << "===============================\n";
}

// ============================================================================
// Scene selection helper
// ============================================================================
void createTestScene(SceneManager* pScene, const std::string& sceneName)
{
    if (sceneName == "quad") {
        createQuadScene(pScene);
    }
    else if (sceneName == "bunny") {
        createBunnyScene(pScene);
    }
    else if (sceneName == "flux") {
        createFluxStressScene(pScene);
    }
    else if (sceneName == "occlusion") {
        createOcclusionStressScene(pScene);
    }
    else {
        std::cerr << "Unknown scene: " << sceneName << "\n";
        std::cerr << "Available scenes: quad, bunny, flux, occlusion\n";
    }
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

    bool genGroundTruth = false;
    bool testMode = false;
    std::string strategyStr = "vis";
    std::string sceneName = "bunny"; // Default scene

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--gt") {
            genGroundTruth = true;
        } else if (arg == "--test-samples" || arg == "-t") {
            testMode = true;
        } else if (arg == "--strategy" || arg == "-s") {
            if (i + 1 < argc) {
                strategyStr = argv[++i];
            }
        } else if (arg == "--scene" || arg == "-sc") {
            if (i + 1 < argc) {
                sceneName = argv[++i];
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -t, --test-samples           Run sampling test" << std::endl;
            std::cout << "  -s, --strategy <name>        uniform, areaimportance, hierarchical, visibilityaware" << std::endl;
            std::cout << "  -sc, --scene <name>          flux, bunny, quad, occlusion (default: quad)" << std::endl;
            std::cout << "  --gt                         Generate ground truth" << std::endl;
            std::cout << "  -h, --help                   Show this help" << std::endl;
            return 0;
        }
    }

    // Create scene based on argument
    auto pScene = std::make_unique<SceneManager>();
    std::cout << "Loading scene: " << sceneName << std::endl;

    if (sceneName == "flux" || sceneName == "fluxstress") {
        createFluxStressScene(pScene.get());
        std::cout << "Created Flux Stress Scene" << std::endl;
    } else if (sceneName == "bunny") {
        createBunnyScene(pScene.get());
        std::cout << "Created Bunny Scene" << std::endl;
    } else if (sceneName == "quad") {
        createQuadScene(pScene.get());
        std::cout << "Created Quad Scene" << std::endl;
    } else if (sceneName == "occlusion" || sceneName == "occlusionstress") {
        createOcclusionStressScene(pScene.get());
        std::cout << "Created Occlusion Stress Scene" << std::endl;
    } else {
        std::cerr << "Unknown scene: " << sceneName << std::endl;
        std::cerr << "Available scenes: flux, bunny, quad, occlusion" << std::endl;
        return -1;
    }

    pScene->commit();
    setupScene(*pScene);

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

    // Pass scene name to renderer
    pRenderer->setSceneName(sceneName);

    auto pCamera = std::make_unique<Camera>(
        glm::vec3{0.f, 3.f, -9.f}, 45.f,
        static_cast<float>(WIDTH) / static_cast<float>(HEIGHT)
    );

    // CRITICAL: Set strategy AFTER scene creation
    pRenderer->setAreaLightStrategy(*pScene, strategy);

    if (genGroundTruth) {
        // Generate ground truth with scene-specific path
        std::string gtDir = "tests/" + sceneName + "/ground_truth";
        std::cout << "Generating canonical ground truth into: " << gtDir << std::endl;
        if (sceneName == "flux") {
            pRenderer->generateGroundTruth(*pCamera, *pScene,
                                       SamplingStrategy::HierarchicalFlux,
                                       /*samples=*/2000,
                                       gtDir);
        } else if (sceneName == "occlusion") {
            pRenderer->generateGroundTruth(*pCamera, *pScene,
                                      SamplingStrategy::VisibilityAwareHierarchical,
                                      /*samples=*/2000,
                                      gtDir);
        }
        else {
            pRenderer->generateGroundTruth(*pCamera, *pScene,
                                      SamplingStrategy::Uniform,
                                      /*samples=*/20000,
                                      gtDir);
        }

        return 0; // exit after GT
    }

    //MISValidation::runAllTests();

    if (testMode) {
        std::cout << "\n==================================" << std::endl;
        std::cout << "Scene: " << sceneName << std::endl;
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