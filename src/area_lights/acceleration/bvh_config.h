//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_BVH_CONFIG_H
#define PROMETHEUS_BVH_CONFIG_H

/// BVH construction configuration
struct BVHConfig {
    int maxTrianglesPerLeaf{4};   // leaf creation threshold
    int maxDepth{20};              // maximum tree depth

    // split heuristic weights
    // inspired by Surface Area Orientation Heuristic (SAOH) from Estevez & Kulla
    float spatialWeight{1.0f};     // weight for spatial extent
    float fluxWeight{1.0f};        // weight for flux distribution
};

#endif //PROMETHEUS_BVH_CONFIG_H