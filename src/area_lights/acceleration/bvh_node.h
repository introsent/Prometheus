//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_BVH_NODE_H
#define PROMETHEUS_BVH_NODE_H
#include <vector>
#include <glm/vec3.hpp>


/// BVH node structure
// stores bounding box, flux, and child information
struct BVHNode {
    // bounding box
    glm::vec3 bboxMin{0.0f};
    glm::vec3 bboxMax{0.0f};

    // accumulated light properties
    float totalFlux{0.0f};   // sum of flux from all triangles in subtree
    float totalArea{0.0f};   // sum of triangle areas in subtree

    // tree structure
    int leftChild{-1};       // index of left child (-1 if leaf)
    int rightChild{-1};      // index of right child (-1 if leaf)
    bool isLeaf{false};

    // leaf node data
    int startTri{0};         // first triangle index (for construction)
    int endTri{0};           // one past last triangle index
    std::vector<int> triangleIndices;  // actual triangle indices in leaf

    BVHNode() = default;

    // get centroid of bounding box
    [[nodiscard]] glm::vec3 getCentroid() const {
        return (bboxMin + bboxMax) * 0.5f;
    }

    // get diagonal vector of bounding box
    [[nodiscard]] glm::vec3 getDiagonal() const {
        return bboxMax - bboxMin;
    }
};

#endif //PROMETHEUS_BVH_NODE_H