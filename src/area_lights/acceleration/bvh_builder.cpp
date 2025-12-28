//
// Created by ivans on 28/12/2025.
//

#include "bvh_builder.h"
#include "math_helpers.h"
#include <algorithm>
#include <iostream>
#include <cfloat>
#include <functional>

int BVHBuilder::build(
    std::vector<BVHNode>& nodes,
    std::vector<TriangleBuildInfo>& triangles) {

    nodes.clear();

    if (triangles.empty()) {
        return -1;
    }

    std::cout << "Building BVH for " << triangles.size() << " triangles...\n";

    // build tree recursively
    const int rootIndex = buildNode(nodes, triangles, 0,
        static_cast<int>(triangles.size()), 0);

    // calculate flux values bottom-up
    calculateNodeFlux(nodes, rootIndex);

    printStatistics(nodes, rootIndex);

    return rootIndex;
}

int BVHBuilder::buildNode(
    std::vector<BVHNode>& nodes,
    std::vector<TriangleBuildInfo>& triangles,
    const int start, const int end,
    const int depth) {

    // create new node
    const int nodeIndex = static_cast<int>(nodes.size());
    nodes.emplace_back();
    BVHNode& node = nodes[nodeIndex];

    // compute bounding box for all triangles in this node
    glm::vec3 bboxMin(FLT_MAX);
    glm::vec3 bboxMax(-FLT_MAX);

    for (int i = start; i < end; ++i) {
        MathHelpers::expandBoundingBox(bboxMin, bboxMax,
            triangles[i].bboxMin);
        MathHelpers::expandBoundingBox(bboxMin, bboxMax,
            triangles[i].bboxMax);
    }

    node.bboxMin = bboxMin;
    node.bboxMax = bboxMax;

    const int numTriangles = end - start;

    // check termination criteria
    if (numTriangles <= m_config.maxTrianglesPerLeaf ||
        depth >= m_config.maxDepth) {

        // create leaf node
        node.isLeaf = true;
        node.startTri = start;
        node.endTri = end;
        node.totalArea = 0.0f;

        // store original triangle indices
        node.triangleIndices.reserve(numTriangles);
        for (int i = start; i < end; ++i) {
            node.triangleIndices.push_back(triangles[i].originalIndex);
            node.totalArea += triangles[i].area;
        }

        return nodeIndex;
    }

    // choose split axis (longest bounding box dimension)
    // this follows the Surface Area Heuristic (SAH) principle
    const glm::vec3 diagonal = node.getDiagonal();
    const int splitAxis = chooseSplitAxis(diagonal);

    // sort triangles by centroid along split axis
    // this implements spatial median split
    std::sort(
        triangles.begin() + start,
        triangles.begin() + end,
        [splitAxis](const TriangleBuildInfo& a, const TriangleBuildInfo& b) {
            return a.centroid[splitAxis] < b.centroid[splitAxis];
        }
    );

    // split at median
    const int mid = start + numTriangles / 2;

    // recursively build children
    const int leftChild = buildNode(nodes, triangles, start, mid, depth + 1);
    const int rightChild = buildNode(nodes, triangles, mid, end, depth + 1);

    // update node with child information
    // note: must access node again as vector may have reallocated
    nodes[nodeIndex].leftChild = leftChild;
    nodes[nodeIndex].rightChild = rightChild;
    nodes[nodeIndex].isLeaf = false;

    return nodeIndex;
}

void BVHBuilder::calculateNodeFlux(
    std::vector<BVHNode>& nodes,
    int nodeIndex) {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(nodes.size())) {
        return;
    }

    BVHNode& node = nodes[nodeIndex];

    if (node.isLeaf) {
        // leaf nodes have their flux computed during construction
        // this is done by summing triangle fluxes
        // (actual implementation depends on triangle data structure)
        return;
    }

    // internal node: recursively calculate children first
    calculateNodeFlux(nodes, node.leftChild);
    calculateNodeFlux(nodes, node.rightChild);

    // validate child indices
    if (node.leftChild < 0 || node.rightChild < 0 ||
        node.leftChild >= static_cast<int>(nodes.size()) ||
        node.rightChild >= static_cast<int>(nodes.size())) {

        node.totalFlux = 0.0f;
        node.totalArea = 0.0f;
        std::cerr << "ERROR: Invalid child indices for node " << nodeIndex << "\n";
        return;
    }

    // sum children's flux and area
    const BVHNode& left = nodes[node.leftChild];
    const BVHNode& right = nodes[node.rightChild];

    node.totalFlux = left.totalFlux + right.totalFlux;
    node.totalArea = left.totalArea + right.totalArea;
}

int BVHBuilder::chooseSplitAxis(const glm::vec3& diagonal) {
    // choose axis with largest extent
    // this minimizes overlap in child bounding boxes

    int axis = 0;  // x-axis

    if (diagonal.y > diagonal.x) {
        axis = 1;  // y-axis
    }

    if (diagonal.z > diagonal[axis]) {
        axis = 2;  // z-axis
    }

    return axis;
}

void BVHBuilder::printStatistics(
    const std::vector<BVHNode>& nodes,
    int rootIndex) {

    if (rootIndex < 0 || rootIndex >= static_cast<int>(nodes.size())) {
        std::cout << "Invalid BVH (empty or invalid root)\n";
        return;
    }

    int leafCount = 0;
    int maxDepth = 0;
    const float totalFlux = nodes[rootIndex].totalFlux;

    // traverse tree to collect statistics
    std::function<void(int, int)> traverse = [&](int index, int depth) {
        if (index < 0 || index >= static_cast<int>(nodes.size())) {
            return;
        }

        maxDepth = std::max(maxDepth, depth);

        const BVHNode& node = nodes[index];

        if (node.isLeaf) {
            leafCount++;
        } else {
            traverse(node.leftChild, depth + 1);
            traverse(node.rightChild, depth + 1);
        }
    };

    traverse(rootIndex, 0);

    std::cout << "BVH Statistics:\n"
              << "  Total nodes: " << nodes.size() << "\n"
              << "  Leaf nodes: " << leafCount << "\n"
              << "  Internal nodes: " << (nodes.size() - leafCount) << "\n"
              << "  Max depth: " << maxDepth << "\n"
              << "  Total flux: " << totalFlux << "\n";
}

// ============================================================================
// hierarchical sampler
// uses BVH for importance sampling based on flux and distance
//
// references:
// - "Importance Sampling of Many Lights" (Estevez & Kulla, 2018)
//   importance = flux / distance²
// ============================================================================

class HierarchicalLightSampler {
public:
    // select a leaf node using hierarchical importance sampling
    // returns: leaf node index, accumulated pdf for the path
    static int selectLeafNode(
        const std::vector<BVHNode>& nodes,
        int rootIndex,
        const glm::vec3& shadingPoint,
        float randomNumber,
        float& outPathPdf);

private:
    // recursively traverse tree, choosing child based on importance
    static int selectNodeRecursive(
        const std::vector<BVHNode>& nodes,
        int nodeIndex,
        const glm::vec3& shadingPoint,
        float u,
        float& pdf);

    // calculate importance of a node from shading point
    // importance = flux / distance²
    static float calculateImportance(
        const BVHNode& node,
        const glm::vec3& shadingPoint);
};

int HierarchicalLightSampler::selectLeafNode(
    const std::vector<BVHNode>& nodes,
    int rootIndex,
    const glm::vec3& shadingPoint,
    float randomNumber,
    float& outPathPdf) {

    outPathPdf = 1.0f;

    if (rootIndex < 0 || rootIndex >= static_cast<int>(nodes.size())) {
        return -1;
    }

    return selectNodeRecursive(nodes, rootIndex, shadingPoint,
        randomNumber, outPathPdf);
}

int HierarchicalLightSampler::selectNodeRecursive(
    const std::vector<BVHNode>& nodes,
    int nodeIndex,
    const glm::vec3& shadingPoint,
    float u,
    float& pdf) {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(nodes.size())) {
        return -1;
    }

    const BVHNode& node = nodes[nodeIndex];

    // base case: reached leaf
    if (node.isLeaf) {
        return nodeIndex;
    }

    // get children
    const BVHNode& left = nodes[node.leftChild];
    const BVHNode& right = nodes[node.rightChild];

    // calculate importance for each child
    // importance combines flux with inverse square distance falloff
    const float importanceLeft = calculateImportance(left, shadingPoint);
    const float importanceRight = calculateImportance(right, shadingPoint);

    const float totalImportance = importanceLeft + importanceRight;

    if (totalImportance <= 0.0f) {
        return -1;  // no valid children
    }

    // probability of choosing left child
    const float probLeft = importanceLeft / totalImportance;

    // make stochastic decision and recurse
    if (u < probLeft) {
        pdf *= probLeft;
        const float newU = (probLeft > 0.0f) ? u / probLeft : 0.0f;
        return selectNodeRecursive(nodes, node.leftChild, shadingPoint,
            newU, pdf);
    } else {
        const float probRight = 1.0f - probLeft;
        pdf *= probRight;
        const float newU = (probRight > 0.0f) ? (u - probLeft) / probRight : 0.0f;
        return selectNodeRecursive(nodes, node.rightChild, shadingPoint,
            newU, pdf);
    }
}

float HierarchicalLightSampler::calculateImportance(
    const BVHNode& node,
    const glm::vec3& shadingPoint) {

    // compute distance from shading point to node centroid
    const glm::vec3 toNode = node.getCentroid() - shadingPoint;
    float distSq = glm::dot(toNode, toNode);

    // add small epsilon to prevent division by zero
    distSq = std::max(distSq, 1e-4f);

    // importance = flux / distance²
    // this approximates the contribution of the light cluster
    // following inverse square law
    return node.totalFlux / distSq;
}