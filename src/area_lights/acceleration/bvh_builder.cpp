#include "bvh_builder.h"
#include "math_helpers.h"
#include <algorithm>
#include <iostream>
#include <cfloat>
#include <functional>
#include <cmath>

int BVHBuilder::build(
    std::vector<BVHNode>& nodes,
    std::vector<TriangleBuildInfo>& triangles) {

    nodes.clear();

    if (triangles.empty()) {
        return -1;
    }

    std::cout << "Building BVH for " << triangles.size() << " triangles...\n";

    const int rootIndex = buildNode(nodes, triangles, 0,
        static_cast<int>(triangles.size()), 0);

    // calculates flux, orientation bounds, and variance
    calculateNodeProperties(nodes, rootIndex);

    printStatistics(nodes, rootIndex);

    return rootIndex;
}

int BVHBuilder::buildNode(
    std::vector<BVHNode>& nodes,
    std::vector<TriangleBuildInfo>& triangles,
    const int start, const int end,
    const int depth) {

    const int nodeIndex = static_cast<int>(nodes.size());
    nodes.emplace_back();
    BVHNode& node = nodes[nodeIndex];

    // compute bounding box
    glm::vec3 bboxMin(FLT_MAX);
    glm::vec3 bboxMax(-FLT_MAX);

    for (int i = start; i < end; ++i) {
        MathHelpers::expandBoundingBox(bboxMin, bboxMax, triangles[i].bboxMin);
        MathHelpers::expandBoundingBox(bboxMin, bboxMax, triangles[i].bboxMax);
    }

    node.bboxMin = bboxMin;
    node.bboxMax = bboxMax;

    const int numTriangles = end - start;

    // leaf node creation
    if (numTriangles <= m_config.maxTrianglesPerLeaf || depth >= m_config.maxDepth) {
        node.isLeaf = true;
        node.startTri = start;
        node.endTri = end;
        node.totalArea = 0.0f;
        node.totalFlux = 0.0f;
        node.fluxVariance = 0.0f;

        // initialize orientation bounds for leaf
        bool firstTriangle = true;
        std::vector<float> fluxValues;
        fluxValues.reserve(numTriangles);

        node.triangleIndices.reserve(numTriangles);
        for (int i = start; i < end; ++i) {
            node.triangleIndices.push_back(triangles[i].originalIndex);
            node.totalArea += triangles[i].area;
            node.totalFlux += triangles[i].flux;
            fluxValues.push_back(triangles[i].flux);

            // build orientation cone from triangle normals
            OrientationCone triCone(triangles[i].normal, static_cast<float>(M_PI) / 2.0f);

            if (firstTriangle) {
                node.orientationBounds = triCone;
                firstTriangle = false;
            } else {
                node.orientationBounds = mergeOrientationCones(
                    node.orientationBounds, triCone);
            }
        }

        // calculate flux variance for adaptive splitting
        if (!fluxValues.empty()) {
            float meanFlux = node.totalFlux / static_cast<float>(fluxValues.size());
            for (float f : fluxValues) {
                float diff = f - meanFlux;
                node.fluxVariance += diff * diff;
            }
            node.fluxVariance /= static_cast<float>(fluxValues.size());
        }

        return nodeIndex;
    }

    // internal node: split and recurse
    const glm::vec3 diagonal = node.getDiagonal();
    const int splitAxis = chooseSplitAxis(diagonal);

    std::sort(
        triangles.begin() + start,
        triangles.begin() + end,
        [splitAxis](const TriangleBuildInfo& a, const TriangleBuildInfo& b) {
            return a.centroid[splitAxis] < b.centroid[splitAxis];
        }
    );

    const int mid = start + numTriangles / 2;

    const int leftChild = buildNode(nodes, triangles, start, mid, depth + 1);
    const int rightChild = buildNode(nodes, triangles, mid, end, depth + 1);

    nodes[nodeIndex].leftChild = leftChild;
    nodes[nodeIndex].rightChild = rightChild;
    nodes[nodeIndex].isLeaf = false;

    return nodeIndex;
}

void BVHBuilder::calculateNodeProperties(
    std::vector<BVHNode>& nodes,
    int nodeIndex) {

    if (nodeIndex < 0 || nodeIndex >= static_cast<int>(nodes.size())) {
        return;
    }

    BVHNode& node = nodes[nodeIndex];

    if (node.isLeaf) {
        // leaf properties are computed during construction
        return;
    }

    // recursively process children first (bottom-up)
    calculateNodeProperties(nodes, node.leftChild);
    calculateNodeProperties(nodes, node.rightChild);

    // validate child indices
    if (node.leftChild < 0 || node.rightChild < 0 ||
        node.leftChild >= static_cast<int>(nodes.size()) ||
        node.rightChild >= static_cast<int>(nodes.size())) {

        node.totalFlux = 0.0f;
        node.totalArea = 0.0f;
        std::cerr << "ERROR: Invalid child indices for node " << nodeIndex << "\n";
        return;
    }

    const BVHNode& left = nodes[node.leftChild];
    const BVHNode& right = nodes[node.rightChild];

    // sum flux and area from children
    node.totalFlux = left.totalFlux + right.totalFlux;
    node.totalArea = left.totalArea + right.totalArea;

    // merge orientation bounds from children
    node.orientationBounds = mergeOrientationCones(
        left.orientationBounds,
        right.orientationBounds
    );

    // compute flux variance between children for adaptive splitting
    // this measures how unevenly flux is distributed in the subtree
    float meanFlux = node.totalFlux / 2.0f;
    float leftDiff = left.totalFlux - meanFlux;
    float rightDiff = right.totalFlux - meanFlux;

    // variance combines child variance with between-child variance
    // this helps identify clusters with high variance that should be split
    node.fluxVariance = (left.fluxVariance + right.fluxVariance) / 2.0f  // average child variance
                      + (leftDiff * leftDiff + rightDiff * rightDiff);   // between-child variance
}

int BVHBuilder::chooseSplitAxis(const glm::vec3& diagonal) {
    int axis = 0;
    if (diagonal.y > diagonal.x) {
        axis = 1;
    }
    if (diagonal.z > diagonal[axis]) {
        axis = 2;
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

    // track orientation bounds statistics
    float maxThetaO = 0.0f;

    std::function<void(int, int)> traverse = [&](int index, int depth) {
        if (index < 0 || index >= static_cast<int>(nodes.size())) {
            return;
        }

        maxDepth = std::max(maxDepth, depth);

        const BVHNode& node = nodes[index];
        maxThetaO = std::max(maxThetaO, node.orientationBounds.thetaO);

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
              << "  Total flux: " << totalFlux << "\n"
              << "  Root theta_o: " << nodes[rootIndex].orientationBounds.thetaO
              << " rad (" << nodes[rootIndex].orientationBounds.thetaO * 180.0f / M_PI << " deg)\n"
              << "  Root flux variance: " << nodes[rootIndex].fluxVariance << "\n";
}