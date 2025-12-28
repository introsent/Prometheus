//
// Created by ivans on 28/12/2025.
//

#include "hierarchical_light_sampler.h"
#include <glm/glm.hpp>

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