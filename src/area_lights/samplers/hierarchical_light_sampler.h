//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_HIERARCHICAL_LIGHT_SAMPLER_H
#define PROMETHEUS_HIERARCHICAL_LIGHT_SAMPLER_H
#include <glm/vec3.hpp>
#include "acceleration/bvh_node.h"


/// Hierarchical sampler
// uses BVH for importance sampling based on flux and distance

// ref: "Importance Sampling of Many Lights" (Estevez & Kulla, 2018)
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
    // importance = flux / distance^2
    static float calculateImportance(
        const BVHNode& node,
        const glm::vec3& shadingPoint);
};


#endif //PROMETHEUS_HIERARCHICAL_LIGHT_SAMPLER_H