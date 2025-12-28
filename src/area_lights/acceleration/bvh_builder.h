//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_BVH_BUILDER_H
#define PROMETHEUS_BVH_BUILDER_H
#include "bvh_config.h"
#include "bvh_node.h"
#include "triangle_build_info.h"
#include <vector>

/// BVH builder class
// constructs hierarchical structure for efficient light sampling

// algorithm:
// 1. compute bounding boxes and centroids for all triangles
// 2. recursively partition triangles using spatial median split
// 3. choose split axis based on longest bounding box dimension
// 4. accumulate flux values bottom-up after construction
class BVHBuilder {
public:
    explicit BVHBuilder(const BVHConfig& config = BVHConfig())
        : m_config(config)
    {}

    // build BVH from triangle information
    // returns: root node index in nodes vector
    int build(
        std::vector<BVHNode>& nodes,
        std::vector<TriangleBuildInfo>& triangles);

    // print tree statistics for debugging
    static void printStatistics(
        const std::vector<BVHNode>& nodes,
        int rootIndex);

private:
    BVHConfig m_config;

    // recursive BVH construction
    int buildNode(
        std::vector<BVHNode>& nodes,
        std::vector<TriangleBuildInfo>& triangles,
        int start, int end,
        int depth);

    // compute flux for each node bottom-up
    static void calculateNodeFlux(
        std::vector<BVHNode>& nodes,
        int nodeIndex);

    // choose best split axis based on bounding box extents
    static int chooseSplitAxis(const glm::vec3& diagonal) ;
};



#endif //PROMETHEUS_BVH_BUILDER_H