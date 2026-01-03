//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_BVH_BUILDER_H
#define PROMETHEUS_BVH_BUILDER_H

#include "bvh_config.h"
#include "bvh_node.h"
#include "triangle_build_info.h"
#include <vector>

class BVHBuilder {
public:
    explicit BVHBuilder(const BVHConfig& config = BVHConfig())
        : m_config(config)
    {}

    int build(
        std::vector<BVHNode>& nodes,
        std::vector<TriangleBuildInfo>& triangles);

    static void printStatistics(
        const std::vector<BVHNode>& nodes,
        int rootIndex);

private:
    BVHConfig m_config;

    int buildNode(
        std::vector<BVHNode>& nodes,
        std::vector<TriangleBuildInfo>& triangles,
        int start, int end,
        int depth);

    // RENAMED: Now calculates flux, orientation bounds, AND variance
    static void calculateNodeProperties(
        std::vector<BVHNode>& nodes,
        int nodeIndex);

    static int chooseSplitAxis(const glm::vec3& diagonal);
};

#endif