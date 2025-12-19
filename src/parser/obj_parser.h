//
// Created by ivans on 19/12/2025.
//

#ifndef PROMETHEUS_OBJ_PARSER_H
#define PROMETHEUS_OBJ_PARSER_H

#include <string>
#include <vector>
#include <cstdint>
#include "vertex.h"

bool ParseOBJ(
    const std::string& filename,
    std::vector<Vertex>& outVertices,
    std::vector<uint32_t>& outIndices
);

#endif //PROMETHEUS_OBJ_PARSER_H