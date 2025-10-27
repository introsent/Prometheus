//
// Created by minaj on 10/27/2025.
//

#ifndef HIT_RESULT_H
#define HIT_RESULT_H
#include "embree4/rtcore_common.h"
#include <glm/glm.hpp>

struct HitResult {
    bool didHit = false;
    float distance = 0.f;
    unsigned geomID = RTC_INVALID_GEOMETRY_ID;
    unsigned primID = 0;
    float u = 0.f, v = 0.f;
    glm::vec3 origin = glm::vec3(0.f);
};

#endif //HIT_RESULT_H
