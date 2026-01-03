//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_BVH_NODE_H
#define PROMETHEUS_BVH_NODE_H
#include <vector>
#include <glm/glm.hpp>
#include <cmath>
#include <corecrt_math_defines.h>

/// Orientation cone for light clusters
/// Based on Conty & Kulla "Importance Sampling of Many Lights" Section 4.1
struct OrientationCone {
    glm::vec3 axis{0.0f, 1.0f, 0.0f};  // Cone axis direction
    float thetaO{static_cast<float>(M_PI)};                  // Normal bounding angle (0 = flat, π = full sphere)
    float thetaE{static_cast<float>(M_PI / 2.f)};          // Emission profile angle (π/2 for diffuse emitters)

    OrientationCone() = default;

    explicit OrientationCone(const glm::vec3& normal, float emissionAngle = M_PI / 2.0f)
        : axis(glm::normalize(normal))
        , thetaO(0.0f)  // Single flat emitter
        , thetaE(emissionAngle)
    {}

    // Check if this cone covers entire sphere of directions
    [[nodiscard]] bool isFullSphere() const {
        return thetaO >= M_PI;
    }

    // compute orientation measure M_omega
    // this measures the "angular extent" of the cluster
    [[nodiscard]] float computeMeasure() const {
        if (thetaO >= M_PI) {
            return 4.0f * static_cast<float>(M_PI);  // full sphere
        }

        float thetaW = std::min(thetaO + thetaE, static_cast<float>(M_PI));
        float cosThetaO = std::cos(thetaO);
        float sinThetaO = std::sin(thetaO);

        // M_omega = 2pi[(1 - cos theta_o) + integral term]
        // simplified closed form from paper
        float baseAngle = 2.0f * static_cast<float>(M_PI) * (1.0f - cosThetaO);

        if (thetaE > 0.0f && thetaW > thetaO) {
            // add emission profile contribution
            float integralTerm = static_cast<float>(M_PI) / 2.0f * (
                2.0f * thetaW * sinThetaO
                - std::cos(thetaO - 2.0f * thetaW)
                - 2.0f * thetaO * sinThetaO
                + cosThetaO
            );
            baseAngle += integralTerm;
        }

        return std::max(baseAngle, 0.0f);
    }
};

/// merge two orientation cones
inline OrientationCone mergeOrientationCones(
    const OrientationCone& a,
    const OrientationCone& b) {

    OrientationCone result;

    // handle degenerate cases
    if (a.thetaO >= M_PI) return a;
    if (b.thetaO >= M_PI) return b;

    // ensure a has larger theta_o
    const OrientationCone* coneA = &a;
    const OrientationCone* coneB = &b;
    if (b.thetaO > a.thetaO) {
        std::swap(coneA, coneB);
    }

    // angle between axes
    float dotProduct = glm::clamp(glm::dot(coneA->axis, coneB->axis), -1.0f, 1.0f);
    float thetaD = std::acos(dotProduct);

    // new emission angle is max of both
    result.thetaE = std::max(coneA->thetaE, coneB->thetaE);

    // check if cone A already covers cone B
    if (std::min(thetaD + coneB->thetaO, static_cast<float>(M_PI)) <= coneA->thetaO) {
        result.axis = coneA->axis;
        result.thetaO = coneA->thetaO;
        return result;
    }

    // compute new cone covering both
    result.thetaO = (coneA->thetaO + thetaD + coneB->thetaO) / 2.0f;

    if (result.thetaO >= M_PI) {
        result.axis = coneA->axis;
        result.thetaO = static_cast<float>(M_PI);
        return result;
    }

    // rotate axis towards b to center the new cone
    float thetaR = result.thetaO - coneA->thetaO;

    if (thetaD > 1e-6f) {
        // compute rotation axis
        glm::vec3 rotAxis = glm::cross(coneA->axis, coneB->axis);
        float rotAxisLen = glm::length(rotAxis);

        if (rotAxisLen > 1e-6f) {
            rotAxis /= rotAxisLen;

            // Rodrigues rotation formula
            float cosR = std::cos(thetaR);
            float sinR = std::sin(thetaR);

            result.axis = coneA->axis * cosR
                        + glm::cross(rotAxis, coneA->axis) * sinR
                        + rotAxis * glm::dot(rotAxis, coneA->axis) * (1.0f - cosR);
            result.axis = glm::normalize(result.axis);
        } else {
            result.axis = coneA->axis;
        }
    } else {
        result.axis = coneA->axis;
    }

    return result;
}

/// BVH node structure with orientation bounds
struct BVHNode {
    // bounding box
    glm::vec3 bboxMin{0.0f};
    glm::vec3 bboxMax{0.0f};

    // accumulated light properties
    float totalFlux{0.0f};
    float totalArea{0.0f};

    // orientation bounds
    OrientationCone orientationBounds;

    // variance estimate for adaptive splitting
    float fluxVariance{0.0f};  // variance in emitter flux within cluster

    // tree structure
    int leftChild{-1};
    int rightChild{-1};
    bool isLeaf{false};

    // leaf node data
    int startTri{0};
    int endTri{0};
    std::vector<int> triangleIndices;

    BVHNode() = default;

    [[nodiscard]] glm::vec3 getCentroid() const {
        return (bboxMin + bboxMax) * 0.5f;
    }

    [[nodiscard]] glm::vec3 getDiagonal() const {
        return bboxMax - bboxMin;
    }

    // compute bounding sphere radius (for theta_u calculation)
    [[nodiscard]] float getBoundingSphereRadius() const {
        return glm::length(getDiagonal()) * 0.5f;
    }
};

#endif //PROMETHEUS_BVH_NODE_H