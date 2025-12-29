//
// Created by ivans on 29/12/2025.
//

#ifndef PROMETHEUS_MIS_VALIDATION_H
#define PROMETHEUS_MIS_VALIDATION_H
#include <iostream>
#include <vector>
#include <cmath>
#include <glm/glm.hpp>
#include "mis/mis_weights.h"
#include "mis/bsdf_sampler.h"

namespace MISValidation {
    /// Test case result
    struct TestResult {
        std::string testName;
        bool passed;
        std::string message;
        float expectedValue;
        float actualValue;
        float relativeError;
    };

    /// Validate MIS weight properties
    std::vector<TestResult> validateMISWeights() {
        std::vector<TestResult> results;

        // test 1: MIS weights should sum to 1
        {
            TestResult test;
            test.testName = "MIS weights sum to 1";

            float pdf1 = 0.3f;
            float pdf2 = 0.7f;

            float weight1 = MISWeightCalculator::balanceHeuristic(pdf1, pdf2);
            float weight2 = MISWeightCalculator::balanceHeuristic(pdf2, pdf1);

            float sum = weight1 + weight2;
            test.expectedValue = 1.0f;
            test.actualValue = sum;
            test.relativeError = std::abs(sum - 1.0f);
            test.passed = (test.relativeError < 1e-5f);
            test.message = test.passed ? "OK" : "FAILED: Weights don't sum to 1";

            results.push_back(test);
        }

        // test 2: MIS weight should be 1 when other PDF is 0
        {
            TestResult test;
            test.testName = "MIS weight = 1 when other PDF = 0";

            float weight = MISWeightCalculator::calculateWeight(0.5f, 0.0f);

            test.expectedValue = 1.0f;
            test.actualValue = weight;
            test.relativeError = std::abs(weight - 1.0f);
            test.passed = (test.relativeError < 1e-5f);
            test.message = test.passed ? "OK" : "FAILED: Should return 1 when other PDF is 0";

            results.push_back(test);
        }

        // test 3: MIS weight should be 0 when own PDF is 0
        {
            TestResult test;
            test.testName = "MIS weight = 0 when own PDF = 0";

            float weight = MISWeightCalculator::calculateWeight(0.0f, 0.5f);

            test.expectedValue = 0.0f;
            test.actualValue = weight;
            test.relativeError = std::abs(weight);
            test.passed = (test.relativeError < 1e-5f);
            test.message = test.passed ? "OK" : "FAILED: Should return 0 when own PDF is 0";

            results.push_back(test);
        }

        // test 4: Balance heuristic symmetry
        {
            TestResult test;
            test.testName = "Balance heuristic is symmetric";

            float pdf1 = 0.4f;
            float pdf2 = 0.6f;

            float weight1 = MISWeightCalculator::balanceHeuristic(pdf1, pdf2);
            float weight2 = MISWeightCalculator::balanceHeuristic(pdf2, pdf1);

            // weight1/pdf1 should equal weight2/pdf2 for balance heuristic
            float ratio1 = weight1 / pdf1;
            float ratio2 = weight2 / pdf2;

            test.expectedValue = ratio1;
            test.actualValue = ratio2;
            test.relativeError = std::abs(ratio1 - ratio2) / std::max(ratio1, ratio2);
            test.passed = (test.relativeError < 1e-5f);
            test.message = test.passed ? "OK" : "FAILED: Balance heuristic not symmetric";

            results.push_back(test);
        }

        // test 5: Power heuristic increases with beta
        {
            TestResult test;
            test.testName = "Power heuristic increases with beta";

            float pdfHigh = 0.8f;
            float pdfLow = 0.2f;

            float weightBeta1 = MISWeightCalculator::powerHeuristic(pdfHigh, pdfLow, 1.0f);
            float weightBeta2 = MISWeightCalculator::powerHeuristic(pdfHigh, pdfLow, 2.0f);

            test.expectedValue = weightBeta2;
            test.actualValue = weightBeta1;
            test.passed = (weightBeta2 > weightBeta1);
            test.message = test.passed ? "OK" : "FAILED: Power heuristic should increase with beta for higher PDF";

            results.push_back(test);
        }

        return results;
    }

    /// Validate BSDF sampling properties
    std::vector<TestResult> validateBSDFSampling() {
        std::vector<TestResult> results;

        // test 1: PDF integrates to 1 (Monte Carlo check)
        {
            TestResult test;
            test.testName = "BSDF PDF integrates to ~1";

            const int numSamples = 10000;
            float sum = 0.0f;
            glm::vec3 normal(0, 1, 0);

            for (int i = 0; i < numSamples; ++i) {
                float u1 = static_cast<float>(i) / numSamples;
                float u2 = static_cast<float>(i * 73) / numSamples;  // Simple decorrelation

                auto sample = BSDFSampler::sampleDiffuse(normal, u1, u2);
                sum += 1.0f / sample.pdf;  // Weight by inverse PDF
            }

            float average = sum / numSamples;
            test.expectedValue = 2.0f * glm::pi<float>();  // Hemisphere solid angle
            test.actualValue = average;
            test.relativeError = std::abs(average - test.expectedValue) / test.expectedValue;
            test.passed = (test.relativeError < 0.05f);  // 5% tolerance for MC
            test.message = test.passed ? "OK" : "FAILED: PDF doesn't integrate to hemisphere solid angle";

            results.push_back(test);
        }

        // test 2: Samples respect cosine distribution
        {
            TestResult test;
            test.testName = "BSDF samples follow cosine distribution";

            const int numSamples = 10000;
            float avgCosTheta = 0.0f;
            glm::vec3 normal(0, 1, 0);

            for (int i = 0; i < numSamples; ++i) {
                float u1 = static_cast<float>(i) / numSamples;
                float u2 = static_cast<float>(i * 73) / numSamples;

                auto sample = BSDFSampler::sampleDiffuse(normal, u1, u2);
                float cosTheta = glm::dot(sample.direction, normal);
                avgCosTheta += cosTheta;
            }

            avgCosTheta /= numSamples;
            test.expectedValue = 2.0f / glm::pi<float>();  // Expected value of cosine-weighted hemisphere
            test.actualValue = avgCosTheta;
            test.relativeError = std::abs(avgCosTheta - test.expectedValue) / test.expectedValue;
            test.passed = (test.relativeError < 0.05f);  // 5% tolerance
            test.message = test.passed ? "OK" : "FAILED: Samples don't follow cosine distribution";

            results.push_back(test);
        }

        // test 3: PDF evaluation matches sampling PDF
        {
            TestResult test;
            test.testName = "PDF evaluation consistent with sampling";

            glm::vec3 normal(0, 1, 0);
            float u1 = 0.5f;
            float u2 = 0.3f;

            auto sample = BSDFSampler::sampleDiffuse(normal, u1, u2);
            float evaluatedPdf = BSDFSampler::pdfDiffuse(normal, sample.direction);

            test.expectedValue = sample.pdf;
            test.actualValue = evaluatedPdf;
            test.relativeError = std::abs(evaluatedPdf - sample.pdf) / std::max(sample.pdf, 1e-6f);
            test.passed = (test.relativeError < 1e-4f);
            test.message = test.passed ? "OK" : "FAILED: Evaluated PDF doesn't match sample PDF";

            results.push_back(test);
        }

        return results;
    }

    /// Validate combined MIS + sampling workflow
    std::vector<TestResult> validateCombinedWorkflow() {
        std::vector<TestResult> results;

        // test: Variance reduction with MIS

        TestResult test;
        test.testName = "MIS reduces variance vs single strategy";

        const int numSamples = 1000;
        glm::vec3 normal(0, 1, 0);

        // Simulate simple case: comparing MIS to light-only sampling
        float varianceMIS = 0.0f;
        float varianceLightOnly = 0.0f;
        float meanMIS = 0.0f;
        float meanLightOnly = 0.0f;

        for (int i = 0; i < numSamples; ++i) {
            float u = static_cast<float>(i) / numSamples;

            // Simulate light PDF and BSDF PDF
            float lightPdf = 0.5f + 0.5f * u;  // Varying PDF
            float bsdfPdf = 0.3f + 0.3f * std::cos(u * glm::pi<float>());

            // MIS weight
            float misWeight = MISWeightCalculator::balanceHeuristic(lightPdf, bsdfPdf);

            // Simulated contribution
            float contrib = 1.0f / lightPdf;
            float contribMIS = contrib * misWeight;

            meanMIS += contribMIS;
            meanLightOnly += contrib;
        }

        meanMIS /= numSamples;
        meanLightOnly /= numSamples;

        // Calculate variance
        for (int i = 0; i < numSamples; ++i) {
            float u = static_cast<float>(i) / numSamples;
            float lightPdf = 0.5f + 0.5f * u;
            float bsdfPdf = 0.3f + 0.3f * std::cos(u * glm::pi<float>());
            float misWeight = MISWeightCalculator::balanceHeuristic(lightPdf, bsdfPdf);

            float contrib = 1.0f / lightPdf;
            float contribMIS = contrib * misWeight;

            varianceMIS += (contribMIS - meanMIS) * (contribMIS - meanMIS);
            varianceLightOnly += (contrib - meanLightOnly) * (contrib - meanLightOnly);
        }

        varianceMIS /= numSamples;
        varianceLightOnly /= numSamples;

        test.expectedValue = varianceLightOnly;
        test.actualValue = varianceMIS;
        test.passed = (varianceMIS < varianceLightOnly);
        test.message = test.passed ?
            "OK: MIS variance lower than single strategy" :
            "WARNING: MIS didn't reduce variance (may be acceptable in some cases)";

        results.push_back(test);
        return results;
    }

    /// run all validation tests
    inline void runAllTests() {
        std::cout << "\n=== MIS Implementation Validation ===" << std::endl;
        std::cout << "=====================================" << std::endl;

        auto misTests = MISValidation::validateMISWeights();
        auto bsdfTests = MISValidation::validateBSDFSampling();
        auto workflowTests = MISValidation::validateCombinedWorkflow();

        int totalTests = 0;
        int passedTests = 0;

        auto printResults = [&](const std::vector<TestResult>& tests, const std::string& category) {
            std::cout << "\n" << category << ":" << std::endl;
            for (const auto& test : tests) {
                totalTests++;
                if (test.passed) passedTests++;

                std::cout << "  [" << (test.passed ? "PASS" : "FAIL") << "] "
                          << test.testName << std::endl;
                std::cout << "    Expected: " << test.expectedValue
                          << ", Actual: " << test.actualValue
                          << ", Error: " << (test.relativeError * 100.0f) << "%" << std::endl;
                std::cout << "    " << test.message << std::endl;
            }
        };

        printResults(misTests, "MIS Weight Tests");
        printResults(bsdfTests, "BSDF Sampling Tests");
        printResults(workflowTests, "Combined Workflow Tests");

        std::cout << "\n=====================================" << std::endl;
        std::cout << "Results: " << passedTests << "/" << totalTests << " tests passed" << std::endl;

        if (passedTests == totalTests) {
            std::cout << "✓ All tests passed! MIS implementation looks correct." << std::endl;
        } else {
            std::cout << "✗ Some tests failed. Please review the implementation." << std::endl;
        }
        std::cout << "=====================================" << std::endl;
    }
}
#endif //PROMETHEUS_MIS_VALIDATION_H