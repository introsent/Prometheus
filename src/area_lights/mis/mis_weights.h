//
// Created by ivans on 28/12/2025.
//

#ifndef PROMETHEUS_MIS_WEIGHTS_H
#define PROMETHEUS_MIS_WEIGHTS_H

#include <cmath>

/// MIS heuristics
// different strategies for combining sampling techniques
enum class MISHeuristic {
    Balance,    // w_i = pdf_i / sum(pdf_j)
    Power,      // w_i = pdf_i^2 / sum(pdf_j^2)  [default, variance optimal]
    Maximum     // w_i = 1 if argmax(pdf_j), else 0
};

/// MIS weight calculator
// implements Veach & Guibas MIS formulas
class MISWeightCalculator {
public:
    // calculate MIS weight for a sample
    // thisPdf: pdf for the technique that generated this sample
    // otherPdf: pdf for the other technique
    static float calculateWeight(
        float thisPdf,
        float otherPdf,
        MISHeuristic heuristic = MISHeuristic::Balance) {

        if (thisPdf <= 0.0f) return 0.0f;
        if (otherPdf <= 0.0f) return 1.0f;

        switch (heuristic) {
            case MISHeuristic::Balance:
                return balanceHeuristic(thisPdf, otherPdf);
            case MISHeuristic::Power:
                return powerHeuristic(thisPdf, otherPdf);
            case MISHeuristic::Maximum:
                return thisPdf > otherPdf ? 1.0f : 0.0f;
            default:
                return balanceHeuristic(thisPdf, otherPdf);
        }
    }

    // balance heuristic: w_i = n_i * pdf_i / sum(n_j * pdf_j)
    // for one sample each: w_i = pdf_i / (pdf_i + pdf_j)
    static float balanceHeuristic(float thisPdf, float otherPdf) {
        return thisPdf / (thisPdf + otherPdf);
    }

    // power heuristic: w_i = (pdf_i)^beta / sum(pdf_j)^beta
    // beta=2 is variance optimal (Veach thesis)
    static float powerHeuristic(float thisPdf, float otherPdf, float beta = 2.0f) {
        const float thisWeight = std::pow(thisPdf, beta);
        const float otherWeight = std::pow(otherPdf, beta);
        return thisWeight / (thisWeight + otherWeight);
    }
};

#endif //PROMETHEUS_MIS_WEIGHTS_H