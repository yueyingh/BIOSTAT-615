#include "admm_utils.h"

// Soft thresholding function.
std::vector<double> soft_threshold(const std::vector<double>& v, double lambda) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] > lambda) {
            result[i] = v[i] - lambda;
        } else if (v[i] < -lambda) {
            result[i] = v[i] + lambda;
        } else {
            result[i] = 0;
        }
    }
    return result;
}

// Implementations of other utility functions go here.
