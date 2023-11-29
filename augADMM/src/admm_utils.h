#ifndef ADMM_UTILS_H
#define ADMM_UTILS_H

#include <vector>

// Soft thresholding function for the ADMM z-update.
std::vector<double> soft_threshold(const std::vector<double>& v, double lambda);

// Add any other utility function declarations here.

#endif // ADMM_UTILS_H
