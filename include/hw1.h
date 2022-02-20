#ifndef AP_HW1_H
#define AP_HW1_H
#endif // AP_HW1_H

// Defining libraries
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using Matrix = std::vector<std::vector<double>>;

namespace algebra {
    // Prototypes
    Matrix zeros(size_t n, size_t m);
    Matrix ones(size_t n, size_t m);
    Matrix random(size_t n, size_t m, double min, double max);
    void show(const Matrix& matrix);
    Matrix multiply(const Matrix& matrix, double c);
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2);

} // end of namespace algebra