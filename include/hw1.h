/**
 * This code is list of headers that are used in hw1.cpp file.
 * @file hw1.h
 * @author Erfan Rasti
 * @version 1.0.2
 */

#ifndef AP_HW1_H
#define AP_HW1_H

// Defining libraries
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using Matrix = std::vector<std::vector<double>>;

namespace algebra {
// Prototypes

// Header of zeros function
Matrix zeros(size_t n, size_t m);

// Header of ones function
Matrix ones(size_t n, size_t m);

// Header of random function
Matrix random(size_t n, size_t m, double min, double max);

// Header of show function
void show(const Matrix& matrix);

// Header of multiply by scalar number function
Matrix multiply(const Matrix& matrix, double c);

// Header of multiply by matrix function
Matrix multiply(const Matrix& matrix1, const Matrix& matrix2);

// Hedaer of sum by scalar number function
Matrix sum(const Matrix& matrix, double c);

// Header of sum by matrix function
Matrix sum(const Matrix& matrix1, const Matrix& matrix2);

// Header of transpose function
Matrix transpose(const Matrix& matrix);

// Header of minor function
Matrix minor(const Matrix& matrix, size_t n, size_t m);

// Header of determinant function
double determinant(const Matrix& matrix);

// Header of inverse function
Matrix inverse(const Matrix& matrix);

// Header of concatenate function
Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis = 0);

// Header of elementary row operation swap function
Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2);

// Header of elementary row operation multiply function
Matrix ero_multiply(const Matrix& matrix, size_t r, double c);

// Header of elementary row operation sum function
Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2);

// Header of upper triangular function
Matrix upper_triangular(const Matrix& matrix);

} // end of namespace algebra

#endif // AP_HW1_H
