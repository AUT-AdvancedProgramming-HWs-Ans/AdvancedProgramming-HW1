/**
 * In this code we are going to implement a Linear Algebra library.
 * @file hw1.c
 * @author Erfan Rasti
 * @version 1.0.0
 */

// Adding libraries
#include "hw1.h"

// Defining zeros function
Matrix algebra::zeros(size_t n, size_t m)
{
    /**
     * @brief This function returns a matrix of zeros.
     * @param n The number of rows.
     * @param m The number of columns.
     * @return A matrix of zeros.
     */

    return Matrix(n, std::vector<double>(m, 0));
} // end of zeros

// Defining ones function
Matrix algebra::ones(size_t n, size_t m)
{
    /**
     * @brief This function returns a matrix of ones.
     * @param n The number of rows.
     * @param m The number of columns.
     * @return A matrix of ones.
     */

    return Matrix(n, std::vector<double>(m, 1));
} // end of ones

// Defining random function
Matrix algebra::random(size_t n, size_t m, double min, double max)
{
    /**
     * @brief This function returns a matrix of random numbers.
     * @param n The number of rows.
     * @param m The number of columns.
     * @param min The minimum value of the random numbers.
     * @param max The maximum value of the random numbers.
     * @return A matrix of random numbers.
     */

    
    if (min >= max) {
        // Throw an error
        throw std::invalid_argument("Caution: min cannot be greater than max");
    }

    // Creating a random number generator
    std::random_device rd;
    // Creating a random number engine
    std::mt19937 gen(rd());
    
    // Creating a uniform distribution
    std::uniform_real_distribution<double> dist(min, max);

    // Creating a matrix of random numbers
    return Matrix(n, std::vector<double>(m, dist(gen)));

} // end of random

// Defining show function
void algebra::show(const Matrix& matrix)
{
    /**
     * @brief This function prints a matrix.
     * @param matrix The matrix to print.
     */

    // Printing the matrix
    for (size_t i {}; i < matrix.size(); i++) {

        // Printing each row
        for (size_t j {}; j < matrix[i].size(); j++) {
            std::cout << std::setw(10) << matrix[i][j];
        }
        // New line
        std::cout << std::endl;
    }
} // end of show

// Defining multiply by constant function
Matrix algebra::multiply(const Matrix& matrix, double c)
{
    /**
     * @brief This function returns a matrix multiplied by a constant.
     * @param matrix The matrix to multiply.
     * @param c The constant.
     * @return A matrix multiplied by a constant.
     */

    // Creating a matrix equal to the input matrix
    Matrix temp = matrix;

    // Multiplying the matrix by the constant
    for (size_t i {}; i < temp.size(); i++) {

        // Multiplying each row by the constant
        for (size_t j {}; j < temp[i].size(); j++) {
            temp[i][j] = temp[i][j] * c;
        }
    }

    // Returning the temp
    return temp;

} // end of multiply

// Defining multiply by matrix function
// Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2)
// {
//     /**
//      * @brief This function returns a matrix multiplied by a matrix.
//      * @param matrix1 The first matrix.
//      * @param matrix2 The second matrix.
//      * @return A matrix multiplied by a matrix.
//      */

//     // Creating a matrix with the dimenson of the multiplication
//     Matrix temp = algebra::zeros(matrix1.size(), matrix2[0].size());

//     // Calculating transpose of matrix2
//     Matrix matrix2_transpose = algebra::transpose(matrix2);


//     // Returning the temp
//     return temp;

// } // end of multiply

