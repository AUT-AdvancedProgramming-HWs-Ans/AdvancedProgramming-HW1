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

    if (n <= 0 || m <= 0) {
        // Throw an exception if n and m are not greater than 0.
        throw std::logic_error("n and m must be greater than 0.");
    }

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

    if (n <= 0 || m <= 0) {
        // Throw an exception if n and m are not greater than 0.
        throw std::logic_error("n and m must be greater than 0.");
    }

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

    if (n <= 0 || m <= 0) {
        // Throw an exception if n and m are not greater than 0.
        throw std::logic_error("n and m must be greater than 0.");
    }

    if (min >= max) {
        // Throw an exception if min is greater than or equal to max.
        throw std::logic_error("min cannot be greater than max");
    }

    // Creating a random number engine
    std::random_device rd;
    // Creating a random number generator
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

// Defining multiply by scalar number function
Matrix algebra::multiply(const Matrix& matrix, double c)
{
    /**
     * @brief This function returns a matrix multiplied by a scalar number.
     * @param matrix The matrix to multiply.
     * @param c The scalar number.
     * @return A matrix multiplied by a scalar number.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    // Creating a matrix equal to the input matrix
    Matrix multiplied_matrix { matrix };

    // Multiplying the matrix by the scalar number
    for (size_t i {}; i < multiplied_matrix.size(); i++) {

        // Multiplying each row by the scalar number
        for (size_t j {}; j < multiplied_matrix[i].size(); j++) {
            multiplied_matrix[i][j] = multiplied_matrix[i][j] * c;
        }
    }

    // Returning the multiplied matrix
    return multiplied_matrix;

} // end of multiply by scalar number

// Defining multiply by matrix function
Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2)
{
    /**
     * @brief This function returns a matrix multiplied by a matrix.
     * @param matrix1 The first matrix.
     * @param matrix2 The second matrix.
     * @return A matrix multiplied by a matrix.
     */

    // Returing an empty matrix if the input matrices are empty
    if (matrix1.empty() && matrix2.empty()) {

        return Matrix {};
    }

    if (matrix1.size() <= 0 || matrix2.size() <= 0 || matrix1[0].size() <= 0 || matrix2[0].size() <= 0) {
        // Throw an exception if the number of rows or columns of matrix1
        throw std::logic_error("The number of rows and columns of the matrices must be greater than zero");
    }

    if (matrix1[0].size() != matrix2.size()) {
        // Throw an exception if the number of columns of matrix1
        // is not equal to the number of rows of matrix2.
        throw std::logic_error("matrices with wrong dimensions cannot be multiplied");
    }

    // Creating a zeros matrix with the dimenson of the multiplication matrix
    Matrix multiplied_matrix { algebra::zeros(matrix1.size(), matrix2[0].size()) };

    // Multiplying the matrices
    for (size_t i {}; i < multiplied_matrix.size(); i++) {

        // Multiplying each row of the first matrix by the  corresponding column of the second matrix
        for (size_t j {}; j < multiplied_matrix[i].size(); j++) {

            // Multiplying each element by the corresponding element in the second matrix
            for (size_t k {}; k < matrix1[0].size(); k++) {
                multiplied_matrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    // Returning the multiplied matrix
    return multiplied_matrix;

} // end of multiply by matrix

// Defining sum by scalar number function
Matrix algebra::sum(const Matrix& matrix, double c)
{
    /**
     * @brief This function returns a matrix summed by a scalar number.
     * @param matrix The matrix to sum.
     * @param c The scalar number.
     * @return A matrix summed by a scalar number.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    // Creating a matrix equal to the input matrix
    Matrix summed_matrix { matrix };

    // Summing the matrix by the scalar number
    for (size_t i {}; i < summed_matrix.size(); i++) {

        // Summing each row by the scalar number
        for (size_t j {}; j < summed_matrix[i].size(); j++) {
            summed_matrix[i][j] = summed_matrix[i][j] + c;
        }
    }

    // Returning the summed matrix
    return summed_matrix;

} // end of sum by scalar number

// Defining sum by matrix function
Matrix algebra::sum(const Matrix& matrix1, const Matrix& matrix2)
{
    /**
     * @brief This function returns a matrix summed by a matrix.
     * @param matrix1 The first matrix.
     * @param matrix2 The second matrix.
     * @return A matrix summed by a matrix.
     */

    // Returinig an empty matrix if the input matrices are empty
    if (matrix1.empty() && matrix2.empty()) {
        return Matrix {};
    }

    if (matrix1.size() <= 0 || matrix2.size() <= 0 || matrix1[0].size() <= 0 || matrix2[0].size() <= 0) {
        // Throw an exception if the number of rows or columns of matrix1
        throw std::logic_error("The number of rows and columns of the matrices must be greater than zero");
    }

    if ((matrix1.size() != matrix2.size()) || (matrix1[0].size() != matrix2[0].size())) {
        // Throw an exception if the number of rows and columns of matrix1
        // are not equal to the number of rows and columns of matrix2.
        throw std::logic_error("matrices with wrong dimensions cannot be summed");
    }

    // Creating a zeros matrix equal to the dimensions of the input matrices
    Matrix summed_matrix { algebra::zeros(matrix1.size(), matrix1[0].size()) };

    // Summing the matrices
    for (size_t i {}; i < summed_matrix.size(); i++) {

        // Summing each row of the first matrix by the corresponding row of the second matrix
        for (size_t j {}; j < summed_matrix[i].size(); j++) {
            summed_matrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }

    // Returning the summed matrix
    return summed_matrix;

} // end of sum by matrix

// Defining transpose function
Matrix algebra::transpose(const Matrix& matrix)
{
    /**
     * @brief This function returns the transpose of a matrix.
     * @param matrix The matrix to transpose.
     * @return The transpose of a matrix.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    // Creating a matrix with the dimensions of transpose matrix
    Matrix transposed_matrix { algebra::zeros(matrix[0].size(), matrix.size()) };

    // Transposing the matrix
    for (size_t i {}; i < transposed_matrix.size(); i++) {

        // Transposing each row
        for (size_t j {}; j < transposed_matrix[i].size(); j++) {
            transposed_matrix[i][j] = matrix[j][i];
        }
    }

    // Returning the transposed matrix
    return transposed_matrix;

} // end of transpose

// Defining minor function
Matrix algebra::minor(const Matrix& matrix, size_t n, size_t m)
{
    /**
     * @brief This function returns the minor of a matrix.
     * @param matrix The matrix to get the minor.
     * @param n The row of the minor.
     * @param m The column of the minor.
     * @return The minor of a matrix.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    if (n >= matrix.size() || m >= matrix[0].size()) {
        // Throw an exception if the row or column of the minor is greater
        // than the number of rows or columns of the matrix
        throw std::logic_error("n and m should be less than the number of rows and columns of the matrix");
    }

    // Creating a matrix with the dimensions of the minor
    Matrix minor_matrix { algebra::zeros(matrix.size() - 1, matrix[0].size() - 1) };

    // Getting the minor
    for (size_t i {}; i < minor_matrix.size(); i++) {

        // Getting the minor of each row
        for (size_t j {}; j < minor_matrix[i].size(); j++) {

            // Getting the minor of each element
            minor_matrix[i][j] = matrix[i < n ? i : i + 1][j < m ? j : j + 1];
        }
    }

    // Returning the minor matrix
    return minor_matrix;

} // end of minor

// Defining determinant function
double algebra::determinant(const Matrix& matrix)
{
    /**
     * @brief This function returns the determinant of a matrix.
     * @param matrix The matrix to get the determinant.
     * @return The determinant of a matrix.
     */

    // Returinig 1 if the input matrix is empty
    if (matrix.empty()) {
        return 1;
    }

    if (matrix.size() != matrix[0].size()) {
        // Throw an exception if the number of rows and columns of matrix
        // are not equal.
        throw std::logic_error("non-square matrices have no determinant");
    }

    // Creating a double equal to the determinant
    if (matrix.size() == 1) {
        return matrix[0][0];
    }

    // Getting the determinant
    double det {};

    // Getting the determinant of first row
    for (size_t i {}; i < matrix.size(); i++) {

        // Getting the determinant of each element
        det += pow(-1, i) * matrix[0][i] * algebra::determinant(algebra::minor(matrix, 0, i));
    }

    // Returning the determinant
    return det;

} // end of determinant

// Defining inverse function
Matrix algebra::inverse(const Matrix& matrix)
{
    /**
     * @brief This function returns the inverse of a matrix.
     * @param matrix The matrix to get the inverse.
     * @return The inverse of a matrix.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    if (matrix.size() != matrix[0].size()) {
        // Throw an exception if the number of rows and columns of matrix
        // are not equal.
        throw std::logic_error("non-square matrices have no inverse");
    }

    // Defining determinant of the matrix
    double det { algebra::determinant(matrix) };

    if (det == 0) {
        // Throw an exception if the determinant of matrix is 0.
        throw std::logic_error("singular matrices have no inverse");
    }

    // Creating a matrix with the dimensions of the inverse
    Matrix inversed_matrix { algebra::zeros(matrix.size(), matrix[0].size()) };

    // Getting the inverse
    for (size_t i {}; i < inversed_matrix.size(); i++) {

        // Getting the inverse of each row
        for (size_t j {}; j < inversed_matrix[i].size(); j++) {

            // Calculating cofactor matrix
            inversed_matrix[j][i] = pow(-1, i + j) * algebra::determinant(algebra::minor(matrix, i, j));
        }
    }

    // Returning the inverse matrix
    return algebra::multiply(inversed_matrix, (1 / det));

} // end of inverse

// Defining concatenate function
Matrix algebra::concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis)
{
    /**
     * @brief This function returns the concatenation of two matrices.
     * @param matrix1 The first matrix to concatenate.
     * @param matrix2 The second matrix to concatenate.
     * @param axis The axis to concatenate.
     * @return The concatenation of two matrices.
     */

    if (axis == 0) {

        // Returinig an empty matrix if the input matrix is empty
        if (matrix1.empty() && matrix2.empty()) {
            return Matrix {};
        }

        if (matrix1[0].size() != matrix2[0].size()) {

            throw ::std::logic_error("matrices with wrong dimensions cannot be concatenated");
        }

        // Creating concatenated matrix
        Matrix concatenated_matrix { matrix1 };

        // Appending the second matrix to the first
        for (size_t i {}; i < matrix2.size(); i++) {
            concatenated_matrix.push_back(matrix2[i]);
        }

        // Returning the concatenated matrix
        return concatenated_matrix;

    } else if (axis == 1) {
        // Returinig an empty matrix if the input matrix is empty
        if (matrix1.empty() && matrix2.empty()) {
            return Matrix {};
        }

        if (matrix1.size() != matrix2.size()) {

            throw ::std::logic_error("matrices with wrong dimensions cannot be concatenated");
        }

        // Creating concatenated matrix
        Matrix concatenated_matrix { algebra::transpose(matrix1) };

        // Appending the second matrix to the first
        for (size_t i {}; i < matrix2[0].size(); i++) {
            concatenated_matrix.push_back(algebra::transpose(matrix2)[i]);
        }

        // Returning the concatenated matrix
        return algebra::transpose(concatenated_matrix);

    } else {
        // Throw an exception if the axis is not 0 or 1.
        throw std::logic_error("axis must be 0 or 1");
    }

} // end of concatenate

// Defining elementary row operation swap function
Matrix algebra::ero_swap(const Matrix& matrix, size_t r1, size_t r2)
{
    /**
     * @brief This function returns the matrix after swapping two rows.
     * @param matrix The matrix to swap rows.
     * @param r1 The first row to swap.
     * @param r2 The second row to swap.
     * @return The matrix after swapping two rows.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    if (r1 >= matrix.size() || r2 >= matrix.size() || r1 < 0 || r2 < 0) {
        // Throw an exception if the row is greater than the number of rows.
        throw std::logic_error("r1 or r2 inputs are out of range");
    }

    // Creating a matrix with the dimensions of the input matrix
    Matrix swaped_matrix { matrix };

    // Swapping the rows
    swaped_matrix[r2].swap(swaped_matrix[r1]);

    // Returning the swaped matrix
    return swaped_matrix;

} // end of ero_swap

// Defining elementary row operation multiply function
Matrix algebra::ero_multiply(const Matrix& matrix, size_t r, double c)
{
    /**
     * @brief This function returns the matrix after multiplying a row by a constant.
     * @param matrix The matrix to multiply a row by a constant.
     * @param r The row to multiply.
     * @param c The constant to multiply by.
     * @return The matrix after multiplying a row by a constant.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    if (r >= matrix.size()) {
        // Throw an exception if the row is greater than the number of rows.
        throw std::logic_error("r input is out of range");
    }

    // Creating a matrix with the dimensions of the input matrix
    Matrix multiplied_matrix { matrix };

    // Multiplying the row by the constant
    for (size_t i {}; i < multiplied_matrix[r].size(); i++) {
        multiplied_matrix[r][i] *= c;
    }

    // Returning the multiplied matrix
    return multiplied_matrix;

} // end of ero_multiply

// Defining elementary row operation sum function
Matrix algebra::ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
{
    /**
     * @brief This function returns the matrix after summing two rows.
     * @param matrix The matrix to sum two rows.
     * @param r1 The first row to sum.
     * @param c The constant to multiply first row by.
     * @param r2 The second row to sum.
     * @return The matrix after summing two rows.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    if (r1 >= matrix.size() || r2 >= matrix.size() || r1 < 0 || r2 < 0) {
        // Throw an exception if the row is greater than the number of rows.
        throw std::logic_error("r1 or r2 inputs are out of range");
    }

    // Creating a matrix with the dimensions of the input matrix
    Matrix summed_matrix { matrix };

    // Summing the rows
    for (size_t i {}; i < summed_matrix[r1].size(); i++) {
        summed_matrix[r2][i] += c * summed_matrix[r1][i];
    }

    // Returning the summed matrix
    return summed_matrix;

} // end of ero_sum

// Defining of upper triangular function
Matrix algebra::upper_triangular(const Matrix& matrix)
{
    /**
     * @brief This function returns the upper triangular matrix of a matrix.
     * @param matrix The matrix to find the upper triangular matrix of.
     * @return The upper triangular matrix of a matrix.
     */

    // Returinig an empty matrix if the input matrix is empty
    if (matrix.empty()) {
        return Matrix {};
    }

    if (matrix.size() != matrix[0].size()) {
        // Throw an exception if the matrix is not square.
        throw std::logic_error("non-square matrices have no upper triangular form");
    }

    if (matrix.size() == 1) {
        // Return the upper triangular matrix if the matrix is 1x1.
        return matrix;
    }

    // Creating a matrix with the dimensions of the input matrix
    Matrix upper_triangular_matrix { matrix };

    // Changing the place of zero elements in the diagonal
    for (size_t i {}; i < upper_triangular_matrix.size(); i++) {

        if (upper_triangular_matrix[i][i] == 0) {
            // Iterating through the rows to find a non-zero element in the diagonal

            // Iterating between the next rows to find a non-zero element
            for (size_t j { i + 1 }; j < upper_triangular_matrix.size(); j++) {
                if ((upper_triangular_matrix[j][i] != 0) && (upper_triangular_matrix[i][j] != 0)) {

                    // Swapping the rows if a non-zero element is found
                    upper_triangular_matrix = algebra::ero_swap(upper_triangular_matrix, i, j);
                    break;
                }
            }
            if (upper_triangular_matrix[i][i] == 0) {
                // Iterating between first rows to the current row to find a non-zero element
                for (size_t j {}; j < i + 1; j++) {
                    if ((upper_triangular_matrix[j][i] != 0) && (upper_triangular_matrix[i][j] != 0)) {

                        // Swapping the rows if a non-zero element is found
                        upper_triangular_matrix = algebra::ero_swap(upper_triangular_matrix, i, j);
                        break;
                    }
                }
            }
        }
    }

    // Looping through the rows
    for (size_t j {}; j < upper_triangular_matrix.size(); j++) {

        // Looping through the columns
        for (size_t i { j + 1 }; i < matrix.size(); i++) {

            // Checking if the element is not zero
            if (upper_triangular_matrix[i][j] != 0) {

                // Summing the rows to get the upper triangular matrix
                upper_triangular_matrix = algebra::ero_sum(upper_triangular_matrix,
                    j,
                    (-upper_triangular_matrix[i][j] / upper_triangular_matrix[j][j]),
                    i);
            }
        }
    }

    // Returning the upper triangular matrix
    return upper_triangular_matrix;
}