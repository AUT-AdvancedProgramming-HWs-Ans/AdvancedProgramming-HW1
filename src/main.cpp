/**
 * This code tests the functions in hw1.cpp file.
 * @file main.cpp
 * @author Erfan Rasti
 * @version 1.0.1
 */

#include "hw1.h"
#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char** argv)
{
    if (false) // make false to run unit-tests
    {
        // debug section

        // // test zeros
        // Matrix m1 = algebra::zeros(3, 0);

        // // test show
        //  Matrix matrix(-1, std::vector<double>(-2, 1));
        // algebra::show(matrix);

        // // test multiply
        // Matrix matrix1{algebra::multiply(Matrix{}, Matrix{})};
        // std::cout<<matrix1.size()<<std::endl;
        // algebra::show(matrix1);

        // Matrix matrix2 {};
        // Matrix matrix3 = algebra::ones(3, 4);
        // std::cout<<matrix2.size()<<std::endl;
        // Matrix matrix4 = algebra::multiply(matrix2, matrix3);

        // // test upper triangular
        // Matrix matrix_test_upper_triangular{{1, 1, 1}, {0, 0, 1}, {1, 0, 1}};
        // Matrix matrix_final_upper_triangular = algebra::upper_triangular(matrix_test_upper_triangular);
        // algebra::show(matrix_final_upper_triangular);

    } else {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret { RUN_ALL_TESTS() };
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}