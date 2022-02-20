
#include "hw1.h"
#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char** argv)
{
    if (false) // make false to run unit-tests
    {
        // debug section

        // test zeros
        // Matrix m = Algebra::zeros(3, 4);

        // test random

        // Matrix m = algebra::zeros(1.1,3);

        // Matrix m2 =

        //  Matrix matrix(-1, std::vector<double>(-2, 1));

        // algebra::show(matrix);

        // algebra::show(m);

        // Matrix matrix{algebra::multiply(Matrix{}, Matrix{})};

        // Matrix matrix1{};
        // int a{matrix1.size()};
        // std::cout<<a<<std::endl;

        // algebra::show(matrix1);

        // Matrix matrix1 {};
        // Matrix matrix2 = algebra::ones(3, 4);
        // double a { matrix1.size() };
        // std::cout<<matrix1.size()<<std::endl;
        // Matrix matrix3 = algebra::multiply(matrix1, matrix2);

        // algebra::show(matrix3);

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