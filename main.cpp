#include "include\matrix.hpp"
#include "include\simplex.hpp"
#include <iostream>
#include <vector>

// compile command g++ -std=c++17 -Iinclude main.cpp src/matrix.cpp src/simplex.cpp -o simplex_solver
// or              g++ --% -std=c++17 -Iinclude main.cpp src/matrix.cpp src/simplex.cpp -o simplex_solver

int main() {
    /*
    Matrix<double> A(2, 2);
    std::vector<double> b = {12, 16};
    std::vector<double> c = {8, 5};

    A(0, 0) = 1;
    A(0, 1) = 1;
    A(1, 0) = 2;
    A(1, 1) = 1;
    std::vector<double> b = {12, 16};
    std::vector<double> c = {8, 5};
    */
    Matrix<double> A(3, 3);
    std::vector<double> b = {2, 5, 6};
    std::vector<double> c = {3, 1, 3};

    A(0, 0) = 2;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 0) = 1;
    A(1, 1) = 2;
    A(1, 2) = 3;
    A(2, 0) = 2;
    A(2, 1) = 2;
    A(2, 2) = 1;

    std::cout << "Matrix A:" << std::endl;
    A.print();

    bool isMaximization = true;

    Simplex simplex(A, b, c, isMaximization);

    std::cout << "Initial Tableau:" << std::endl;
    simplex.printTableau();

//    simplex.run_simplex_step(isMaximization);
//    std::cout << "Tableau after one Simplex step:" << std::endl;
//    simplex.printTableau();

//    simplex.run_simplex_step(isMaximization);
//    std::cout << "Tableau after two Simplex steps:" << std::endl;
    simplex.run_simplex_algorithm(isMaximization);
    std::cout << "Result:" << std::endl;
    simplex.printTableau();

    return 0;
}
