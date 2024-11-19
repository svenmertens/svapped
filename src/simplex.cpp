/**
 * @file simplex.cpp
 * @brief Contains the Simplex class with related functions and the Simplex method.
 * @details This file uses the matrix class to create tableaus for the Simplex method and allows returning a finished tableau
 * after running Simplex steps.
 */
#include "..\include\simplex.hpp"
#include <algorithm>
#include <iostream>

// ===== Helper Function Declarations =====

/**
 * @brief Initializes the tableau for the Simplex algorithm.
 * @param A The matrix of constraint coefficients.
 * @param b The right-hand side vector of the constraints.
 * @param c The objective function coefficients.
 * @param isMaximization Indicates whether the problem is a maximization problem.
 */
void Simplex::initializeTableau(const Matrix<double>& A, const std::vector<double>& b, const std::vector<double>& c, bool isMaximization) {
    size_t m = A.getRows(); // Get the amount of defined constraints
    size_t n = A.getCols(); // Get the amount of defined variables

    // Populate the template with the variables by looping through the rows and then the columns
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            tableau(i, j) = A(i, j); // Adding the variables to the template
        }
        tableau(i, n + i) = 1.0; // Adding the slack variables to the template, adding it one column further to the right than in the previous column
        tableau(i, n + m) = b[i]; // Adding the right-hand side to the template at the rightmost column
    }

    for (size_t j = 0; j < n; ++j) {
        tableau(m, j) = isMaximization ? c[j] : -c[j]; /* Add the objective function's values to the last row with the sign depending on whether the problem 
        is a minimization or maximization problem */
    }
}

// ===== Constructor =====

/**
 * @brief Constructs a Simplex object and initializes the tableau.
 * @param A The matrix of constraint coefficients.
 * @param b The right-hand side vector of the constraints.
 * @param c The objective function coefficients.
 * @param isMaximization Indicates whether the problem is a maximization problem.
 */
Simplex::Simplex(const Matrix<double>& A, const std::vector<double>& b, const std::vector<double>& c, bool isMaximization)
    : tableau(A.getRows() + 1, A.getCols() + A.getRows() + 1) {
    initializeTableau(A, b, c, isMaximization);
}

// ===== Helper Function Definitions =====

size_t Simplex::getMinElementIndex(const std::vector<double>& vector) const {
    return std::distance(vector.begin(), std::min_element(vector.begin(), vector.end()));
}

size_t Simplex::getPivotColumn(bool isMaximization) {
    size_t lastRow = tableau.getRows() - 1;
    for (size_t col = 0; col < tableau.getCols(); col++) {
        if (isMaximization == true) {
            if (tableau(lastRow, col) > 0) {
                return col;
            }
        }
        if (isMaximization == false) {
            if (tableau(lastRow, col) < 0) {
                return col;
            }
        }
    }
    return std::numeric_limits<size_t>::max();
}

size_t Simplex::getPivotRow2(size_t pivotColumnIndex) {
    size_t rowCount = tableau.getRows();
    size_t colCount = tableau.getCols();
    std::vector<double> RHSColumnRatio(rowCount - 1);
    for (size_t row = 0; row < rowCount - 1; row++) {
        double columnValue = tableau(row, pivotColumnIndex);
        double rightHandSide = tableau(row, colCount - 1);
        RHSColumnRatio[row] = columnValue != 0 ? rightHandSide / columnValue : HIGHVALUE;
//        std::cout << "Righthandside: " << rightHandSide << std::endl;
//        std::cout << "RowCount: " << rowCount << "  " << "ColCount: " << colCount << std::endl;
//        std::cout << "Element from row " << row << " and column 1: " << tableau(row, 1) << std::endl;
    }
    return getMinElementIndex(RHSColumnRatio);
}

size_t Simplex::getPivotRow(size_t pivotColumnIndex) {
    size_t rowCount = tableau.getRows();
    size_t colCount = tableau.getCols();
    size_t minIndex = std::numeric_limits<size_t>::max();
    double minValue = std::numeric_limits<double>::max();

    for (size_t row = 0; row < rowCount - 1; row++) {
        double columnValue = tableau(row, pivotColumnIndex);
        double rightHandSide = tableau(row, colCount - 1);
        if (columnValue > 0) { // Only consider valid ratios
            double ratio = rightHandSide / columnValue;
            if (ratio < minValue || (ratio == minValue && row < minIndex)) {
                minValue = ratio;
                minIndex = row; // Tie-breaking by smallest row index
            }
        }
    }
    return minIndex;
}

void Simplex::standardizePivotRow(size_t pivotRow, size_t pivotColumn) {
    double pivotElement = tableau(pivotRow, pivotColumn);
    for (size_t col = 0; col < tableau.getCols(); col++) {
        tableau(pivotRow, col) /= pivotElement;
    }
}

void Simplex::tableauPivotOperations(size_t pivotRow, size_t pivotColumn) {
    for (size_t row = 0; row < tableau.getRows(); row++) {
        if (row != pivotRow) {
            double factor = tableau(row, pivotColumn);
            for (size_t col = 0; col < tableau.getCols(); col++) {
                tableau(row, col) -= factor * tableau(pivotRow, col);
            }
        }
    }
}

size_t Simplex::run_simplex_step(bool isMaximization, size_t pivotColumnIndex) {
    size_t pivotRowIndex = getPivotRow(pivotColumnIndex);
//    std::cout << "Pivot row index: " << pivotRowIndex << std::endl << "Pivot col index: " << pivotColumnIndex << std::endl;
    standardizePivotRow(pivotRowIndex, pivotColumnIndex);
    tableauPivotOperations(pivotRowIndex, pivotColumnIndex);
    pivotColumnIndex = getPivotColumn(isMaximization);
    return pivotColumnIndex;
}

void Simplex::run_simplex_algorithm(bool isMaximization) {
    size_t maxIterations = 10;
    size_t iterationCount = 0;
    size_t pivotColumnIndex = getPivotColumn(isMaximization);

    while (pivotColumnIndex != std::numeric_limits<size_t>::max()) {
        if (++iterationCount > maxIterations) {
            std::cout << "Tableau:" << std::endl;
            tableau.print();
            throw std::runtime_error("Simplex algorithm exceeded maximum iterations. Possible cycling detected.");
        }
        pivotColumnIndex = run_simplex_step(isMaximization, pivotColumnIndex);
    }
}

std::vector<double> Simplex::getSolution() {
    size_t m = tableau.getRows();
    size_t n = tableau.getCols() - 1;
    std::vector<double> solution(n - 1, 0.0);

    for (size_t i = 0; i < m - 1; ++i) {
        for (size_t j = 0; j < n - 1; ++j) {
            if (tableau(i, j) == 1 && tableau(m - 1, j) == 0) {
                solution[j] = tableau(i, n - 1);
                break;
            }
        }
    }

    return solution;
}

double Simplex::getObjectiveValue() {
    size_t m = tableau.getRows();
    size_t n = tableau.getCols();
    return tableau(m - 1, n - 1);
}

void Simplex::printTableau() const {
    tableau.print();
}
