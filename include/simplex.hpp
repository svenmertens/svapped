// simplex.hpp
#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include "matrix.hpp"
#include <vector>
#include <limits>

const int HIGHVALUE = 100000;

class Simplex {
private:
    Matrix<double> tableau;
    std::vector<int> basicVariables;
    std::vector<int> nonBasicVariables;

public:
    Simplex(const Matrix<double>& A, const std::vector<double>& b, const std::vector<double>& c, bool isMaximization);
    
    void initializeTableau(const Matrix<double>& A, const std::vector<double>& b, const std::vector<double>& c, bool isMaximization);
    size_t getMinElementIndex(const std::vector<double>& vector) const;
    size_t getPivotColumn(bool isMaximization);
    size_t getPivotRow2(size_t pivotColumnIndex);
    size_t getPivotRow(size_t pivotColumnIndex);
    void standardizePivotRow(size_t pivotRow, size_t pivotColumn);
    void tableauPivotOperations(size_t pivotRow, size_t pivotColumn);
    size_t run_simplex_step(bool isMaximization, size_t pivotColumnIndex);
    void run_simplex_algorithm(bool isMaximization);
    std::vector<double> getSolution();
    double getObjectiveValue();
    void printTableau() const;
};

#endif  // SIMPLEX_HPP
