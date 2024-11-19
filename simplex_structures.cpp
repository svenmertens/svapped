#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>


int HIGHVALUE = 100000;

template <typename T>
class Matrix {
private:
    std::vector<T> data;  // 1D vector to store matrix elements
    size_t rows, cols;    // number of rows and columns

public:
    // Constructor to initialize matrix with a given size and initial value
    Matrix(size_t r, size_t c, T initial = T()) : rows(r), cols(c), data(r * c, initial) {}

    // Overload () operator to access elements with matrix(i, j)
    T& operator()(size_t i, size_t j) {
        return data[i * cols + j];
    }

    // Const version for read-only access
    const T& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    // Get number of rows
    size_t getRows() const {
        return rows;
    }

    // Get number of columns
    size_t getCols() const {
        return cols;
    }

    // Print the matrix for debugging purposes
    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
};


class Simplex {
private:
    Matrix<double> tableau;              // Tableau matrix
    std::vector<int> basicVariables;     // Basic variables (tracking indices of basic vars)
    std::vector<int> nonBasicVariables;  // Non-basic variables (tracking indices)
    
public:
    /* Constructor that initializes the tableau with A, b, and c 
    A represents the lefthand side of the constraints
    b represents the righthand side of the constraints
    c represents the objective function */
    Simplex(const Matrix<double>& A, const std::vector<double>& b, const std::vector<double>& c)
        : tableau(A.getRows() + 1, A.getCols() + A.getRows() + 1) {  // tableau includes slack variables
        initializeTableau(A, b, c);
    }

    // Function to initialize the tableau
    void initializeTableau(const Matrix<double>& A, const std::vector<double>& b, const std::vector<double>& c) {
        size_t m = A.getRows();  // Number of constraints
        size_t n = A.getCols();  // Number of variables

        // Initialize tableau by copying A, b, and c
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                tableau(i, j) = A(i, j);  // Copy A matrix
            }
            tableau(i, n + i) = 1.0;       // Add slack variables (identity matrix)
            tableau(i, n + m) = b[i];      // Right-hand side (b)
        }

        // Add objective function in the last row
        for (size_t j = 0; j < n; ++j) {
            tableau(m, j) = -c[j];         // Minimization requires negating c
        }
    }

    size_t getMinElementIndex(std::vector<double> vector) const {
        auto minElementIt = std::min_element(vector.begin(), vector.end());
        int index = std::distance(vector.begin(), minElementIt);
        return index;
    }

    size_t getPivotColumn() {
        size_t lastRow = tableau.getRows() - 1;
        for (size_t col = 0; col < tableau.getCols(); col++) {
            double objectiveFunctionVariable = tableau(lastRow, col);
            if (objectiveFunctionVariable < 0) {
                size_t pivotColumnIndex = col;
                return pivotColumnIndex;
            }
        }
        return std::numeric_limits<size_t>::max();
    }

    size_t getPivotRow(size_t pivotColumnIndex) {
        size_t rowCount = tableau.getRows();
        size_t colCount = tableau.getCols();
        std::vector<double> RHSColumnRatio(rowCount);
        for (size_t row = 0; row < rowCount; row++) {
            double columnValue = tableau(row, pivotColumnIndex);
            double rightHandSide = tableau(row, colCount);
            if (columnValue != 0) {
                RHSColumnRatio[row] = rightHandSide / columnValue;
            } else {
                RHSColumnRatio[row] = HIGHVALUE;
            }
        }
        size_t pivotRow = getMinElementIndex(RHSColumnRatio);
        return pivotRow;
    }

    void standardizePivotRow(size_t pivotRow, size_t pivotColumn) {
        double pivotElement = tableau(pivotRow, pivotColumn);
        if (pivotElement != 1) {
            for (size_t col = 0; col < tableau.getCols(); col++) {
                tableau(pivotRow, col) = tableau(pivotRow, col) / pivotElement;
            }
        }
    }

    void tableauPivotOperations(size_t pivotRow, size_t pivotColumn) {
        for (size_t row = 0; row < tableau.getRows(); row++) {
            double pivotColumnElement = tableau(row, pivotColumn);
            if (pivotColumnElement != 0 && row != pivotRow) {
                double additionFactor = tableau(row, pivotColumn);
                for (size_t col = 0; col < tableau.getCols(); col++) {
                    tableau(row, col) = tableau(row, col) - additionFactor * tableau(pivotRow, col);
                }
            }
        }
    }

    // Main Simplex algorithm function
    void run_simplex_step() {
        size_t pivotColumnIndex = getPivotColumn();
        if (pivotColumnIndex == std::numeric_limits<size_t>::max()) {
            return;
        }
        size_t pivotRowIndex = getPivotRow(pivotColumnIndex);
        standardizePivotRow(pivotRowIndex, pivotColumnIndex);
        tableauPivotOperations(pivotRowIndex, pivotColumnIndex);
    }

    // Extract the solution from the final tableau
    std::vector<double> getSolution() {
        size_t m = tableau.getRows();
        size_t n = tableau.getCols() - 1; // Last column is the RHS (b column)

        std::vector<double> solution(n - 1, 0.0); // Solution vector excluding slack variables

        // Find the basic variables in the tableau
        for (size_t i = 0; i < m - 1; ++i) {
            for (size_t j = 0; j < n - 1; ++j) {
                if (tableau(i, j) == 1 && tableau(m - 1, j) == 0) {
                    solution[j] = tableau(i, n - 1);  // If it's a basic variable, assign the RHS value
                    break;
                }
            }
        }

        return solution;
    }

    // Get the value of the objective function
    double getObjectiveValue() {
        size_t m = tableau.getRows();
        size_t n = tableau.getCols();

        // Return the objective function value (bottom-right corner of the tableau)
        return tableau(m - 1, n - 1);
    }

    // Optional: Print the tableau for debugging purposes
    void printTableau() const {
        tableau.print();
    }
};



int main() {
    Matrix<double> A(2, 2);
    std::vector<double> b = {12, 16};

    // Define vector c for objective function coefficients (1x2)
    std::vector<double> c = {8, 5};

    A(0, 0) = 1;
    A(0, 1) = 1;
    A(1, 0) = 2;
    A(1, 1) = 1;

    // Print the matrix
    A.print();

    std::cout << "b: ";
    for (double val : b) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "c: ";
    for (double val : c) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    int matrixRows = A.getRows();
    std::cout << "Matrix rows: " << matrixRows << std::endl; 

    Simplex tableau(A, b, c);

    std::cout << "Print Tableau: " << std::endl;
    tableau.printTableau();

    size_t col = tableau.getPivotColumn();
    std::cout << "Pivot column: " << col << std::endl;

    tableau.run_simplex_step();

    std::cout << "Print Tableau after Simplex: " << std::endl;
    tableau.printTableau();

    return 0;
}
