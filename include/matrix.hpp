// matrix.hpp
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>

template <typename T>
class Matrix {
private:
    std::vector<T> data;  // 1D vector to store matrix elements
    size_t rows, cols;    // number of rows and columns

public:
    Matrix(size_t r, size_t c, T initial = T()) : rows(r), cols(c), data(r * c, initial) {}

    T& operator()(size_t i, size_t j) {
        return data[i * cols + j];
    }

    const T& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }

    void print() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                double printedVariable = (*this)(i, j);
                std::cout << printedVariable << " ";
            }
            std::cout << std::endl;
        }
    }
};

#endif  // MATRIX_HPP
