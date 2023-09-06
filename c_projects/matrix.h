#pragma once

#include <iostream>
#include <complex>
#include <memory>

class matrix
{
  // Friends
  friend std::ostream & operator<<(std::ostream &os, matrix& mat);
protected:
  //std::complex<double> *matrix_data {nullptr};
  std::unique_ptr<std::complex<double>[]> matrix_data {nullptr};
  size_t rows{0};
  size_t columns{0};
public:
  // Constructors and destructor
  matrix()= default;
  // Parameterized  and conpy constructors
  matrix(size_t n, size_t m);
  matrix(const matrix &other);
  // destructor
  ~matrix()= default;

  // Access functions
  int get_rows() const {return rows;} // Return number of rows
  int get_cols() const {return columns;} // Return number of columns
  int index(size_t m, size_t n) const; // Return position in array of element (m,n)
  std::complex<double> & operator() (size_t m, size_t n) const {return matrix_data[index(m,n)];} // Returen element (m,n)

  // operators
  matrix & operator=(const matrix &arr2d);
  matrix operator*(matrix mat) const;

  // tensor product function
  matrix tensor_product(const matrix& mat2) const;
};