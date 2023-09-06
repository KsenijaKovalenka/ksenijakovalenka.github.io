# include "matrix.h"
#include<cmath>

// Parameterized constructor
matrix::matrix(size_t n, size_t m)
{
if(n<1 || m<1)
	{
	std::cerr<<"Error: trying to declare a matrix with one of the dimentions of size < 1"<<std::endl;
	exit(1);
	}
rows = n;
columns = m;
matrix_data = std::make_unique<std::complex<double>[]>(rows * columns);
for(size_t i{}; i<rows*columns; i++) {
	matrix_data[i]= 0;
}
}

//copy constructor  
matrix::matrix(const matrix &other)
: matrix_data(std::make_unique<std::complex<double>[]>(other.rows * other.columns)),
rows(other.rows),
columns(other.columns)
{
	// Copy values from the other matrix
	for (size_t i = 0; i < rows * columns; i++) {
		matrix_data[i] = other.matrix_data[i];
	}
}


// get an array index from the row and column ordering
int matrix::index(size_t m, size_t n) const // Return position in array of element (m,n)
{
	if(m>0 && m<=rows && n>0 && n<=columns) return (n-1)+(m-1)*columns;
	else {std::cerr<<"Error: out of range"<<std::endl; exit(1);}
}

// assignment
matrix & matrix::operator=(const matrix &arr2d)
{
	if(&arr2d == this) return *this; // no self assignment
	// First delete this object's array
	matrix_data.reset(); rows=0; columns=0;
	// Now copy size and declare new array
	rows=arr2d.get_rows();
	columns=arr2d.get_cols();
	if(rows*columns>0){
		matrix_data=std::make_unique<std::complex<double>[]>(rows * arr2d.columns);
		// Copy values into new array
		for(size_t i{1};i<=rows;i++) {
			for(size_t j{1};j<=columns;j++) {
			matrix_data[index(i,j)] = arr2d(i,j);
			}
		}
	}
return *this; // Special pointer!!!
}

// multiplication
matrix matrix::operator*(matrix mat) const
{
	// check dimensions
	if (columns==mat.rows){
		// create sum matrix and initialise data pointer
		matrix product{rows,mat.columns};
		product.matrix_data=std::make_unique<std::complex<double>[]>(rows * mat.columns);
		// incert values into sum matrix
		for(size_t i{1};i<=rows;i++) {
			for(size_t j{1};j<=mat.columns;j++) {
				for(size_t k{1};k<=columns;k++) {
				// perform multiplication 
				product(i,j) = product(i,j) + matrix_data[index(i,k)] * mat(k,j);
				}
			}
		}
		return product;
	} else {std::cerr<<"Error: matricies don't have same dimention of multiplication"<<std::endl;exit(1);}
}


// insertion to output stream
std::ostream & operator<<(std::ostream &os, matrix& mat) //const matrix doesnt work
{
  for(size_t i{1};i<=mat.rows;i++) {
    for(size_t j{1};j<=mat.columns;j++) {
      os<<mat(i,j)<<" ";
    }
    os<<std::endl;
  }
  return os;
}

matrix matrix::tensor_product(const matrix& mat2) const
{
    const size_t rows1 = this->get_rows();
    const size_t cols1 = this->get_cols();
    const size_t rows2 = mat2.get_rows();
    const size_t cols2 = mat2.get_cols();

    matrix result(rows1 * rows2, cols1 * cols2);

    for (size_t i{1}; i<=rows1; i++) {
        for (size_t j{1}; j<=cols1; j++) {
            for (size_t k{1}; k<=rows2; k++) {
                for (size_t l{1}; l<=cols2; l++) {
                    result((i-1) * rows2 + k, (j-1) * cols2 + l) = (*this)(i, j) * mat2(k, l);
                }
            }
        }
    }

    return result;
}

