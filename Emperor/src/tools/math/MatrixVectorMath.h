/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
#ifndef MATRIXVECTORMATH_H_
#define MATRIXVECTORMATH_H_

#include <fstream>
#include <vector>
#include <cstdlib>
#include <map>
#include <vector>
#include <assert.h>
#include <typeinfo>
#include <cmath>
#include "Message.h"
#include "AuxiliaryParameters.h"
// Including Eigen
#ifdef USE_EIGEN
#include "EigenAdapter.h"
// Including intel MKL
#elif USE_INTEL_MKL
#include "IntelMKLAdaptor.hpp"
#endif

namespace EMPIRE {
namespace MathLibrary {

// Variables
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Methods

/***********************************************************************************************
 * \brief Copy dense vector vec1 <- vec2
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \author Stefan Sicklinger
 ***********/
void copyDenseVector(double *vec1, const double *vec2, const int elements);

/***********************************************************************************************
 * \brief Compute Euclidean norm of vector
 * \param[in] vec1 the 1st vector
 * \param[in] elements number of elements in vec1
 * \return Euclidean norm
 * \author Stefan Sicklinger
 ***********/
double vector2norm(const double *vec1, const int elements);

/***********************************************************************************************
 * \brief Computes a vector-scalar product and adds the result to a vector. vec1 <- a*vec1 + vec2
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \param[in] a    scalar
 * \param[in] elements number of elements in vec1
 * \author Stefan Sicklinger
 ***********/
void computeDenseVectorAddition(double *vec1, const double *vec2 ,const double a, const int elements);

/***********************************************************************************************
 * \brief Computes vector scales by scalar vec1 <- vec1*a
 * \param[in] vec1 the 1st vector
 * \param[in] a    scalar
 * \param[in] elements number of elements in vec1
 * \author Stefan Sicklinger
 ***********/
void computeDenseVectorMultiplicationScalar(double *vec1 ,const double a, const int elements);

/***********************************************************************************************
 * \brief Compute the square of the Euclidean distance of two points in n-D space.
 * \param[out] The square of the Euclidean distance between two points in the n-D space
 * \param[in] _length The length of the n-dimensional space
 * \param[in] _Pi The first point
 * \param[in] _Pj The second point
 * \author Andreas Apostolatos
 ***********/
double computeDenseEuclideanNorm(int, double*, double*);

/***********************************************************************************************
 * \brief Compute the dot product between two vectors in the n-D space
 * \param[out] The dot product between two vectors in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vecI The 1st vector
 * \param[in] _vecJ The 2nd vector
 * \author Andreas Apostolatos
 * \edit Aditya Ghantasala (mixing Stefan's implementation)
 ***********/
double computeDenseDotProduct(int, double*, double*);

/***********************************************************************************************
* \brief Compute the cross product between two vectors in the 3-D space
* \param[in/out] _product The product of vector1 and vector 2
* \param[in] _v1 The 1st vector
* \param[in] _v2 The 2nd vector
* \author Chenshen Wu
***********/
void crossProduct(double* _product, double* _v1, double* _v2);

/***********************************************************************************************
 * \brief Solve a 2x2 linear system by close form formula
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[out] b the solution is written to b
 * \return whether the determinant is zero or not
 * \author Tianyang Wang
 ***********/
bool solve2x2LinearSystem(const double *A, double *b, double EPS = 1E-13);

/***********************************************************************************************
 * \brief Solve a 3x3 linear system by close form formula, the system has one row which has all entries 1.
 * \brief Therefore, it cannot solve general 3x3 linear system.
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[out] b the solution is written to b
 * \author Tianyang Wang
 ***********/
void solve3x3LinearSystem(const double *A, int planeToProject, double *b);

/***********************************************************************************************
 * \brief Compute the matrix product between two matrices (for mortar, not so general)
 * \param[in] n n
 * \param[in] m m
 * \param[in] A the matrix with size n*n
 * \param[in] B the matrix with size n*m
 * \param[out] B B is overwritten by A*B (n*m)
 * \author Tianyang Wang
 ***********/
void computeMatrixProduct(int n, int m, const double *A, double *B);

/***********************************************************************************************
 * \brief Compute the matrix product between two matrices (general)
 * \param[in] n The number of rows of matrix A
 * \param[in] p The number of columns of matrix A and the number of rows of matrix B
 * \param[in] m The number of columns of matrix B
 * \param[in] A the matrix with size n*p
 * \param[in] B the matrix with size p*m
 * \param[out] C the product matrix with size n*p
 * \author Andreas Apostolatos
 ***********/
void computeMatrixProduct(int n, int p, int m, const double *A, double *B, double *C);

/***********************************************************************************************
 * \brief Compute the transpose matrix product between two matrices (general)
 * \param[in] p The number of rows of matrix A and B
 * \param[in] n The number of columns of matrix A
 * \param[in] m The number of columns of matrix B
 * \param[in] A the matrix with size p*n
 * \param[in] B the matrix with size p*m
 * \param[out] C The transpose product with size n*m
 * \author Andreas Apostolatos
 ***********/
void computeTransposeMatrixProduct(int p, int n, int m, const double *A, const double *B, double *C);

/***********************************************************************************************
 * \brief My own sparse matrix (csr format, non-symmetric) vector multiplication routine (A * x = y).
 *        The interface is simplified and compatible with mkl_dcsrmv.
 *        This routine supports only one-based indexing of the input arrays.
 * \param[in] trans 'N' --- no transpose, 'T' --- transpose
 * \param[in] numRows number of rows of matrix A
 * \param[in] numCols number of columns of matrix A
 * \param[in] A sparse matrix A in CSR format
 * \param[in] JA JA of A
 * \param[in] IA IA of A
 * \param[in] x vector x
 * \param[in] y vector y
 * \author Tianyang Wang
 ***********/
void dcsrmv(char trans, int numRows, int numCols, const double *A, const int *JA, const int *IA,
        const double *x, double *y);

/***********************************************************************************************
 * \brief My own sparse matrix (csr format, symmetric) vector multiplication routine (A * x = y).
 *        The interface is simplified and compatible with mkl_dcsrsymv.
 *        This routine supports only one-based indexing of the input arrays.
 * \param[in] n size of matrix A
 * \param[in] A sparse matrix A in CSR format (symmetric)
 * \param[in] IA IA of A
 * \param[in] JA JA of A
 * \param[in] x vector x
 * \param[in] y vector y
 * \author Tianyang Wang
 ***********/
void dcsrsymv(int n, const double *A, const int *IA, const int *JA, const double *x, double *y);

/***********************************************************************************************
 * \brief Solve 3x3 Linear system, Ax = b
 * \param[in] _A Square 3x3 matirx
 * \param[in/out] _b Right hand side vector which stores also the solution
 * \param[in] _EPS tolerance up to which matrix _A is regular
 * \param[out] Flag on the solvability of the system
 * \author Chenshen Wu
 ***********/
bool solve3x3LinearSystem(const double* _A, double* _b, double _EPS);

/***********************************************************************************************
 * \brief Computes the Determinant of a given 3x3 matrix
 * \param[in] _A Matrix
 * \param[out] The Determinant of the given matrix
 * \author Chenshen Wu
 ***********/
double det3x3(const double* _A);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/********//**
 * \brief This is a template class does compressed sparse row matrix computations: CSR Format (3-Array Variation)
 *        http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/GUID-9FCEB1C4-670D-4738-81D2-F378013412B0.htm
 * \author Stefan Sicklinger
 * \edit Aditya Ghantasala
 **************************************************************************************************/
template<class T>
class SparseMatrix {
public:
#ifdef USE_INTEL_MKL
    typedef std::vector<std::map<size_t, T> > mat_t;
    typedef size_t row_iter;
    typedef std::map<size_t, T> col_t;
    typedef typename col_t::iterator col_iter;
#endif

    static int clearCount;
    /***********************************************************************************************
     * \brief Constructor for symmetric matrices
     * \param[in] _m is the number of rows & columns
     * \param[in] _isSymmetric only the upper triangular form is stored
     * \author Stefan Sicklinger
     ***********/
    SparseMatrix(const size_t _m, const bool _isSymmetric=false) {
        m = _m;
        n = _m;
        isSquare = true;
        isSymmetric = _isSymmetric;
        isFactorized = false;
        isDetermined = false;
        if (!((typeid(T) == typeid(double)) || (typeid(T) == typeid(float)))) {
            assert(0);
        }

#ifdef USE_INTEL_MKL
        mat = new mat_t(m);
        rowIndex = new std::vector<int>(m + 1);
        intelMKL = new PardisoAdapter(isSymmetric);
#elif USE_EIGEN
        eigenMat = new EigenAdapter<T>(m,n);
        isFull = false;
#endif
    }
    /***********************************************************************************************
     * \brief Constructor for unsymmetric matrices
     * \param[in] _m is the number of rows
     * \param[in] _n is the number of columns
     * \author Stefan Sicklinger
     ***********/
    SparseMatrix(const size_t _m, const size_t _n) {
        m = _m;
        n = _n;
        isSquare = false;
        isSymmetric = false;
        isDetermined = false;
        isFactorized = false;
#ifdef USE_INTEL_MKL
        mat = new mat_t(m);
        rowIndex = new std::vector<int>(m + 1);
        intelMKL = new PardisoAdapter(isSymmetric);
#elif USE_EIGEN
        eigenMat = new EigenAdapter<T>(m,n);
#endif
    }
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~SparseMatrix() {
#ifdef USE_INTEL_MKL
    	//std::cout<<"Cleaning Pardiso"<<std::endl;
    	if(clearCount > 0){
    	 intelMKL->cleanPardiso(&values[0], &((*rowIndex)[0]), &columns[0]);
    	 clearCount++;
    	}

        if(intelMKL != NULL){
        	//std::cout<<"Deleting intelMKL :: "<<intelMKL<<std::endl;
        	//std::cout<<""<<std::endl;
           delete intelMKL;
        }
        delete mat;
        delete rowIndex;
#elif USE_EIGEN
       delete eigenMat;
#endif
    }


    /***********************************************************************************************
     * \brief Getters for the dimensions of the matrix
     * \param[out] number of rows and columns
     * \author Aditya Ghantasala
     ***********/
    inline size_t getNumberOfRows(){return m;};
    inline size_t getNumberOfColumns(){return n;};

    /***********************************************************************************************
     * \brief Operator overloaded () for assignment of value e.g. A(i,j)=10
     * \param[in] i is the number of rows
     * \param[in] j is the number of columns
     * \author Stefan Sicklinger
     ***********/
    inline T& operator()(size_t i, size_t j) {
        if (i >= m || j >= n)
            assert(0);

#ifdef USE_INTEL_MKL
        if (i > j && isSymmetric == true)
            assert(0);
        //not allowed
        isDetermined=false;
        return (*mat)[i][j];
#elif USE_EIGEN
        isDetermined=false;
        return (*eigenMat)(i,j);
#endif
    }
    /***********************************************************************************************
     * \brief Operator overloaded () for assignment of value e.g. A(i,j)=10
     * \param[in] i is the number of rows
     * \param[in] j is the number of columns
     * \author Stefan Sicklinger
     ***********/
    inline T operator()(size_t i, size_t j) const {
        if (i >= m || j >= n)
            assert(0);
#ifdef USE_INTEL_MKL
        if (i > j && isSymmetric == true)
            assert(0);
        //not allowed
        return (*mat)[i][j];
#elif USE_EIGEN
        isDetermined=false;
        return (*eigenMat)(i,j);
#endif
    }


    /***********************************************************************************************
     * \brief This fills the three vectors of the CSR format (one-based)
     * \author Stefan Sicklinger
     ***********/
    void determineCSR() {
    	// Check if this function is already called once.
    	if(isDetermined)
    		return;

#ifdef USE_INTEL_MKL
    	columns.clear();
    	values.clear();
    	row_iter ii;
    	col_iter jj;
    	size_t ele_row = 0; //elements in current row
    	std::cout << std::scientific;

    	for (ii = 0; ii < m; ii++) {
    		(*rowIndex)[ii] = (ele_row + 1);
    		for (jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++) {
    			columns.push_back(((*jj).first) + 1);
    			values.push_back((*jj).second);
    			ele_row++;
    		}

    	}
    	(*rowIndex)[m] = (ele_row + 1);

    	// Tell the matrix that it is determined once.
    	isDetermined = true;
#elif USE_EIGEN
    	makeFullMatrix();
    	isFull = true;
    	eigenMat->determineCSR();
    	// Tell the matrix that it is determined once.
    	isDetermined = true;
#endif
    }


    /***********************************************************************************************
     * \brief This function is a fast alternative to the operator overloading alternative
     * \param[in] 	-- transpose 	Bool flag specifying if a transpose of the matrix should be multiplied or not.
     * \param[in] 	-- x 			vector to be multiplied
     * \param[out] 	-- y 			result vector
     * \param[in] 	-- elements 	are the number of entries in the resultant vector (number of rows of the matrix m)
     * \author Stefan Sicklinger
     * \edit   Aditya Ghantasala
     ***********/
    void mulitplyVec(bool transpose, T* vec, T* resultVec, size_t elements) { //Computes y=A*x

    	// Checking the possibility of multiplication.
    	assert(elements == m);
    	assert(vec != NULL);
    	assert(resultVec != NULL);

    	// Formulating the vectors of the sparse matrix
    	determineCSR();

#ifdef USE_INTEL_MKL
    	T sum;
    	size_t iter;
    	for (iter = 0; iter < elements; iter++) {
    		resultVec[iter] = 0;
    	}

    	row_iter ii;
    	col_iter jj;

    	for (ii = 0; ii < m; ii++) {
    		sum = 0;
    		for (jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++) {
    			sum += (*jj).second * vec[(*jj).first];
    			if ((ii != (*jj).first) && isSymmetric) { //not on the main diagonal
    				//   std::cout << (*ii).first << " ssss "<< (*jj).second <<" tttt "<< x[(*ii).first] << "uuuu" << (*jj).first << std::endl;
    				resultVec[(*jj).first] += (*jj).second * vec[ii];
    			}
    		}
    		resultVec[ii] = sum;
    	}
#elif USE_EIGEN
    	eigenMat->mulitplyVec(false, vec, resultVec, elements);
#endif
    }

    /***********************************************************************************************
     * \brief This function is a fast alternative to the operator overloading alternative
     * \param[in] 	-- x 			Vector to be multiplied
     * \param[out] 	-- y 			Result vector
     * \param[in] 	-- elements 	are the number of entries in the resultant vector (number of rows of the matrix m)
     * \author Chenshen Wu
     ***********/
    void transposeMulitplyVec(T* x, T* y, const size_t elements) { //Computes y=A*x

    	assert(x != NULL);
    	assert(y != NULL);
    	assert(elements >= 0);

    	// Formulating the vectors of the sparse matrix
    	determineCSR();
#ifdef USE_INTEL_MKL
    	if (this->m != elements)
    		assert(0);
    	if (isSymmetric) {
    		mulitplyVec(true, x, y, elements);
    		return;
    	}

    	size_t iter;
    	for (iter = 0; iter < this->n; iter++) {
    		y[iter] = 0;
    	}

    	row_iter ii;
    	col_iter jj;

    	for (ii = 0; ii < m; ii++)
    		for (jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++)
    			y[(*jj).first] += (*jj).second * x[ii];
#elif USE_EIGEN
    	eigenMat->mulitplyVec(true, x, y, elements);
#endif

    }


    /***********************************************************************************************
     * \brief This function returns the sum of a requested row of the sparse matrix
     * \param[in] 	-- row 			Row number of the sparse matrix for which sum should be obtained
     * \param[out] 	-- sum 			Sum of the row'th row.
     * \author Aditya Ghantasala
     ***********/
    T getRowSum(size_t row) { //Computes y=A*x
    	// Checking if the row requested is with in the limits
    	assert(row <= this->m);
    	T sum = 0.0;

    	// Formulating the vectors of the sparse matrix
    	determineCSR();

#ifdef USE_INTEL_MKL
    	if (isSymmetric) {
    		// TODO Check how to return the right value for a symmetric matrix.
    	}
    	col_iter jj;

   		for (jj = (*mat)[row].begin(); jj != (*mat)[row].end(); jj++)
    		sum += (*jj).second;
#elif USE_EIGEN
   		eigenMat->getRowSum(row);
#endif
   		return sum;
    }

    bool isRowEmpty(size_t row) {
#ifdef USE_INTEL_MKL
    	if((*mat)[row].begin() == (*mat)[row].end())
    		return true;
		return false;
#elif USE_EIGEN
		return eigenMat->isRowEmpty(row);
#endif
    }


    /***********************************************************************************************
     * \brief This function deletes or resets a whole row in sparse matrix.
     * \param[in] 	-- row 			Row number of the sparse matrix for which should be deleted.
     * \author Aditya Ghantasala
     * \ edit Altug Emiroglu : when a row is deleted isFactorized flag is set to false
     ***********/
    void deleteRow(size_t row){
#ifdef USE_INTEL_MKL
    	(*mat)[row].clear();
#elif USE_EIGEN
    	eigenMat->deleteRow(row);
#endif
    }


//    /***********************************************************************************************
//     * \brief This function resizes the sparse matrix to given sizes
//     * \param[in] startRow 			- Start row of the sub matrix required.
//     * \param[in] startCol 			- Start Column of the sub matrix required.
//     * \param[in] numRows 			- Number of rows from the startRow.
//     * \param[in] numColumns 		- Number of columns from the startCol
//     * \author Aditya Ghantasala
//     ***********/
//    void resize(long int startRow, long int startCol, long int numRows, long int numColumns) {
//
//    	for(int i=startRow; i<startRow+numRows; i++){
//    		for(int j= startCol; j<startCol+numColumns; j++){
//    			col_iter jj;
//    			jj = (*mat)[i].find(j);
//    			(*mat)[i].erase(jj);
//    		}
//    	}
//    }


    /***********************************************************************************************
     * \brief This function multiplies the whole row of the sparse matrix with a given number
     * \param[in] 	-- row 			Row number of the sparse matrix for which should be multiplied with.
     * \param[in] 	-- face 		factor with which the row should be multiplied.
     * \author Aditya Ghantasala
     ***********/
    void multiplyRowWith(size_t row, T fact) { //Computes y=A*x
    	// Checking if the row requested is with in the limits
    	assert(row <= this->m);
    	// Formulating the vectors of the sparse matrix
    	determineCSR();

#ifdef USE_INTEL_MKL
    	if (isSymmetric) {
    		// TODO Check how to return the right value for a symmetric matrix.
    	}
    	col_iter jj;
    	T dum;
   		for (jj = (*mat)[row].begin(); jj != (*mat)[row].end(); jj++){
    		dum = (*jj).second;
    		(*jj).second = dum * fact;
   		}
#elif USE_EIGEN
   		eigenMat->multiplyRowWith(row, fact);
#endif

    }

    /***********************************************************************************************
     * \brief This function performs the prepare of a solution
     * \param[in]  	-- x pointer to solution vector
     * \param[out] 	-- b pointer to rhs vector
     * \return std vectordetermineCSR
     * \author Stefan Sicklinger
     * \ȩdit Aditya Ghantasala
     ***********/
    void solve(T* x, T* b) { //Computes x=A^-1 *b
	
	assert(x != NULL);
	assert(b != NULL);

    	if(!isFactorized){
#ifdef USE_INTEL_MKL
    		factorize();
#elif USE_EIGEN
    		eigenMat->factorize();
#endif
    		isFactorized = true;
    	}

    	// Constructing the sparse matrix entities
    	determineCSR();
#ifdef USE_INTEL_MKL
    	intelMKL->solve(isSymmetric, m, &values[0], &((*rowIndex)[0]), &columns[0], x, b);
#elif USE_EIGEN
    	eigenMat->solve(x, b);
#endif
    }

    /***********************************************************************************************
     * \brief This function factorizes and prepares for a solution
     * \author Stefan Sicklinger
     * \ȩdit Aditya Ghantasala
     ***********/
    void factorize() {
    	// Constructing the sparse matrix entities
#ifdef USE_INTEL_MKL
   	determineCSR();
  	intelMKL->factorize(isSymmetric, m, &values[0], &((*rowIndex)[0]), &columns[0]);
#elif USE_EIGEN
    determineCSR();
  	eigenMat->factorize();
#endif

    }


    /***********************************************************************************************
     * \brief This function resizes the sparse matrix to given sizes
     * \param[in] startRow 			- Start row of the sub matrix required.
     * \param[in] startCol 			- Start Column of the sub matrix required.
     * \param[in] numRows 			- Number of rows from the startRow.
     * \param[in] numColumns 		- Number of columns from the startCol
     * \author Aditya Ghantasala
     ***********/
    void resize(long int startRow, long int startCol, long int numRows, long int numColumns) {
    	// Constructing the sparse matrix entities
#ifdef USE_INTEL_MKL
    	for(int i=startRow; i<startRow+numRows; i++){
    		for(int j= startCol; j<startCol+numColumns; j++){
    			col_iter jj;
    			jj = (*mat)[i].find(j);
    			(*mat)[i].erase(jj);
    		}
    	}
#elif USE_EIGEN
  	eigenMat->resize(startRow, startCol, numRows, numColumns);
#endif
    }


	/***********************************************************************************************
	 * \brief This function clean Pardiso
	 * \author Stefan Sicklinger
	 * \edit Aditya Ghantasala
	 ***********/
    void reset(){
#ifdef USE_INTEL_MKL
    	intelMKL->resetPardiso(&values[0], &((*rowIndex)[0]), &columns[0] );
		values.clear();
		columns.clear();
		(*rowIndex).clear();
#elif USE_EIGEN
		// TODO Checkk what to do in case of Eigen
#endif
    }

    /***********************************************************************************************
     * \brief This prints the matrix in CSR style i j value
     * \author Stefan Sicklinger
     * \edit Aditya Ghantasala
     ***********/
    void printCSR() {
#ifdef USE_INTEL_MKL
        row_iter ii;
        col_iter jj;
        size_t ele_row; //elements in current row
        std::cout << std::scientific;

        for (ii = 0; ii < m; ii++) {
            for (jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++) {
                std::cout << ii << ' ';
                std::cout << (*jj).first << ' ';
                std::cout << (*jj).second << std::endl;
            }
        }
        std::cout << std::endl;
#elif USE_EIGEN
       print();
#endif

    }

    /***********************************************************************************************
     * \brief This returns the flag on if the factorization is performed
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    inline bool isFactorization() { return isFactorized; }

    /***********************************************************************************************
     * \brief This sets the flag on the factorization
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    inline void setFactorization(bool _flag) { isFactorized = _flag; }

    /***********************************************************************************************
     * \brief This prints the matrix in full style
     * \author Stefan Sicklinger
     ***********/
    void print() {

#ifdef USE_INTEL_MKL
        size_t ii_counter;
        size_t jj_counter;

        std::cout << std::scientific;
        for (ii_counter = 0; ii_counter < m; ii_counter++) {
            for (jj_counter = 0; jj_counter < n; jj_counter++) {
                if ((*mat)[ii_counter].find(jj_counter) != (*mat)[ii_counter].end()) {
                    std::cout << '\t' << (*mat)[ii_counter].find(jj_counter)->second;
                } else {
                    std::cout << '\t' << 0.0;
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
#elif USE_EIGEN
        eigenMat->print();
#endif


    }
    /***********************************************************************************************
     * \brief This prints the matrix in full style in a file ready to be imported in matlab via dlmread('filename',' ')
     * \author Fabien Pean
     ***********/

    void printFullToFile(std::string filename) {
#ifdef USE_INTEL_MKL
        size_t ii_counter;
        size_t jj_counter;

        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ofstream::out);
        ofs << std::scientific;
        for (ii_counter = 0; ii_counter < m; ii_counter++) {
            for (jj_counter = 0; jj_counter < n; jj_counter++) {
            	if(jj_counter!=0) ofs<<" ";
                if(isSymmetric) {
                    if(ii_counter<=jj_counter) {
                        ofs<<(this->operator()(ii_counter,jj_counter));
                    } else {
                        ofs<<(this->operator()(jj_counter,ii_counter));
                    }
                }else{
                    ofs<<(this->operator()(ii_counter,jj_counter));
                }
            }
            ofs<<std::endl;
        }
        ofs.close();
#elif USE_EIGEN
        eigenMat->printToFile(filename);
#endif
    }
    /***********************************************************************************************
     * \brief This prints the matrix in sparse format (i,j)->v in a file.
     *  	The offset can be set to transfer from C-indexing[0] to Matlab[1]
     *  	Can be read in Matlab by M=dlmread('filename'); followed by M=spconvert(M);
     * \author Fabien Pean
     ***********/
    void printCSRToFile(std::string filename, int offset=0) {
#ifdef USE_INTEL_MKL
        std::ofstream ofs;
        ofs.open(filename.c_str(), std::ofstream::out);
        ofs << std::scientific;
        for (row_iter ii = 0; ii < m; ii++) {
            for (col_iter jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++) {
                ofs << ii + offset << ' ';
                ofs << (*jj).first + offset << ' ';
                ofs << (*jj).second << std::endl;
                if(isSymmetric && ii != (*jj).first) {
                	ofs << (*jj).first + offset << ' ';
                	ofs << ii + offset << ' ';
                	ofs << (*jj).second << std::endl;
                }
            }
        }
        // Print last value if zero anyway to get right matrix size exported
        if((*mat)[m-1][n-1] == 0) {
            ofs << m-1 + offset << ' ';
            ofs << n-1 + offset << ' ';
            ofs << (*mat)[m-1][n-1] << std::endl;
        }
        ofs << std::endl;
        ofs.close();
#elif use_EIGEN
        eigenMat->printToFile(filename);
#endif
    }

private:
    /// true if a square matrix should be stored
    bool isSquare;
    /// true if a symmetric matrix should be stored
    bool isSymmetric;
    /// true if the determineCSR function is used once
    bool isDetermined;
    /// true if the matrix is factorized by intel
    bool isFactorized;

    /// number of rows
    size_t m;
    /// number of columns
    size_t n;

#ifdef USE_INTEL_MKL
    /// pointer to the vector of maps
    mat_t* mat;
    /// A real array that contains the non-zero elements of a sparse matrix. The non-zero elements are mapped into the values array using the row-major upper triangular storage mapping.
    std::vector<T> values;
    /// Element i of the integer array columns is the number of the column that contains the i-th element in the values array.
    std::vector<int> columns;
    /// Element j of the integer array rowIndex gives the index of the element in the values array that is first non-zero element in a row j.
    std::vector<int>* rowIndex;

    /// Object to access intel MKL functions
    PardisoAdapter* intelMKL;
#elif USE_EIGEN
    EigenAdapter<T>* eigenMat;
    bool isFull;
#endif

#ifdef USE_EIGEN
    void makeFullMatrix(){
    	if(!isFull){
    	if(isSymmetric){
    		for(int i=0; i<m; i++){
    			for(int j=0; j<n; j++){
    				(*eigenMat)(j,i) = (*eigenMat)(i,j);
    			}
    		}
    	}
    	}
    }

#endif


};

template<class T>
int EMPIRE::MathLibrary::SparseMatrix<T>::clearCount = 0;


}
}


#endif /* MATRIXVECTORMATH_H_ */
