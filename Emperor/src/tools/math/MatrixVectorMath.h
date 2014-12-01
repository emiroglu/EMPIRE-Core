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
#include "mkl.h"


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



// Classes
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/********//**
 * \brief This is a template class does compressed sparse row matrix computations: CSR Format (3-Array Variation)
 *        http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/GUID-9FCEB1C4-670D-4738-81D2-F378013412B0.htm
 * \author Stefan Sicklinger
 **************************************************************************************************/
template<class T>
class SparseMatrix {
public:
    typedef std::vector<std::map<size_t, T> > mat_t;
    typedef size_t row_iter;
    typedef std::map<size_t, T> col_t;
    typedef typename col_t::iterator col_iter;
    /***********************************************************************************************
     * \brief Constructor for symmetric matrices
     * \param[in] _m is the number of rows & columns
     * \param[in] _isSymmetric only the upper triangular form is stored
     * \author Stefan Sicklinger
     ***********/
    SparseMatrix(const size_t _m, const bool _isSymmetric) {
        m = _m;
        n = _m;
        isSquare = true;
        isSymmetric = _isSymmetric;
        if (!((typeid(T) == typeid(double)) || (typeid(T) == typeid(float)))) {
            assert(0);
        }
        mat = new mat_t(m);
        rowIndex = new std::vector<int>(m + 1);

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
        mat = new mat_t(m);
        rowIndex = new std::vector<int>(m + 1);
    }
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~SparseMatrix() {
        delete mat;
    }
    ;
    /***********************************************************************************************
     * \brief Operator overloaded () for assignment of value e.g. A(i,j)=10
     * \param[in] i is the number of rows
     * \param[in] j is the number of columns
     * \author Stefan Sicklinger
     ***********/
    inline T& operator()(size_t i, size_t j) {
        if (i >= m || j >= n)
            assert(0);
        if (i > j && isSymmetric == true)
            assert(0);
        //not allowed
        return (*mat)[i][j];
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
        if (i > j && isSymmetric == true)
            assert(0);
        //not allowed
        return (*mat)[i][j];
    }
    /***********************************************************************************************
     * \brief Operator overloaded * for Matrix vector multiplication
     * \return std vector
     * \author Stefan Sicklinger
     ***********/
    std::vector<T> operator*(const std::vector<T>& x) { //Computes y=A*x
        if (this->m != x.size())
            assert(0);
        //not allowed

        std::vector<T> y(this->m, 0);
        T sum;
        T sumSym;

        row_iter ii;
        col_iter jj;

        for (ii = 0; ii < m; ii++) {
            sum = 0;
            for (jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++) {
                sum += (*jj).second * x[(*jj).first];

                if ((ii != (*jj).first) && isSymmetric) { //not on the main diagonal
                    //   std::cout << (*ii).first << " ssss "<< (*jj).second <<" tttt "<< x[(*ii).first] << "uuuu" << (*jj).first << std::endl;
                    y[(*jj).first] += (*jj).second * x[ii];
                }

            }
            y[ii] += sum;
        }

        return y;
    }
    /***********************************************************************************************
     * \brief This function is a fast alternative to the operator overloading alternative
     * \param[in] x vector to be multiplied
     * \param[out] y result vector
     * \param[in] elements are the number of entries in the vector
     * \author Stefan Sicklinger
     ***********/
    void mulitplyVec(const T* x, T* y, const size_t elements) { //Computes y=A*x
        if (this->m != elements)
            assert(0);
        //not allowed
        T sum;
        size_t iter;
        for (iter = 0; iter < elements; iter++) {
            y[iter] = 0;
        }

        row_iter ii;
        col_iter jj;

        for (ii = 0; ii < m; ii++) {
            sum = 0;
            for (jj = (*mat)[ii].begin(); jj != (*mat)[ii].end(); jj++) {
                sum += (*jj).second * x[(*jj).first];
                if ((ii != (*jj).first) && isSymmetric) { //not on the main diagonal
                    //   std::cout << (*ii).first << " ssss "<< (*jj).second <<" tttt "<< x[(*ii).first] << "uuuu" << (*jj).first << std::endl;
                    y[(*jj).first] += (*jj).second * x[ii];
                }
            }
            y[ii] += sum;
        }
    }

    /***********************************************************************************************
     * \brief This function is a fast alternative to the operator overloading alternative
     * \param[in] x vector to be multiplied
     * \param[out] y result vector
     * \param[in] elements are the number of entries in the vector
     * \author Chenshen Wu
     ***********/
    void transposeMulitplyVec(const T* x, T* y, const size_t elements) { //Computes y=A*x
        if (this->m != elements)
	    assert(0);
	if (isSymmetric) {
	    mulitplyVec(x, y, elements);
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
    }


/***********************************************************************************************
 * \brief This function analysis and factorize the matrix
 * \author Stefan Sicklinger
 ***********/
    void factorize(){
#ifdef USE_INTEL_MKL
    	this->determineCSR();
        if (isSymmetric) {
            pardiso_mtype = 2;  // real symmetric matrix postive definite matrix
        } else {
            pardiso_mtype = 11; // real and unsymmetric matrix
        }
        // set pardiso default parameters
        pardisoinit(pardiso_pt, &pardiso_mtype, pardiso_iparm);
        pardiso_iparm[2] = 3; //The parallel (OpenMP) version of the nested dissection algorithm
        pardiso_maxfct = 1; // max number of factorizations
        pardiso_mnum = 1; // which factorization to use
        pardiso_msglvl = 0; // do NOT print statistical information
        pardiso_neq = m; // number of rows of C_BB
        pardiso_nrhs = 1; // number of right hand side
        pardiso_phase = 12; // analysis and factorization
        pardiso_error = 0;
        //cout<<"factorizing"<<endl;
        mkl_set_num_threads(1);
        pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
                &pardiso_neq, &values[0], &((*rowIndex)[0]), &columns[0], &pardiso_idum,
                &pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, &pardiso_ddum, &pardiso_ddum,
                &pardiso_error);
        if (pardiso_error != 0) {
            ERROR_OUT() << "Error pardiso factorization failed with error code: " << pardiso_error
                    << std::endl;
            exit(EXIT_FAILURE);
        }
#endif
    }

    /***********************************************************************************************
     * \brief This function performs the prepare of a solution
     * \param[in]  pointer to rhs vector
     * \param[out] pointer to solution vector
     * \return std vector
     * \author Stefan Sicklinger
     ***********/
    void solve(T* x, T* b) { //Computes x=A^-1 *y
#ifdef USE_INTEL_MKL
        // pardiso forward and backward substitution
        pardiso_phase = 33; // forward and backward substitution
        pardiso_error = 0;
        pardiso_iparm[5] = 0; // write solution to b if true otherwise to x
        mkl_set_num_threads(1); // set number of threads to 1 for mkl call only
        pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
                &pardiso_neq, &values[0], &((*rowIndex)[0]), &columns[0], &pardiso_idum,
                &pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, b, x, &pardiso_error);
#endif
    }

    /***********************************************************************************************
     * \brief This function clean Pardiso
     * \author Stefan Sicklinger
     ***********/
    void resetPardiso(){
#ifdef USE_INTEL_MKL
        // clean pardiso
        pardiso_phase = 0; //Release internal memory for L and U matrix number mnum
        pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
                &pardiso_neq, &values[0], &((*rowIndex)[0]), &columns[0], &pardiso_idum,
                &pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, &pardiso_ddum, &pardiso_ddum,
                &pardiso_error);
        pardiso_phase = -1; //Release all internal memory for all matrices
        values.clear();
        columns.clear();
        (*rowIndex).clear();
#endif
    }

    /***********************************************************************************************
     * \brief This function clean Pardiso
     * \author Stefan Sicklinger
     ***********/
    void cleanPardiso(){
#ifdef USE_INTEL_MKL
        // clean pardiso
        pardiso_phase = -1; // deallocate memory
        pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
                &pardiso_neq, &values[0], &((*rowIndex)[0]), &columns[0], &pardiso_idum,
                &pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, &pardiso_ddum, &pardiso_ddum,
                &pardiso_error);
        if (pardiso_error != 0) {
            ERROR_OUT() << "Error deallocation of pardiso failed with error code: " << pardiso_error
                    << std::endl;
            exit(EXIT_FAILURE);
        }
#endif
    }

    /***********************************************************************************************
     * \brief This prints the matrix in CSR style i j value
     * \author Stefan Sicklinger
     ***********/
    void printCSR() {
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
    }
    /***********************************************************************************************
     * \brief This prints the matrix in full style
     * \author Stefan Sicklinger
     ***********/
    void print() {
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

    }
    /***********************************************************************************************
     * \brief This prints the matrix in full style in a file ready to be imported in matlab via dlmread('filename',' ')
     * \author Fabien Pean
     ***********/
    void printToFile(std::string filename) {
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
    }
private:
    /// pointer to the vector of maps
    mat_t* mat;
    /// true if a square matrix should be stored
    bool isSquare;
    /// true if a symmetric matrix should be stored
    bool isSymmetric;
    /// number of rows
    size_t m;
    /// number of columns
    size_t n;
    /// A real array that contains the non-zero elements of a sparse matrix. The non-zero elements are mapped into the values array using the row-major upper triangular storage mapping.
    std::vector<T> values;
    /// Element i of the integer array columns is the number of the column that contains the i-th element in the values array.
    std::vector<int> columns;
    /// Element j of the integer array rowIndex gives the index of the element in the values array that is first non-zero element in a row j.
    std::vector<int>* rowIndex;
    /// pardiso variable
    void *pardiso_pt[64]; // this is related to internal memory management, see PARDISO manual
    /// pardiso variable
    int pardiso_iparm[64];
    /// pardiso variable
    int pardiso_mtype;
    /// pardiso variable
    int pardiso_maxfct;
    /// pardiso variable
    int pardiso_mnum;
    /// pardiso variable
    int pardiso_msglvl;
    /// pardiso variable
    int pardiso_neq;
    /// pardiso variable
    int pardiso_nrhs;
    /// pardiso variable
    int pardiso_phase;
    /// pardiso variable
    double pardiso_ddum;
    /// pardiso variable
    int pardiso_idum;
    /// pardiso variable
    int pardiso_error;
    /***********************************************************************************************
     * \brief This fills the three vectors of the CSR format (one-based)
     * \author Stefan Sicklinger
     ***********/
    void determineCSR() {
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
    }
};



}
}


#endif /* MATRIXVECTORMATH_H_ */
