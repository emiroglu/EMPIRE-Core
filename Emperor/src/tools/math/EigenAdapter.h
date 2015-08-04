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


#ifndef EIGENADAPTER_H_
#define EIGENADAPTER_H_

// Including the Eigen header
#include <Sparse>
namespace EMPIRE {
namespace MathLibrary {


/********//**
 * \brief This class acts as the interface to INTEL the linear algebra libraries used. (INTEL MKL)
 * \author Aditya Ghantasala
 **************************************************************************************************/
template<class T>
class EigenAdapter {

	// Defining the definition of sparse matrix.
	typedef Eigen::SparseMatrix<T,Eigen::RowMajor> SpMat; // declares a column-major sparse matrix type of double
	//typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

#ifdef EIGEN_ITERATIVE
	typedef Eigen::BiCGSTAB<SpMat> SpCgSolver;
#else
	typedef Eigen::SparseQR< SpMat, Eigen::NaturalOrdering<int> > SpSolver;
#endif();

private:
	// Number of rows
	int m;
	// Number of columns
	int n;

	// Sparse matrix from Eigen.
	SpMat* A;
	
	// Bool values to save the state
	bool isFactorized, isAnalyzed, isCompressed;

#ifdef EIGEN_ITERATIVE
	SpCgSolver *solver;
#else
	// QR Solver for the sparse matrix system.
	// QR is chosen to accommodate rectangular matrix.
	SpSolver *solver;
#endif

public:
	/***********************************************************************************************
	 * \brief Constructor for unsymmetric matrices
	 * \param[in] _m is the number of rows
	 * \param[in] _n is the number of columns
	 * \author Aditya Ghantasala
	 ***********/
	EigenAdapter(size_t i_m, size_t i_n):m(i_m),n(i_n){
		A = new SpMat(m,n);
		determineCSR(); // The matrix should be in compressed mode for the solver to initiate.
#ifdef EIGEN_ITERATIVE
		solver = new SpCgSolver();
#else
		solver = new SpSolver((*A));
#endif
		isFactorized = false;
		isAnalyzed = false;
		isCompressed = false;
	}

    /***********************************************************************************************
     * \brief Destructor
     * \author Aditya Ghantasala
     ***********/
    virtual ~EigenAdapter() {
    	delete A;
    	delete solver;
    }

	/***********************************************************************************************
	 * \brief Operator overloaded () for assignment of value e.g. A(i,j)=10
	 * \param[in] i is the number of rows
	 * \param[in] j is the number of columns
	 * \author Aditya Ghantasala
	 ***********/
	inline T& operator()(size_t i, size_t j) {
		return A->coeffRef(i,j);
	}

	/***********************************************************************************************
	 * \brief Operator overloaded () for assignment of value e.g. k = A(i,j)
	 * \param[in] i is the number of rows
	 * \param[in] j is the number of columns
	 * \author Aditya Ghantasala
	 ***********/
	inline T& operator()(size_t i, size_t j) const {
		return A->coeffRef(i,j);
	}


	/***********************************************************************************************
	 * \brief Method to make the the CSR format of the matrix. This creates necessary data for Eigen.
	 * \author Aditya Ghantasala
	 ***********/
	void determineCSR(){
		A-> makeCompressed();
		isCompressed = true;
	}


	/***********************************************************************************************
	 * \brief Method resizes the sparse matrix with given start and end col and rows
	 * \author Aditya Ghantasala
	 ***********/
	void resize(long int startRow, long int startCol, long int numRows, long int numColumns){

	}


	/***********************************************************************************************
	 * \brief This function is a fast alternative to the operator overloading alternative
	 * \param[in] 	-- transpose 	Bool flag specifying if a transpose of the matrix should be multiplied or not.
	 * \param[in] 	-- x 			Vector to be multiplied
	 * \param[out] 	-- y 			Result vector
	 * \param[in] 	-- elements 	are the number of entries in the resultant vector (number of rows of the matrix m)
	 * \edit   Aditya Ghantasala
	 ***********/
	void mulitplyVec(bool transpose, T* vec, T* resultVec, size_t elements) { //Computes y=A*x
		assert(elements == m);
		SpMat * b;
		SpMat * y;	

		
		// TODO This is a very naive implementation. Improve it.

		if(transpose == false){
			b = new SpMat(n,1); // Vector to multiply.
			for(int i=0; i<n; i++){
				(*b).coeffRef(i,0) = vec[i];
			}

			y = new SpMat(m,1); //

			(*y) = (*A) * (*b);

			for(int i=0; i<m; i++){
				resultVec[i] = (*y).coeff(i,0);
			}

			delete b;
			delete y;
			return;
		}else{

			b = new SpMat(m,1); // Vector to multiply.
			for(int i=0; i<m; i++){
				(*b).coeffRef(i,0) = vec[i];
			}

			y = new SpMat(n,1); //			

			(*y) = SpMat(A->transpose()) * (*b);
			for(int i=0; i<n; i++){
				resultVec[i] = (*y).coeff(i,0);
			}

			delete b;
			delete y;
			return;
		}
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
    }

	/***********************************************************************************************
	 * \brief This function returns the sum of a requested row of the sparse matrix
	 * \param[in] 	-- row 			Row number of the sparse matrix for which sum should be obtained
	 * \param[out] 	-- sum 			Sum of the row'th row.
	 * \author Aditya Ghantasala
	 ***********/
	T getRowSum(size_t row){
		T sum = 0;

		for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it((*A), row); it; ++it){
			sum += it.value();
		}

		return sum;
	}

	/***********************************************************************************************
         * \brief This function returns boolean whether row is empty or not
         * \param[in]   -- row                  Row number of the sparse matrix to check
         * \author Fabien Pean
         ***********/
	bool isRowEmpty(size_t row) {
		return A->row(row).norm()==0;
	}

	/***********************************************************************************************
	 * \brief This function deletes or resets a whole row in sparse matrix.
	 * \param[in] 	-- row 			Row number of the sparse matrix for which should be deleted.
	 * \author Aditya Ghantasala
	 ***********/
	void deleteRow(size_t row){
		for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it((*A), row); it; ++it){
			A->coeffRef((int)it.row(),(int)it.col()) = 0;
		}
		A->makeCompressed();
	}

    /***********************************************************************************************
     * \brief This function multiplies the whole row of the sparse matrix with a given number
     * \param[in] 	-- row 			Row number of the sparse matrix for which should be multiplied with.
     * \param[in] 	-- fact 		Factor with which the row should be multiplied.
     * \author Aditya Ghantasala
     ***********/
	void multiplyRowWith(size_t row, T fact) {
		for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it((*A), row); it; ++it){
			A->coeffRef((int)it.row(),(int)it.col()) = it.value() * fact;
		}
	}

    /***********************************************************************************************
     * \brief This function performs the prepare of a solution
     * \param[in]  	-- x pointer to solution vector
     * \param[out] 	-- b pointer to rhs vector
     * \ȩdit Aditya Ghantasala
     ***********/
	void solve(T* x, T* b) { //Computes x=A^-1 *b
#ifdef EIGEN_ITERATIVE
        if(!isCompressed)
            determineCSR();

        if(!isFactorized) {
            solver->compute(*A);
            isFactorized = 1;
        }

        if(!isAnalyzed)
            solver->analyzePattern((*A));
		//Converting into Eigen format
		Eigen::VectorXd RHS(m);

		for(int i=0; i<m; i++)
			RHS(i) = b[i];

		Eigen::VectorXd x1 = solver->solve(RHS);

		for(int j=0; j<n; j++)
			x[j] = x1(j);
#else
	if(!isCompressed)
		determineCSR();

	if(!isFactorized)
	  factorize();

	if(!isAnalyzed)
	  solver->analyzePattern((*A));


    	//Converting into Eigen format
    	Eigen::VectorXd RHS(m);

    	for(int i=0; i<m; i++)
    		RHS(i) = b[i];

    	Eigen::VectorXd x1 = solver->solve(RHS);

    	for(int j=0; j<n; j++)
    		x[j] = x1(j);

#endif
    }


    /***********************************************************************************************
     * \brief This function factorizes and prepares for a solution
     * \author Aditya Ghantasala
     ***********/
    void factorize() {
    	// Compute the ordering permutation vector from the structural pattern of A
    	solver->analyzePattern((*A));
    	// Compute the numerical factorization
    	solver->factorize((*A));
    	isFactorized = true;
    	isAnalyzed = true;
    }


    void print(){
    	std::cout<<(*A)<<std::endl;
    }

    void printToFile(std::string filename){
    	std::ofstream ofs;
    	ofs.open(filename.c_str(), std::ofstream::out);
    	ofs << std::scientific;
    	ofs<<(*A);
    	ofs<<std::endl;
    	ofs.close();
    }

};

}
}
#endif /* EIGENADAPTER_H_ */
