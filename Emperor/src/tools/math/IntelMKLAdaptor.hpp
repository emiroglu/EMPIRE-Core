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

#ifndef SPARSELINEARALGEBRA_HPP_
#define SPARSELINEARALGEBRA_HPP_

#include "Message.h"
#include <assert.h>


namespace EMPIRE {
namespace MathLibrary {

// Pre defining the class SparseMatrix so that it can be used here
//class SparseMatrix;

// INCLUDES
#ifdef USE_INTEL_MKL
#include "mkl.h"



/********//**
 * \brief This class acts as the interface to INTEL the linear algebra libraries used. (INTEL MKL)
 * \author Aditya Ghantasala
 **************************************************************************************************/
class PardisoAdapter {

private:
	// ##########################
	// Variables for PARDISO (INTEL)
	// ##########################
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
	/// pardiso variable
	//char pardiso_trans;
	/// pardiso variable
	double pardiso_alpha;
	/// pardiso variable
	char *pardiso_descra;
	/// pardiso variable
	double pardiso_beta;
	/// pardiso variable
	int mklSetNumThreads;
	/// pardiso variable
	char up;



public:

	/***********************************************************************************************
	 * \brief Empty Constructor to the class
	 * \author Aditya Ghantasala
	 ***********/
	PardisoAdapter(bool isSymmetric=false){
		pardiso_descra = "G00F"; // general matrix, indexing from 1
		pardiso_alpha = 1.0;
		pardiso_beta = 0.0;
		mklSetNumThreads = 1;  /// OpenMP parallelization variable
		up = 'u';
		// Initializing the pardiso
		if(isSymmetric)
			pardiso_mtype = 2; // real symmetric matrix
		else
			pardiso_mtype = 11; // Real unsymmetric matrix

		// set pardiso default parameters
		pardisoinit(pardiso_pt, &pardiso_mtype, pardiso_iparm);
	}

	/***********************************************************************************************
	 * \brief Destructor to the class
	 * \author Aditya Ghantasala
	 ***********/
	~PardisoAdapter(){

	}



	/***********************************************************************************************
	 * \brief This function clean Pardiso
     * \param[in]  mat_values 	-- Pointer to arrray containing the non-zero values of sparse matrix.
     * \param[in]  rowIndex		-- Pointer to the integer array of which element j gives the index of the element in the values array that is first non-zero element in a row j.
     * \param[in]  columns		-- Pointer to the integer array of which element i is the number of the column that contains the i-th element in the mat_values array.
     * \author Stefan Sicklinger
	 * \edit Aditya Ghantasala
	 ***********/
	void resetPardiso(double *mat_values, int *rowIndex, int *columns ){
		// clean pardiso
		pardiso_phase = 0; //Release internal memory for L and U matrix number mnum
		pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
				&pardiso_neq, mat_values, rowIndex, columns, &pardiso_idum,
				&pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, &pardiso_ddum, &pardiso_ddum,
				&pardiso_error);
		pardiso_phase = -1; //Release all internal memory for all matrices

	}

	/***********************************************************************************************
	 * \brief This function clean Pardiso
     * \param[in]  mat_values 	-- Pointer to arrray containing the non-zero values of sparse matrix.
     * \param[in]  rowIndex		-- Pointer to the integer array of which element j gives the index of the element in the values array that is first non-zero element in a row j.
     * \param[in]  columns		-- Pointer to the integer array of which element i is the number of the column that contains the i-th element in the mat_values array.
	 * \author Stefan Sicklinger
	 * \edit Aditya Ghantasala
	 ***********/
	void cleanPardiso(double *mat_values, int *rowIndex, int *columns){
        // clean pardiso
        pardiso_phase = -1; // deallocate memory
        pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
                &pardiso_neq, mat_values, rowIndex, columns, &pardiso_idum,
                &pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, &pardiso_ddum, &pardiso_ddum,
                &pardiso_error);
        if (pardiso_error != 0) {
            std::cerr << "Error deallocation of pardiso failed with error code: " << pardiso_error
                    << std::endl;
            exit(EXIT_FAILURE);
        }
	}



    /***********************************************************************************************
     * \brief This function performs the multiplication of give sparse CSR matrix with given vector
     * \param[in]  transpose	-- Bool flag if a transpose of the matrix should be multiplied.
     * \param[in]  m 			-- Number of rows of sparse matrix
     * \param[in]  n 			-- Number of columns of sparse matrix
     * \param[in]  mat_values 	-- Pointer to arrray containing the non-zero values of sparse matrix.
     * \param[in]  rowIndex		-- Pointer to the integer array of which element j gives the index of the element in the values array that is first non-zero element in a row j.
     * \param[in]  columns		-- Pointer to the integer array of which element i is the number of the column that contains the i-th element in the mat_values array.
     * \param[in]  vec			-- Pointer to the array which is to be multiplied to the sparse matrix.
     * \param[in]  isSymmetric	-- Flag mentioning if the sparse matrix is symmetric or not.
     * \param[out] resultVec 	-- Pointer to solution vector
     * \author Aditya Ghantasala
     ***********/
	 void multiply(bool transpose, int m, int n, double *mat_values, int *rowIndex, int* columns, double *vec, double *resultVec, bool isSymmetric=false){
        mkl_set_num_threads(mklSetNumThreads);
        char pardiso_trans;
        if(transpose)
        	pardiso_trans = 'T';
        else
        	pardiso_trans = 'N';

        assert(mat_values != NULL);
        assert(rowIndex!= NULL);
        assert(columns!= NULL);
        assert(vec != NULL);
        assert(resultVec != NULL);
        for(int i = 0; i<m; i++)
        	//std::cout<<"Entry ["<<i<<"] is :: "<<columns[i]<<std::endl;

        if(isSymmetric){
        	//mkl_dcsrsymv(&up, &m, mat_values, rowIndex, columns, vec, resultVec);
        }else{
            //std::cout<<"From Intel Adapter Multiplied :: !! "<<std::endl;
        	//mkl_dcsrmv(&pardiso_trans, &m, &n, &pardiso_alpha, pardiso_descra, mat_values, rowIndex, columns, &columns[1], vec, &pardiso_beta, resultVec);
        }
	 }




		/***********************************************************************************************
		 * \brief This function analysis and factorize the matrix
		 * \param[in]  isSymmetric 	-- Flag specifying if the Sparse matrix is symmetric or not.
	     * \param[in]  m 			-- Number of rows of sparse matrix
	     * \param[in]  mat_values 	-- Pointer to arrray containing the non-zero values of sparse matrix.
	     * \param[in]  rowIndex	 	-- Pointer to the integer array of which element j gives the index of the element in the values array that is first non-zero element in a row j.
	     * \param[in]  columns		-- Pointer to the integer array of which element i is the number of the column that contains the i-th element in the mat_values array.
		 * \author Stefan Sicklinger
		 * \edit Aditya Ghantasala
		 ***********/
		void factorize(bool isSymmetric, int m, double *mat_values, int *rowIndex, int* columns){

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
					&pardiso_neq, mat_values, rowIndex, columns, &pardiso_idum,
					&pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, &pardiso_ddum, &pardiso_ddum,
					&pardiso_error);

			if (pardiso_error != 0) {
				ERROR_OUT() << "Error pardiso factorization failed with error code: " << pardiso_error
						<< std::endl;
				exit(EXIT_FAILURE);
			}
		}






	 /***********************************************************************************************
	  * \brief This function performs the prepare of a solution
      * \param[in]  m 			-- Number of rows of sparse matrix
      * \param[in]  mat_values 	-- Pointer to arrray containing the non-zero values of sparse matrix.
      * \param[in]  rowIndex	-- Pointer to the integer array of which element j gives the index of the element in the values array that is first non-zero element in a row j.
      * \param[in]  columns		-- Pointer to the integer array of which element i is the number of the column that contains the i-th element in the mat_values array.
	  * \param[in]  x 			-- Pointer to rhs vector
	  * \param[out] y 			-- Pointer to solution vector
	  * \return std vector
	  * \author Aditya Ghantasala
	  ***********/
	 void solve(bool isSymmetric, int m, double *mat_values, int *rowIndex, int* columns, double* x, double * b) { //Computes x=A^-1 *b

		 // Factorizing in P and Q
		 factorize(isSymmetric, m, mat_values, rowIndex, columns);
		 // pardiso forward and backward substitution
		 pardiso_phase = 33; // forward and backward substitution
		 pardiso_error = 0;
		 pardiso_iparm[5] = 0; // write solution to b if true otherwise to x // TODO 6 for new version of pardiso and 5 for old
		 mkl_set_num_threads(1); // set number of threads to 1 for mkl call only
		 pardiso(pardiso_pt, &pardiso_maxfct, &pardiso_mnum, &pardiso_mtype, &pardiso_phase,
				 &pardiso_neq, mat_values, rowIndex, columns, &pardiso_idum,
				 &pardiso_nrhs, pardiso_iparm, &pardiso_msglvl, b, x, &pardiso_error);

		 // Checking if the solve is successfull or not.
		 if (pardiso_error != 0) {
			 ERROR_BLOCK_OUT("MortarMapper","conservativeMapping","pardiso solver failed!");
		 }
	 }


};
#endif


}
}

#endif /* SPARSELINEARALGEBRA_HPP_ */
