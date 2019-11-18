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
/***********************************************************************************************//**
 * \file GMRES.h
 * This file holds the class GMRES
 * \date 4/27/2012
 **************************************************************************************************/

#ifndef GMRES_H_
#define GMRES_H_

#include "AbstractCouplingAlgorithm.h"
#include <vector>


namespace EMPIRE {
namespace MathLibrary{
template<typename T>
class SparseMatrix;

template<class T>
class SparseMatrix;
}

class ConnectionIO;

/********//**
 * \brief Class GMRES does GMRES based Co-Simulation
 ***********/
class GMRES: public AbstractCouplingAlgorithm {
public:

	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _name 			The name of the coupling algorithm
	 * \param[in] maxOuterItter		The maximum number of outer iterations for GMRES algorithm
	 * \param[in] maxInnerItter 	The maximum number of inner iterations for GMRES algorithm.
	 * \param[in] residualTolerance The tolerance for the GMRES convergence.
	 * \author Aditya Ghantasala
	 ***********/
	GMRES(std::string _name, int maxOuterItter, int maxInnerItter, double residualTolerance);

	/***********************************************************************************************
	 * \brief Destructor
	 * \author Aditya Ghantasala
	 ***********/
	virtual ~GMRES();

	/***********************************************************************************************
	 * \brief Calculate the new value of the output (improved value for next iteration/time step)
	 *        This is the functions which actually does the GMRES.
	 * \author Aditya Ghantasala
	 ***********/
	void calcNewValue();
	/***********************************************************************************************
	 * \brief Calculate the new value of the Residual
	 * \author Aditya Ghantasala
	 ***********/
	void calcCurrentResidual(){}

	/***********************************************************************************************
	 * \brief Init GMRES
	 * \author Aditya Ghantasala
	 ***********/
	void init();

	/***********************************************************************************************
	 * \brief Adds a Client code to get data
	 * \author Aditya Ghantasala
	 ***********/
	void addInputConnection(ConnectionIO *connection);

	/***********************************************************************************************
	 * \brief Adds a Client code to send data
	 * \author Aditya Ghantasala
	 ***********/
	void addOutputConnection(ConnectionIO *connection);


	/*
	 * May we we need the following functions
	 *
	 * ConstructGlobalResidualVector
	 * SolveInterfaceSystemWithGMRES
	 * ReceiveResidualsFromClients (to be modular)
	 * SendVectorToGetEffectFromClients (to be modular)
	 * ReceiveEffectsFromClients (to be modular)   // Basically we put all the send and receive calls in these function calls.
	 *
	 */

private:

	/// Properties that define the GMRES algorithm

	/// Size of the system to be solved with GMRES.
	/// Basically the total number of degrees of freedom on the interface.
	unsigned long int sysSize;
	/// Number of DOFs on the interface of the system.
	unsigned long int oneSize;
	/// Count variable
	unsigned long int count;
	/// Maximum outer iterations
	unsigned long int maxOuterItter;
	/// No of iterations after which GMRES should restart.
	unsigned long int maxInnerItter;
	/// GMRES convergence tolerance.
	double tolerance;
	/// Solution
	double *solutionUpdate;
	/// GMRES Update
	double *gmresUpdate;
	// RHS vector for the GMRES method
	double *rhs;
	// Boolean to find if the RHS is formulated
	bool rhsDone;


	/// I/O functions for communicating with the clients during the update calculation
	std::vector<ConnectionIO*> functionInput;
	std::vector<ConnectionIO*> functionOutput;


	/// friend class in unit test
	friend class TestGMRES;


	/***********************************************************************************************
	 * \brief Method assembles the global residual vector from the clients for GMRES algorithm
	 * \param[out] globalResidualVec 	-A double pointer to the global residual vector.
	 * 									 Memory should be allocated before passing the pointer.
	 * \author Aditya Ghantasala
	 ***********/
	void constructResidualVector(double *globalResidualVec);

	/***********************************************************************************************
	 * \brief Method gets the effect of the system matrix (for GMRES)
	 * \param[out] effect 			-A double pointer to the global residual vector.
	 * 								 Memory should be allocated before passing the pointer.
	 * \param[in]  vec 				-The vector on which the effect should be obtained.
	 * \author Aditya Ghantasala
	 ***********/
	void getEffectFromClients(double *vec, double *effect);

	/***********************************************************************************************
	 * \brief Method updates the right hand side for the GMRES system
	 * \author Aditya Ghantasala
	 ***********/
	void formulateRHS();


	/***********************************************************************************************
	 * \brief Method calculates the rotations necessary for orthogonalization.
	 * \author Aditya Ghantasala
	 ***********/
	void calculateNextRotations(double, double, double*, double*);

};

} // End of namespace EMPIRE



#endif /* GMRES_H_ */
