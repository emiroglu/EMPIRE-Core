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

#include "GMRES.h"
#include "DataField.h"
#include "ConnectionIO.h"
#include "Signal.h"
#include "Residual.h"
#include "MathLibrary.h"

#include <assert.h>
#include <math.h>
#include <sstream>
#include <string.h>


using namespace std;

namespace EMPIRE {

GMRES::GMRES(std::string _name, int maxOuterItter, int maxInnerItter, double tolerance) :
						AbstractCouplingAlgorithm(_name) {

	this->maxOuterItter = maxOuterItter;
	this->maxInnerItter =  maxInnerItter;
	this->tolerance = tolerance;
	this->rhsDone = false;
	this->count = 0;
	this->oneSize = 0;
	this->sysSize = 0;
	this->gmresUpdate = NULL;
	this->rhs = NULL;
	this->solutionUpdate = NULL;
}

GMRES::~GMRES(){
	// Cleaning up the memory
	delete []rhs;
	delete []solutionUpdate;
	delete []gmresUpdate;
}

// TODO :: Check all the indices and initialization of the vectors.
void GMRES::calcNewValue(){

	// Loop for obtaining the GMRES size of the system and computing the residual.

	for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++) {
		this->sysSize += it->second->size;
		oneSize = it->second->size;
		// Computing the residual for the current state.
        it->second->computeCurrentResidual();
	}
	// Adding size for the lambdas
	this->sysSize += oneSize;

	//// Necessary vectors and for GMRES.
	//// Here we are not doing any pre-conditioning for the system

	if(count <= 0){
		// RHS vector
		this->rhs = (double *) calloc(this->sysSize,  sizeof(double));
		// Solution vector
		this->solutionUpdate = (double *) calloc(this->sysSize,  sizeof(double));
		// GMRES update
		gmresUpdate = (double *) calloc(this->sysSize,  sizeof(double));
	}

	// Global residual vector
	double *residual = (double *) calloc(this->sysSize,  sizeof(double));
	// Vector for storing the intermediate orthogonal vectors
	double *w = (double *) calloc(this->sysSize,  sizeof(double));
	// Vector storing the orthogonal Vectors
	double *V = (double *) calloc(this->sysSize * (this->maxInnerItter+1),  sizeof(double));

	// Reduced system system matrix
	EMPIRE::MathLibrary::SparseMatrix<double> H((this->maxInnerItter), (this->maxInnerItter+1));
	// Vectors for storing the rotations
	double *cs = (double *) calloc(this->maxInnerItter,  sizeof(double));
	double *sn = (double *) calloc(this->maxInnerItter,  sizeof(double));
	double *s  = (double *) calloc(this->maxInnerItter+1,  sizeof(double));
	// Vector y
	double *y = (double *) calloc(this->maxInnerItter,  sizeof(double));

	double err = 10.0;

	// Calculating the initial Residue
	constructResidualVector(residual);

	//// Actual GMRES algorithm.
	// GMRES outer Iterations
	int oi;
	for(oi=0; oi<this->maxOuterItter; oi++){ // oi -- Outer Iterations

		// Calculating the norm of initial residue
		double beta = EMPIRE::MathLibrary::vector2norm(residual,this->sysSize);
		s[1] = beta;

		// Defining the initial Krylov vector
		EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(residual,1/beta,this->sysSize);
		// Storing the krylov vector in V
		for(int i=0; i<this->sysSize; i++){
			V[i] = residual[i];
		}

		//// GMRES inner iterations
		// Also includes Gram-Schmidt orthogonalization.
		for(int ii=0; ii<this->maxInnerItter; ii++){ // ii -- Inner Iterations

			getEffectFromClients(residual,w);

			// Orthogonalization
			for (int k=0; k<ii; k++){
				double dum = EMPIRE::MathLibrary::computeDenseDotProduct(this->sysSize,w,&V[this->sysSize*k]);
				H(k,ii) = dum;
				double h_ki = dum;
				EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(&V[this->sysSize*k], h_ki, this->sysSize);
				for(int l=0; l<this->sysSize; l++){
					w[l] = w[l] - V[this->sysSize*k + l];
				}
			}

			// Calculating and Storing the next Krylov vector
			H(ii+1,ii) = EMPIRE::MathLibrary::vector2norm(w,this->sysSize);
			double h_ii = H(ii,ii);
			EMPIRE::MathLibrary::computeDenseVectorMultiplicationScalar(w, 1/h_ii, this->sysSize);
			for(int i=0; i<this->sysSize; i++){
				V[this->sysSize*(ii+1) +i] = w[i];
			}

			// Applying rotations
			for (int k=0; k<ii-1; k++){
				double temp = cs[k]*H(k,ii) + sn[k]*H(k+1,ii);
				H(k+1,ii) = -sn[k]*H(k,ii) + cs[k]*H(k+1,ii);
				H(k,ii) = temp;
			}

			// Calulating the (ii+1)th Rotations
			calculateNextRotations(H(ii,ii), H(ii+1,ii), &cs[ii], &sn[ii]);
			double temp = cs[ii] * s[ii];
			s[ii + 1] = -sn[ii] * s[ii];
			s[ii] = temp;
			H(ii, ii) = cs[ii]*H(ii,ii) + sn[ii]*H(ii+1,ii);
			H(ii+1, ii) = 0.0;
			double err = std::abs(s[ii+1]);

			if(err <= this->tolerance){
				// Before solving, we have to resize the matrix to suit the dimension
				H.resize(0,0,ii,ii);
				// Solve the system and update the solution.
				H.solve(y,s);
				// get update for the update of solution : multiply krylov vector matrix V with y.
				for(unsigned long int n=0; n<this->sysSize; n++)
					for(unsigned long int a=0; a<ii; a++)
						gmresUpdate[n] += V[n*this->sysSize + this->maxInnerItter]*y[a];

				// Applying gmres update to the solution update
				for(unsigned long int n=0; n<this->sysSize; n++)
					solutionUpdate[n] += gmresUpdate[n];

				break;
			}

		}

		if(err <= this->tolerance){
			break;
		}

		// Solve the system and update the solution.
		H.resize(0,0,this->maxInnerItter, this->maxInnerItter);
		H.solve(y,s);
		// get update for the update of solution : multiply krylov vector matrix V with y.
		for(unsigned long int n=0; n<this->sysSize; n++)
			for(unsigned long int a=0; a<this->maxInnerItter; a++)
				gmresUpdate[n] += V[n*this->sysSize + this->maxInnerItter]*y[a];

		// Applying gmres update to the solution update
		for(unsigned long int n=0; n<this->sysSize; n++)
			solutionUpdate[n] += gmresUpdate[n];

		// Calculate residue vector with the new updated solution.
		constructResidualVector(residual);
		s[this->maxInnerItter+1] = EMPIRE::MathLibrary::vector2norm(residual,this->sysSize);

	} // End of GMRES iterations

	/// Information output
    stringstream toOutput;
    toOutput << scientific;
    toOutput << "GMRES iterations converged after ::  "<<oi;
    INDENT_OUT(1, toOutput.str(), infoOut);
    cout.unsetf(ios_base::floatfield);


	//// Updating the output
    //// TODO :: this is wrong. We do not have equal number of residuals and outputs in the case of GMRES.
	assert(outputs.size() == residuals.size());
	int oldResidualSize=0;
	for (map<int, CouplingAlgorithmOutput*>::iterator it = outputs.begin(); it != outputs.end(); it++) {

		CouplingAlgorithmOutput *output = it->second;
		double *newOuput = new double[this->sysSize];

		for(unsigned long int i=0; i<this->oneSize; i++){
			// TODO :: here after all the GMRES iterations, the update is added to the copy at the beginning.
			// There are in total three vectors
			// 1. the outputCopyAtIterationBeginning -> this is from the AbstractCouplingAlgorithm
			// 2. solutionUpdate -> defined in this class -> This will be updated after every GMRES iteration -> this should be added to outputCopyAtIterationBeginning
			// 3. gmresUpdate -> defined in this class -> this will be calculated in the GMRES iterations and will be added to solutionUpdate after every GMRES iteration.

			// TODO :: Now the question is what should this output contain.
			// Should it contain the new velocities/displacements on the
			newOuput[i] = output->outputCopyAtIterationBeginning[i]+this->solutionUpdate[i];  // TODO : Change to actual update
		}
		output->overwrite(newOuput);
	}


	//// Cleaning up the memory
	delete []residual;
	delete []w;
	delete []V;
	delete []cs;
	delete []sn;
	delete []y;
	delete []s;
	count++;
}

void GMRES::init(){
	// Do nothing
}


void GMRES::addInputConnection(ConnectionIO *connection){
	this->functionInput.push_back(connection);
}

void GMRES::addOutputConnection(ConnectionIO *connection){
	this->functionOutput.push_back(connection);
}

void GMRES::constructResidualVector(double *residualVec){
	// formulating the RHS for the GMRES system.
	if(!rhsDone){
		formulateRHS();
		rhsDone = true;
	}

	// Calculating the effect of the the system matrix on the current solution. (Ax)
	double *effect = new double[this->sysSize];
	getEffectFromClients(solutionUpdate, effect);

	// finding the residue
	// r = b - Ax (here Ax is the effect and b is the RHS)
	for(int i=0; i<this->sysSize; i++){
		residualVec[i] = rhs[i] - effect[i];
	}

	// Cleaning up the memory
	delete []effect;
}


void GMRES::formulateRHS(){
	/*
	 * RHS for the GMRES system
	 *
	 * 			 _                     _
	 * 			|  Res_fluid 			|
	 * RHS = 	| 						|
	 * 			|  Res_structure		|
	 * 			|						|
	 * 			|  lambda (u_f - u_s)	|
	 * 			|_					   _|
	 *
	 */

	// Get Fluid and structural solver residue
	// here both functionInput and functionOutput will have same size
	for(int i=0; i<this->functionInput.size(); i++){
		this->functionInput[i]->receive();								//// INCOMING COMMUNICATION
	}


	// Storing the obtained residues into the RHS vector.
	// We have to use them as soon as we receive because we will rewrite them once we receive next entitiy.
	for(int i=0; i<this->functionInput.size(); i++){
		for(int j=0; j<oneSize; j++){
			rhs[i*oneSize+j] = functionInput[i]->dataField->data[j];
		}
	}

	// Get Fluid and Structural solver velocity
	// Here the previous data in the functionInput and functionOutput will be over written.
	for(int i=0; i<this->functionInput.size(); i++){
		this->functionInput[i]->receive();								//// INCOMING COMMUNICATION
	}

	// Calculate lambda as u_f - u_s
	// This is 2 solver FSI specific.
	unsigned long int lambdaPosition = oneSize*functionInput.size();
	for(int i=0; i<oneSize; i++){
		rhs[i+lambdaPosition] = functionInput[functionInput.size() -1]->dataField->data[i] - functionInput[functionInput.size()]->dataField->data[i];
	}
}

void GMRES::getEffectFromClients(double *vec, double *effect){
	// Setting the data in the input
	for(int i=0; i<this->functionInput.size(); i++){
		for(int j=0; j<oneSize; j++){
			functionOutput[i]->dataField->data[j] = vec[i*oneSize+j];
		}
	}
	// Send fluid and structural part of vec to respective solvers
	for(int i=0; i<this->functionInput.size(); i++){
		functionInput[i]->send();										//// OUTGOING COMMUNICATION
	}
	// Receive the effect from fluid and structural solvers.
	for(int i=0; i<this->functionInput.size(); i++){
		functionInput[i]->receive();									//// INCOMING COMMUNICATION
	}
}


void GMRES::calculateNextRotations(double a, double b, double* c, double* s){
	 if ( b == 0.0 ){
	      *c = 1.0;
	      *s = 0.0;
	 }else if ( std::abs(b) > std::abs(a) ){
	      double temp = a / b;
	      *s = 1.0 / std::sqrt( 1.0 + (temp*temp) );
	      *c = temp * (*s);
	 }else{
	      double temp = b / a;
	      *c = 1.0 / std::sqrt( 1.0 + (temp*temp) );
	      *s = temp * (*c);
	 }
}



}

