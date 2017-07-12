/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
// inclusion of standard libraries   (only if really necessary here in *.h)
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <iomanip>

// Inclusion of user-defined libraries
#include "IGAPatchSurface.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class IGAPatchSurface
 ***********/

class TestIGAPatchSurface: public CppUnit::TestFixture {

private:
	IGAPatchSurface* theIGAPatchSurface;
	double Tol;
	double relTol;
	double TolDeriv;

public:
	void setUp() {
		// Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
        Tol = 1e-15;

		// Assign a relaxed tolerance value (corresponding to maximum accuracy provided by MATLAB) for the Newton-Rapson iteration error (accumulative error appears here)
        relTol = Tol*10;

		// Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the derivative functional values
        TolDeriv = relTol*10;

		// Provide an id for the basis
		int id_basis = 3;

		// Provide an id for the patch itself
		int id_patch = 1;

		// The polynomial degrees
		int p = 4;
		int q = 3;

		// Number of knots in both directions
		int uNoKnots = 15;
		int vNoKnots = 13;

		// The knot vectors in each direction
		double* uKnotVector = new double[uNoKnots];
		for (int i = 0; i < 5; i++)
			uKnotVector[i] = -4;
		for (int i = 5; i < 7; i++)
			uKnotVector[i] = 5;
		for (int i = 7; i < 10; i++)
			uKnotVector[i] = 11;
		for (int i = 10; i < 15; i++)
			uKnotVector[i] = 21;
		double* vKnotVector = new double[vNoKnots];
		for (int i = 0; i < 4; i++)
			vKnotVector[i] = -11;
		for (int i = 4; i < 7; i++)
			vKnotVector[i] = -0.1;
		for (int i = 7; i < 9; i++)
			vKnotVector[i] = 0.0;
		for (int i = 9; i < 13; i++)
			vKnotVector[i] = 1.0;

		// The Control Point net
		int uNoControlPoints = uNoKnots - p - 1;
		int vNoControlPoints = vNoKnots - q - 1;
		IGAControlPoint** controlPointNet = new IGAControlPoint*[uNoControlPoints
				* vNoControlPoints];
//		double* controlPointWeights = new double[uNoControlPoints
//				* vNoControlPoints];

// Control Points for a NURBS

// First row
		controlPointNet[0] = new IGAControlPoint(0, 0.0, -1.0, 4.500000000000000,
				1.0);
		controlPointNet[1] = new IGAControlPoint(1, 0.0, -0.500000000000000,
				4.500000000000000, 3.0);
		controlPointNet[2] = new IGAControlPoint(2, 0.0, -0.250000000000000,
				4.500000000000000, 1.0);
		controlPointNet[3] = new IGAControlPoint(3, 0.0, -0.125000000000000,
				4.500000000000000, 4.0);
		controlPointNet[4] = new IGAControlPoint(4, 0.0, 0.0, 4.500000000000000,
				5.0);
		controlPointNet[5] = new IGAControlPoint(5, 0.0, 0.125000000000000,
				4.500000000000000, 1.200000000000000);
		controlPointNet[6] = new IGAControlPoint(6, 0.0, 0.250000000000000,
				4.500000000000000, 1.0);
		controlPointNet[7] = new IGAControlPoint(7, 0.0, 0.500000000000000,
				4.500000000000000, 3.400000000000000);
		controlPointNet[8] = new IGAControlPoint(8, 0.0, 1.0, 4.500000000000000,
				5.600000000000000);

		// Second row
		controlPointNet[9] = new IGAControlPoint(9, 0.450000000000000, -1.0,
				4.500000000000000, 0.560000000000000);
		controlPointNet[10] = new IGAControlPoint(10, 0.450000000000000,
				-0.500000000000000, 4.500000000000000, 1.340000000000000);
		controlPointNet[11] = new IGAControlPoint(11, 0.450000000000000,
				-0.250000000000000, 4.500000000000000, 1.001000000000000);
		controlPointNet[12] = new IGAControlPoint(12, 0.450000000000000,
				-0.125000000000000, 4.500000000000000, 4.890000000000000);
		controlPointNet[13] = new IGAControlPoint(13, 0.450000000000000, 0.0,
				4.500000000000000, 5.210000000000000);
		controlPointNet[14] = new IGAControlPoint(14, 0.450000000000000,
				0.125000000000000, 4.500000000000000, 0.200000000000000);
		controlPointNet[15] = new IGAControlPoint(15, 0.450000000000000,
				0.250000000000000, 4.500000000000000, 1.0);
		controlPointNet[16] = new IGAControlPoint(16, 0.450000000000000,
				0.500000000000000, 4.500000000000000, 3.400000000000000);
		controlPointNet[17] = new IGAControlPoint(17, 0.450000000000000, 1.0,
				4.500000000000000, 5.600000000000000);

		// Third row
		controlPointNet[18] = new IGAControlPoint(9, 0.500000000000000, -1.0,
				4.500000000000000, 2.6);
		controlPointNet[19] = new IGAControlPoint(10, 0.500000000000000,
				-0.500000000000000, 4.500000000000000, 0.3);
		controlPointNet[20] = new IGAControlPoint(11, 0.500000000000000,
				-0.250000000000000, 4.500000000000000, 1.2);
		controlPointNet[21] = new IGAControlPoint(12, 0.500000000000000,
				-0.125000000000000, 4.500000000000000, 5.8);
		controlPointNet[22] = new IGAControlPoint(13, 0.500000000000000, 0.0,
				4.500000000000000, 3.2);
		controlPointNet[23] = new IGAControlPoint(14, 0.500000000000000,
				0.125000000000000, 4.500000000000000, 6.0);
		controlPointNet[24] = new IGAControlPoint(15, 0.500000000000000,
				0.250000000000000, 4.500000000000000, 1.0);
		controlPointNet[25] = new IGAControlPoint(16, 0.500000000000000,
				0.500000000000000, 4.500000000000000, 1.4);
		controlPointNet[26] = new IGAControlPoint(17, 0.500000000000000, 1.0,
				4.500000000000000, 5.0);

		// Fourth row
		controlPointNet[27] = new IGAControlPoint(27, 0.562500000000000, -1.0,
				4.500000000000000, 1.0);
		controlPointNet[28] = new IGAControlPoint(28, 0.562500000000000,
				-0.500000000000000, 4.500000000000000, 3.0);
		controlPointNet[29] = new IGAControlPoint(29, 0.562500000000000,
				-0.250000000000000, 4.500000000000000, 1.0);
		controlPointNet[30] = new IGAControlPoint(30, 0.562500000000000,
				-0.125000000000000, 4.500000000000000, 4.0);
		controlPointNet[31] = new IGAControlPoint(31, 0.562500000000000, 0.0,
				4.500000000000000, 5.0);
		controlPointNet[32] = new IGAControlPoint(32, 0.562500000000000,
				0.125000000000000, 4.500000000000000, 1.2);
		controlPointNet[33] = new IGAControlPoint(33, 0.562500000000000,
				0.250000000000000, 4.500000000000000, 1.0);
		controlPointNet[34] = new IGAControlPoint(34, 0.562500000000000,
				0.500000000000000, 4.500000000000000, 3.4);
		controlPointNet[35] = new IGAControlPoint(35, 0.562500000000000, 1.0,
				4.500000000000000, 5.6);

		// Fifth row
		controlPointNet[36] = new IGAControlPoint(36, 0.642857142857143, -1.0,
				4.500000000000000, 1.5);
		controlPointNet[37] = new IGAControlPoint(37, 0.642857142857143,
				-0.500000000000000, 4.500000000000000, 2.3);
		controlPointNet[38] = new IGAControlPoint(38, 0.642857142857143,
				-0.250000000000000, 4.500000000000000, 1.12);
		controlPointNet[39] = new IGAControlPoint(39, 0.642857142857143,
				-0.125000000000000, 4.500000000000000, 1.89);
		controlPointNet[40] = new IGAControlPoint(40, 0.642857142857143, 0.0,
				4.500000000000000, 2.21);
		controlPointNet[41] = new IGAControlPoint(41, 0.642857142857143,
				0.125000000000000, 4.500000000000000, 3.2);
		controlPointNet[42] = new IGAControlPoint(42, 0.642857142857143,
				0.250000000000000, 4.500000000000000, 4.3);
		controlPointNet[43] = new IGAControlPoint(43, 0.642857142857143,
				0.500000000000000, 4.500000000000000, 3.4);
		controlPointNet[44] = new IGAControlPoint(44, 0.642857142857143, 1.0,
				4.500000000000000, 2.6);

		// Sixth row
		controlPointNet[45] = new IGAControlPoint(45, 0.750000000000000, -1.0,
				4.500000000000000, 3.0);
		controlPointNet[46] = new IGAControlPoint(46, 0.750000000000000,
				-0.500000000000000, 4.500000000000000, 5.0);
		controlPointNet[47] = new IGAControlPoint(47, 0.750000000000000,
				-0.250000000000000, 4.500000000000000, 2.0);
		controlPointNet[48] = new IGAControlPoint(48, 0.750000000000000,
				-0.125000000000000, 4.500000000000000, 5.0);
		controlPointNet[49] = new IGAControlPoint(49, 0.750000000000000, 0.0,
				4.500000000000000, 6.0);
		controlPointNet[50] = new IGAControlPoint(50, 0.750000000000000,
				0.125000000000000, 4.500000000000000, 7.2);
		controlPointNet[51] = new IGAControlPoint(51, 0.750000000000000,
				0.250000000000000, 4.500000000000000, 2.0);
		controlPointNet[52] = new IGAControlPoint(52, 0.750000000000000,
				0.500000000000000, 4.500000000000000, 1.4);
		controlPointNet[53] = new IGAControlPoint(53, 0.750000000000000, 1.0,
				4.500000000000000, 3.6);

		// Seventh row
		controlPointNet[54] = new IGAControlPoint(54, 0.900000000000000, -1.0,
				4.500000000000000, 1.101);
		controlPointNet[55] = new IGAControlPoint(55, 0.900000000000000,
				-0.500000000000000, 4.500000000000000, 3.1232);
		controlPointNet[56] = new IGAControlPoint(56, 0.900000000000000,
				-0.250000000000000, 4.500000000000000, 1.001);
		controlPointNet[57] = new IGAControlPoint(57, 0.900000000000000,
				-0.125000000000000, 4.500000000000000, 4.12);
		controlPointNet[58] = new IGAControlPoint(58, 0.900000000000000, 0.0,
				4.500000000000000, 5.8);
		controlPointNet[59] = new IGAControlPoint(59, 0.900000000000000,
				0.125000000000000, 4.500000000000000, 1.2);
		controlPointNet[60] = new IGAControlPoint(60, 0.900000000000000,
				0.250000000000000, 4.500000000000000, 2.21);
		controlPointNet[61] = new IGAControlPoint(61, 0.900000000000000,
				0.500000000000000, 4.500000000000000, 1.12);
		controlPointNet[62] = new IGAControlPoint(62, 0.900000000000000, 1.0,
				4.500000000000000, 4.8);

		// Eighth row
		controlPointNet[63] = new IGAControlPoint(63, 1.125000000000000, -1.0,
				4.500000000000000, 9.0);
		controlPointNet[64] = new IGAControlPoint(64, 1.125000000000000,
				-0.500000000000000, 4.500000000000000, 3.3);
		controlPointNet[65] = new IGAControlPoint(65, 1.125000000000000,
				-0.250000000000000, 4.500000000000000, 1.1);
		controlPointNet[66] = new IGAControlPoint(66, 1.125000000000000,
				-0.125000000000000, 4.500000000000000, 4.4);
		controlPointNet[67] = new IGAControlPoint(67, 1.125000000000000, 0.0,
				4.500000000000000, 5.5);
		controlPointNet[68] = new IGAControlPoint(68, 1.125000000000000,
				0.125000000000000, 4.500000000000000, 1.0);
		controlPointNet[69] = new IGAControlPoint(69, 1.125000000000000,
				0.250000000000000, 4.500000000000000, 1.1);
		controlPointNet[70] = new IGAControlPoint(70, 1.125000000000000,
				0.500000000000000, 4.500000000000000, 3.0);
		controlPointNet[71] = new IGAControlPoint(71, 1.125000000000000, 1.0,
				4.500000000000000, 5.5);

		// Nineth row
		controlPointNet[72] = new IGAControlPoint(72, 1.500000000000000, -1.0,
				4.500000000000000, 0.23);
		controlPointNet[73] = new IGAControlPoint(73, 1.500000000000000,
				-0.500000000000000, 4.500000000000000, 3.032);
		controlPointNet[74] = new IGAControlPoint(74, 1.500000000000000,
				-0.250000000000000, 4.500000000000000, 1.01);
		controlPointNet[75] = new IGAControlPoint(75, 1.500000000000000,
				-0.125000000000000, 4.500000000000000, 4.40404);
		controlPointNet[76] = new IGAControlPoint(76, 1.500000000000000, 0.0,
				4.500000000000000, 5.0505);
		controlPointNet[77] = new IGAControlPoint(77, 1.500000000000000,
				0.125000000000000, 4.500000000000000, 1.01);
		controlPointNet[78] = new IGAControlPoint(78, 1.500000000000000,
				0.250000000000000, 4.500000000000000, 1.0);
		controlPointNet[79] = new IGAControlPoint(79, 1.500000000000000,
				0.500000000000000, 4.500000000000000, 3.0);
		controlPointNet[80] = new IGAControlPoint(80, 1.500000000000000, 1.0,
				4.500000000000000, 3.0);

		// 10th row
		controlPointNet[81] = new IGAControlPoint(81, 4.500000000000000, -1.0,
				4.500000000000000, 1.02);
		controlPointNet[82] = new IGAControlPoint(82, 4.500000000000000,
				-0.500000000000000, 4.500000000000000, 1.023);
		controlPointNet[83] = new IGAControlPoint(83, 4.500000000000000,
				-0.250000000000000, 4.500000000000000, 0.902);
		controlPointNet[84] = new IGAControlPoint(84, 4.500000000000000,
				-0.125000000000000, 4.500000000000000, 1.102);
		controlPointNet[85] = new IGAControlPoint(85, 4.500000000000000, 0.0,
				4.500000000000000, 3.2);
		controlPointNet[86] = new IGAControlPoint(86, 4.500000000000000,
				0.125000000000000, 4.500000000000000, 4.202);
		controlPointNet[87] = new IGAControlPoint(87, 4.500000000000000,
				0.250000000000000, 4.500000000000000, 1.01);
		controlPointNet[88] = new IGAControlPoint(88, 4.500000000000000,
				0.500000000000000, 4.500000000000000, 3.401);
		controlPointNet[89] = new IGAControlPoint(89, 4.500000000000000, 1.0,
				4.500000000000000, 5.9006);

//		for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
//			controlPointWeights[i] = controlPointNet[i].getW();

// Test just one object of the class (that works pretty also)
		theIGAPatchSurface = new IGAPatchSurface(id_basis, p, uNoKnots, uKnotVector, q, vNoKnots, vKnotVector,
				uNoControlPoints, vNoControlPoints, controlPointNet);

		// nurbsBasis1D->printPolynomialDegree();
		// nurbsBasis1D->printNoKnots();
		// nurbsBasis1D->printKnotVector();
		// nurbsBasis1D->printNoBasisFunctions();
		// nurbsBasis1D->printControlPointNet();
		// NurbsBasis1D nurbsBasis1DCopy = *nurbsBasis1D;
		// delete nurbsBasis1D;

		// Make test on memory leakage (that works pretty fine)
//		for (int j = 1; j < 1e9; j++) {
//
//			// Copy the knot vectors
//			double* uKnotVectorCopy = new double[uNoKnots];
//			for (int i = 0; i < uNoKnots; i++)
//				uKnotVectorCopy[i] = uKnotVector[i];
//			double* vKnotVectorCopy = new double[vNoKnots];
//			for (int i = 0; i < vNoKnots; i++)
//				vKnotVectorCopy[i] = vKnotVector[i];
//
//			// Copy the Control Point net and the Control Point weights
//			IGAControlPoint* controlPointNetCopy =
//					new IGAControlPoint[uNoControlPoints * vNoControlPoints];
//			for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
//				controlPointNetCopy[i] = controlPointNet[i];
//
//			// Create an object of the class IGAPatchSurface
//			IGAPatchSurface* IGAPatchSurfaceCopy = new IGAPatchSurface(id_patch,
//					id_basis, p, uNoKnots, uKnotVectorCopy, q, vNoKnots,
//					vKnotVectorCopy, uNoControlPoints, vNoControlPoints,
//					controlPointNetCopy);
//
//			// Free the memory on the heap from the pointer
//			delete IGAPatchSurfaceCopy;
//		}
	}

	void tearDown() {
		delete theIGAPatchSurface;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the constructor
	 ***********/
	void testConstructor() {
		CPPUNIT_ASSERT(theIGAPatchSurface->getIGABasis()->getId()==3);
		CPPUNIT_ASSERT(
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()==4);
		CPPUNIT_ASSERT(
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()==3);
		CPPUNIT_ASSERT(
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getNoKnots()==15);
		CPPUNIT_ASSERT(
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getNoKnots()==13);
		for (int i = 0; i < 5; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == -4);
		for (int i = 5; i < 7; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 5);
		for (int i = 7; i < 10; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 11);
		for (int i = 10; i < 15; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 21);
		for (int i = 0; i < 4; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == -11);
		for (int i = 4; i < 7; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == -0.1);
		for (int i = 7; i < 9; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == 0.0);
		for (int i = 9; i < 13; i++)
			CPPUNIT_ASSERT(
					theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == 1.0);
		CPPUNIT_ASSERT(theIGAPatchSurface->getUNoControlPoints()==10);
		CPPUNIT_ASSERT(theIGAPatchSurface->getVNoControlPoints()==9);
	}

	/***********************************************************************************************
	 * \brief Test case: Test the knot span
	 ***********/
	void testIGAPatchSurfaceKnotSpan() {

		// The parametric coordinates on the Spline parameter space
		double u = 13.3300000033001;
		double v = -5.00998989;
		int uCorrectknotSpan = 9;
		int vCorrectknotSpan = 3;

		// Check the find knot span function
		CPPUNIT_ASSERT(
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u)==uCorrectknotSpan);
		CPPUNIT_ASSERT(
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v)==vCorrectknotSpan);
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the local B-Spline basis functions
	 ***********/
	void testIGAPatchSurfaceBSplineBasisFunctions() {
		// Compute the non-zero basis functions at another parametric location
		double u = 17.002333009;
		double v = -0.000001234;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);
		int uNoLocalBasisFunctions =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1;
		int vNoLocalBasisFunctions =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
						+ 1;
		int noLocalBasisFunctions = uNoLocalBasisFunctions
				* vNoLocalBasisFunctions;

		// The polynomial degrees
		int p = 4;
		int q = 3;

		// Number of knots in both directions
		int uNoKnots = 15;
		int vNoKnots = 13;

		// The number of the Control Points at each parametric direction
		int uNoControlPoints = uNoKnots - p - 1;
		int vNoControlPoints = vNoKnots - q - 1;

		/*
		 * Make the patch's underlying basis a B-Spline by setting all the weights to 1. First copy everything to a new
		 * instance of the class IGAPatchSurface in order not to lose the stored Control Point information.
		 *
		 *              !!!!!!!!!!!!!!!!!!!!!!!!!
		 *              !!!!!   ACHTUNG     !!!!!
		 *              !!!!!!!!!!!!!!!!!!!!!!!!!
		 *
		 * By changing the Control Point weights to 1.0 we still call the computeLocalBasisFunctions from the NurbsBasis2D
		 * and not this function from BSplineBasis2D cause this is decided in the constructor of the instance. However it should be
		 * computed the B-Spline basis functions but in an inefficient way
		 *
		 */

		for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
			theIGAPatchSurface->getControlPointNet()[i]->setW(1.0);

		double* localBasisFunctions = new double[noLocalBasisFunctions];
		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctions(
				localBasisFunctions, u, uKnotSpan, v, vKnotSpan);

		/* cout << endl;
		 cout << endl;
		 cout << "The non-zero B-Spline basis functions at (u,v) = ( " << u << " , " << v << " ):";
		 cout << endl;
		 for (int i = 0; i < noLocalBasisFunctions; i++) {
		 cout << localBasisFunctions[i] << " ";
		 }
		 cout << endl;
		 cout << endl;*/

		// Value provided by MATLAB
		double CorrectlocalBasisFunctions[] = { 0.000000000000000,
				0.000000000000000, 0.000000000000001, 0.000000000000001,
				0.000000000000000, 0.000000000007292, 0.000000000074448,
				0.000000000157816, 0.000000000157969, 0.000000000059296,
				0.014511603260064, 0.148153855088449, 0.314060685220477,
				0.314366157723487, 0.118001972754625, 0.001451101233533,
				0.014814782213839, 0.031404789640423, 0.031435335646814,
				0.011799716761467 };
		for (int i = 0; i < noLocalBasisFunctions; i++) {
			cout << i << " :  " << localBasisFunctions[i] << "    -    "
					<< CorrectlocalBasisFunctions[i] << "    =     "
					<< localBasisFunctions[i] - CorrectlocalBasisFunctions[i]
					<< endl;
			CPPUNIT_ASSERT(
					fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
		}

		// Clear the heap from the pointer
		delete[] localBasisFunctions;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the local NURBS basis functions
	 ***********/
	void testIGAPatchSurfaceNurbsBasisFunctions() {
		// Compute the non-zero basis functions at another parametric location
		double u = 17.002333009;
		double v = -0.000001234;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);
		int uNoLocalBasisFunctions =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1;
		int vNoLocalBasisFunctions =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
						+ 1;
		int noLocalBasisFunctions = uNoLocalBasisFunctions
				* vNoLocalBasisFunctions;

		// The polynomial degrees
		int p = 4;
		int q = 3;

		// Number of knots in both directions
		int uNoKnots = 15;
		int vNoKnots = 13;

		// The number of the Control Points at each parametric direction
		int uNoControlPoints = uNoKnots - p - 1;
		int vNoControlPoints = vNoKnots - q - 1;

		double* localBasisFunctions = new double[noLocalBasisFunctions];
		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctions(
				localBasisFunctions, u, uKnotSpan, v, vKnotSpan);

		/*cout << endl;
		 cout << endl;
		 cout << "The non-zero 2D NURBS basis functions at (u,v) = ( " << u << " , " << v << " ):";
		 cout << endl;
		 for (int i = 0; i < noLocalBasisFunctions; i++) {
		 cout << localBasisFunctions[i] << " ";
		 }
		 cout << endl;
		 cout << endl;*/

		// Value provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 4.495585113804783e-17, 1.229382471599235e-16, 3.995988326933906e-16,
                                                1.947765238090499e-16, 1.233769054530474e-16, 5.854909519645797e-12,
                                                9.962457868877927e-11, 8.447478718422842e-11, 2.113923798912647e-10,
                                                9.521909622754816e-11, 1.213001775954701e-02, 3.969111309430823e-02,
                                                3.463036222447142e-01, 4.879889914307526e-01, 3.789804072907874e-02,
                                                4.272053085451214e-04, 1.744593267633994e-02, 4.622796996550126e-02,
                                                8.413260691500543e-03, 3.473845603145905e-03};
        for (int i = 0; i < noLocalBasisFunctions; i++) {
            CPPUNIT_ASSERT(
                        fabs(localBasisFunctions[i] - CorrectlocalBasisFunctions[i])<=Tol);
        }

		// Clear the heap from the pointer
		delete[] localBasisFunctions;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the local basis functions and their derivatives (Test 1)
	 ***********/
	void testIGAPatchSurfaceBasisFunctionsAndDerivativesTest1() {
		// Compute the non-zero basis functions and their derivatives at another parametric location
		double u = 5.0;
		double v = -0.1;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// The number of local basis functions
		int noBasisFunctions =
				(theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1)
						* (theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
								+ 1);

		// Return the functional values R, dR/du, dR/dv
		int derivDegree = 1;
		double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
				* (derivDegree + 2) * noBasisFunctions / 2];

		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
				localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v,
				vKnotSpan);

		// Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 1.012658227848101e-01, 2.531645569620253e-01, 6.455696202531646e-01, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              -5.217112642204775e+00, -1.304278160551194e+01, -3.325909309405544e+01, -0.000000000000000e+00,
                                                              -0.000000000000000e+00, 1.088607594936709e+01, 2.582278481012658e+01, 1.481012658227848e+01,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
                                                              -8.688066371129982e-02, -7.655468318823551e-02, 1.634353468995353e-01, 0.000000000000000e+00,
                                                              0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00, 0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
                                                              0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, -0.000000000000000e+00,
                                                              -0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00};

		// The tolerance value is different for the functional and the derivative values due to high gradients which appear
		int indexBasisDeriv = 0;
		int counter = 0;
		for (int i = 0; i <= derivDegree; i++)
			for (int j = 0; j <= derivDegree - i; j++)
				for (int k = 0; k < noBasisFunctions; k++) {
					// Get the index of the basis function
					indexBasisDeriv =
							theIGAPatchSurface->getIGABasis()->indexDerivativeBasisFunction(
									derivDegree, i, j, k);

					// Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(localBasisFunctionsAndDerivatives[indexBasisDeriv] - CorrectlocalBasisFunctionsAndDerivatives[counter])<=Tol*10);
					counter++;
                }

		// Clear the heap from the pointer
		delete[] localBasisFunctionsAndDerivatives;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the local basis functions and their derivatives (Test 2)
	 ***********/
	void testIGAPatchSurfaceBasisFunctionsAndDerivativesTest2() {
		// Compute the non-zero basis functions and their derivatives at another parametric location
		double u = 10.1100000099;
		double v = -.44449811111;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// The number of local basis functions
		int noBasisFunctions =
				(theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1)
						* (theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
								+ 1);

		// Return the functional values R, dR/du, dR/dv, d^2R/du^2
		int derivDegree = 2;

		// Compute the local basis functions and their derivatives
		double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
				* (derivDegree + 2) * noBasisFunctions / 2];
		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
				localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v,
				vKnotSpan);

		// Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 5.665652111327743e-10, 1.369178057561520e-07, 3.752321486533538e-06, 6.103217419530942e-06, 1.443061319105917e-06,
                                                              2.546666413064336e-07, 1.639271770281257e-05, 1.379663989838727e-05, 4.675099871563213e-04, 4.510004760008026e-04,
                                                              5.106287021474073e-06, 5.784376191418766e-04, 2.113661204680182e-03, 2.005445995406023e-02, 2.032170442428270e-02,
                                                              1.955722707982197e-05, 9.846368100513380e-04, 7.339814453205543e-02, 8.192977284885009e-01, 6.226617288113014e-02,
                                                              -5.018423789140440e-09, -1.212766968500929e-06, -3.323666727590472e-05, -5.406002854858067e-05, -1.278210011961427e-05,
                                                              -1.492376309940641e-06, -9.606324342176029e-05, -8.084992378864370e-05, -2.739663215855675e-03, -2.642915549138109e-03,
                                                              -1.461728850095665e-05, -1.655839071177488e-03, -6.050579509565212e-03, -5.740801988736555e-02, -5.817303554455428e-02,
                                                              2.638362534965420e-06, 1.328321678520368e-04, 9.901757231685975e-03, 1.105271428819554e-01, 8.400001546985621e-03,
                                                              3.013060805133030e-08, 7.281450854065949e-06, 1.995528948331980e-04, 3.245763211480399e-04, 7.674370777740269e-05,
                                                              4.452006301068493e-06, 2.865726038171726e-04, 2.411887455934393e-04, 8.172870219495891e-03, 7.884255867357446e-03,
                                                              -1.264502157999187e-06, -1.432421668811061e-04, -5.234192953426450e-04, -4.966212784900777e-03, -5.032392223676647e-03,
                                                              -1.336115307519445e-07, -6.726865260024795e-06, -5.014431956664508e-04, -5.597297776323400e-03, -4.253915260459686e-04,
                                                              -2.538516510637189e-09, -4.418005668735724e-07, -7.250303219664627e-06, 6.176860712049581e-07, 1.149576884081068e-06,
                                                              -1.141043362638443e-06, -5.289532602209136e-05, -2.665811632473955e-05, 4.731511059913412e-05, 3.592776793721455e-04,
                                                              -2.287890900704561e-05, -1.866477969220580e-03, -4.084054282814173e-03, 2.029644321619571e-03, 1.618875188599570e-02,
                                                              -8.762688366471939e-05, -3.177184285439229e-03, -1.418212180187655e-01, 8.291836260621074e-02, 4.960271061017116e-02,
                                                              2.250710674164042e-08, 3.918590706099079e-06, 6.436534581969372e-05, -5.235618994841844e-06, -1.012681491942663e-05,
                                                              6.696478979349807e-06, 3.106056382867790e-04, 1.567522987001755e-04, -2.592238877856792e-04, -2.087998384081443e-03,
                                                              6.569043295342129e-05, 5.365321815980380e-03, 1.177263656854971e-02, -5.035871024040827e-03,-4.555750209782555e-02,
                                                              -1.106627591683007e-05, -3.906053579475987e-04, -1.629882279020058e-02, 4.281502012367986e-02, 9.095422450950559e-03,
                                                              8.602822884368602e-09, 9.264461788580807e-07, 4.967539832674171e-06, -1.193449805947254e-06, 9.239174725680466e-07,
                                                              3.866901755820748e-06, 1.109203481095631e-04, 1.826478847781635e-05, -9.141894595867033e-05, 2.887522618732625e-04,
                                                              7.753473382995434e-05, 3.913963702545520e-03, 2.798186739784572e-03, -3.921536739615213e-03, 1.301093553085535e-02,
                                                              2.969602746005120e-04, 6.662486337676968e-03, 9.716870154006790e-02, -1.602090582499311e-01, 3.986580771942673e-02};

		// The tolerance value is different for the functional and the derivative values due to high gradients which appear
		int indexBasisDeriv = 0;
		int counter = 0;
		for (int i = 0; i <= derivDegree; i++)
			for (int j = 0; j <= derivDegree - i; j++)
				for (int k = 0; k < noBasisFunctions; k++) {
					// Get the index of the basis function
					indexBasisDeriv =
							theIGAPatchSurface->getIGABasis()->indexDerivativeBasisFunction(
									derivDegree, i, j, k);

					// Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(localBasisFunctionsAndDerivatives[indexBasisDeriv]-CorrectlocalBasisFunctionsAndDerivatives[counter])<=relTol);
					counter++;
				}

		// Clear the heap from the pointer
		delete[] localBasisFunctionsAndDerivatives;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the Coordinates of a point on the IGA 2D patch
	 ***********/
	void testIGAPatchSurfacePointOnSurface() {

		// theIGAPatchSurface->printControlPointNet();

		// The parameters on the IGA surface and their knot span indices
		double u = 8.000000000021;
		double v = -2.11111111231;
		// double u =
		//		theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[0];
		// double v =
		//		theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[0];

		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// Initialize the Cartesian Coordinates
		double cartesianCoordinatesPointOnPatch[3];

		// Compute the Cartesian Coordinates of the point with the above surface parameters
		theIGAPatchSurface->computeCartesianCoordinates(
				cartesianCoordinatesPointOnPatch, u, uKnotSpan, v, vKnotSpan);

		// Values provided by MATLAB
        double CorrectcartesianCoordinatesPointOnPatch[] = { 5.389614693204954e-01,
                                                             5.700831281051164e-01,
                                                             4.500000000000000e+00};

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(
                        fabs(cartesianCoordinatesPointOnPatch[i]-CorrectcartesianCoordinatesPointOnPatch[i])<=Tol);
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the Coordinates of a point on the IGA 2D patch by pre-computing the basis functions
	 ***********/
	void testIGAPatchSurfacePointOnSurfaceMethod2() {

		// The parameters on the IGA surface and their knot span indices
		double u = 8.000000000021;
		double v = -2.11111111231;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// Initialize the Cartesian Coordinates
		double cartesianCoordinatesPointOnPatch[3];

		// Compute the local basis functions
		int pDegree =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
		int qDegree =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
		int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
		double* localBasisFunctions = new double[noLocalBasisFunctions];
		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctions(
				localBasisFunctions, u, uKnotSpan, v, vKnotSpan);

		// Compute the Cartesian Coordinates of the point with the above surface parameters
		theIGAPatchSurface->computeCartesianCoordinates(
				cartesianCoordinatesPointOnPatch, localBasisFunctions,
				uKnotSpan, vKnotSpan);

		// Values provided by MATLAB
        double CorrectcartesianCoordinatesPointOnPatch[] = {5.389614693204954e-01,
                                                            5.700831281051164e-01,
                                                            4.500000000000000e+00};

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(
                        fabs(cartesianCoordinatesPointOnPatch[i] - CorrectcartesianCoordinatesPointOnPatch[i])<=Tol);

		// Free the memory from the heap
		delete[] localBasisFunctions;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the base vectors for the surface patch
	 ***********/
	void testIGAPatchSurfaceBaseVectors() {

		// The parameters on the IGA surface and their knot span indices
		double u = 2.222122;
		double v = -3.3339333;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// Modify the Control Point net for getting some curvature into the structure
		theIGAPatchSurface->getControlPointNet()[theIGAPatchSurface->getVNoControlPoints()
				+ 4]->setZ(102.78054 * 4.500000000000000);
		theIGAPatchSurface->getControlPointNet()[4
				* theIGAPatchSurface->getVNoControlPoints()]->setZ(
				102.78054 * 4.500000000000000 / 1.5672);

		// Get the polynomial degree of the basis in each direction
		int pDegree =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
		int qDegree =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
		int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

		// Derivative order of the basis functions needed for the computation of the base vectors
		int derivDegree = 1;

		// Compute the local basis functions and their derivatives
		double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
				* (derivDegree + 2) * noLocalBasisFunctions / 2];
		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
				localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v,
				vKnotSpan);

		// Compute the base vectors at the given surface parameters
		double baseVectors[6];
		theIGAPatchSurface->computeBaseVectors(baseVectors,
				localBasisFunctionsAndDerivatives, uKnotSpan, vKnotSpan);

		// Values provided by MATLAB
        double correctBaseVectors[] = {-3.438137561383751e-03,
                                       3.670146188906259e-02,
                                       1.424647671570531e+01,
                                       1.849566126177724e-02,
                                       2.642578782315075e-02,
                                       -1.876668817294945e+01};

        for (int i = 0; i < 6; i++)
            CPPUNIT_ASSERT( fabs(baseVectors[i]-correctBaseVectors[i]) <= TolDeriv);

		// Free the memory from the heap
		delete[] localBasisFunctionsAndDerivatives;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the base vectors and their derivatives for the surface patch
	 ***********/
	void testIGAPatchSurfaceBaseVectorsAndDerivatives() {

		// Testing if the implemented index works correctly

		// The parameters on the IGA surface and their knot span indices
		double u = 9.0012100121;
		double v = -1.00010001;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// Modify the Control Point net for getting some curvature into the structure
		theIGAPatchSurface->getControlPointNet()[theIGAPatchSurface->getVNoControlPoints()
				+ 4]->setZ(102.78054 * 4.500000000000000);
		theIGAPatchSurface->getControlPointNet()[4
				* theIGAPatchSurface->getVNoControlPoints()]->setZ(
				102.78054 * 4.500000000000000 / 1.5672);

		// Get the polynomial degrees of the patch
		int pDegree =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
		int qDegree =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

		// Decide on absolute degree of the partial derivatives for the base vectors
		int derivDegreeBaseVec = 1;

		// Initialize the index of the base vector
		int indexBaseVct = 0;

		// Indices related to the base vectors
		int indexUBaseVec = 0;
		int indexVBaseVec = 1;

		// Compute the local basis functions and their derivatives
		int noBasisFcts = (pDegree + 1) * (qDegree + 1);
		int derivDegree = derivDegreeBaseVec + 1;
		double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
				* (derivDegree + 2) * noBasisFcts / 2];
		theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
				localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v,
				vKnotSpan);

		// Compute the base vectors and their derivatives
		int noBaseVec = 2;
		int noCoord = 3;
		double* baseVctsAndDerivs = new double[(derivDegreeBaseVec + 1)
				* (derivDegreeBaseVec + 2) * noCoord * noBaseVec / 2];
		theIGAPatchSurface->computeBaseVectorsAndDerivatives(baseVctsAndDerivs,
				localBasisFunctionsAndDerivatives, derivDegreeBaseVec,
				uKnotSpan, vKnotSpan);

		// Tests on the pointer baseVctsAndDerivs

		// Analytical values provided by MATLAB
        double correctBaseVctGuAndDerivs[] = { 2.216785405425909e-03, 9.242381431419883e-02, 5.213713223726086e+00,
                                               -6.291438422062177e-04, -1.907329196461980e-02, 2.570460922816262e+00,
                                               -1.132430305658589e-03, -6.826525391671297e-02, 4.839841982005542e+00};
        double correctBaseVctGvAndDerivs[] = { 9.073264459056658e-03, 4.006608792432891e-02, -3.563132653244571e-01,
                                               -3.137323441257892e-03, -8.431822133857268e-03, 1.444111659374656e+00,
                                               -6.291438422062177e-04, -1.907329196461980e-02, 2.570460922816262e+00};

		// Testing the base vector Gu and its derivatives
		int counterBaseVec = 0;
		for (int i = 0; i <= derivDegreeBaseVec; i++) {
			for (int j = 0; j <= derivDegreeBaseVec - i; j++) {
				for (int k = 0; k < noCoord; k++) {
					// Find the index of the derivatives to the base vector
					indexBaseVct =
							theIGAPatchSurface->indexDerivativeBaseVector(
									derivDegreeBaseVec, i, j, k, indexUBaseVec);

					// Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(baseVctsAndDerivs[indexBaseVct] - correctBaseVctGuAndDerivs[counterBaseVec]) <= relTol);

					// Update counter
					counterBaseVec++;
				}
			}
		}

		// Testing the base vector Gv and its derivatives
		counterBaseVec = 0;
		for (int i = 0; i <= derivDegreeBaseVec; i++) {
			for (int j = 0; j <= derivDegreeBaseVec - i; j++) {
				for (int k = 0; k < noCoord; k++) {
					// Find the index of the derivatives to the base vector
					indexBaseVct =
							theIGAPatchSurface->indexDerivativeBaseVector(
									derivDegreeBaseVec, i, j, k, indexVBaseVec);

					// Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                                fabs(baseVctsAndDerivs[indexBaseVct] - correctBaseVctGvAndDerivs[counterBaseVec])<=Tol);

					// Update counter
					counterBaseVec++;
				}
			}
        }

		// Free the memory from the heap
		delete[] localBasisFunctionsAndDerivatives;
		delete[] baseVctsAndDerivs;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the projection of an arbitrary point on the 2D IGA patch
	 ***********/
	void testProjectionOnIGAPatch() {


		// Modify the Control Point net for getting some curvature into the structure
		theIGAPatchSurface->getControlPointNet()[theIGAPatchSurface->getVNoControlPoints()
				+ 4]->setZ(102.78054 * 4.500000000000000);
		theIGAPatchSurface->getControlPointNet()[4
				* theIGAPatchSurface->getVNoControlPoints()]->setZ(
				102.78054 * 4.500000000000000 / 1.5672);

		// Initial guesses for the Newton-Rapson iteration
		double u = 10.0002;
		double v = -6.2365;

		// The vertex to be projected onto the NURBS patch
        double vertex[] = { .5, -1.72, 1.8 };

		// Flag on the convergence of the Newton-Rapson iterations
		bool flag = 1;

		// Compute the orthogonal projection of the point on the NURBS patch
        flag = theIGAPatchSurface->computePointProjectionOnPatch(u, v, vertex);

		// Compare the values with the ones from MATLAB

		// On the return flag
        CPPUNIT_ASSERT(flag == 1);

		// On the Cartesian components of the orthogonal projection
        double orthogonalProjection[] = {1.337577535013243e-01,
                                         6.359047538272620e-01,
                                         4.525192863798665e+00};
        
	for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(fabs(orthogonalProjection[i] - vertex[i]) <= relTol);

	// On the surface parametric values of the orthogonal projection
        CPPUNIT_ASSERT(fabs(u - 1.638675138407125e+01) <= TolDeriv);
        CPPUNIT_ASSERT(fabs(v + 9.976275205718462e+00) <= TolDeriv);
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the NURBS basis functions for memory leakage
	 ***********/
	void testIGAPatchSurfaceBasisFunctions4Leakage() {

		// Initialize variables
		double u = 17.002333009;
		double v = -0.000001234;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);
		int uNoLocalBasisFunctions =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1;
		int vNoLocalBasisFunctions =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
						+ 1;
		int noLocalBasisFunctions = uNoLocalBasisFunctions
				* vNoLocalBasisFunctions;

		// The polynomial degrees
		int p = 4;
		int q = 3;

		// Number of knots in both directions
		int uNoKnots = 15;
		int vNoKnots = 13;

		// The number of the Control Points at each parametric direction
		int uNoControlPoints = uNoKnots - p - 1;
		int vNoControlPoints = vNoKnots - q - 1;

		for (int i = 1; i < 1e9; i++) {
			// Compute the non-zero basis functions at another parametric location
			double* localBasisFunctions = new double[noLocalBasisFunctions];
			theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctions(
					localBasisFunctions, u, uKnotSpan, v, vKnotSpan);

			// Free the memory from the heap
			delete[] localBasisFunctions;
		}
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (Test 1)
	 ***********/
	void testIGAPatchSurfaceBasisFunctionsAndDerivativesTest14Leakage() {

		// Initialize variables
		double u = 5.0;
		double v = -0.1;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// The number of local basis functions
		int noBasisFunctions =
				(theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1)
						* (theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
								+ 1);

		// Return the functional values R, dR/du, dR/dv
		int derivDegree = 1;

		// Number of iterations
		int noIterations = 1e9;

		for (int i = 0; i < noIterations; i++) {
			// Compute the non-zero basis functions and their derivatives at another parametric location
			double* localBasisFunctionsAndDerivatives = new double[(derivDegree
					+ 1) * (derivDegree + 1) * noBasisFunctions];

			theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
					localBasisFunctionsAndDerivatives, derivDegree, u,
					uKnotSpan, v, vKnotSpan);

			// Clear the heap from the pointer
			delete[] localBasisFunctionsAndDerivatives;
		}
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (Test 2)
	 ***********/
	void testIGAPatchSurfaceBasisFunctionsAndDerivativesTest24Leakage() {

		// Initialize variables
		double u = 10.1100000099;
		double v = -.44449811111;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// The number of local basis functions
		int noBasisFunctions =
				(theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
						+ 1)
						* (theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
								+ 1);

		// Return the functional values R, dR/du, dR/dv, d^2R/du^2
		int derivDegree = 2;

		// Number of iterations
		int noIterations = 1e9;

		for (int i = 0; i < noIterations; i++) {
			// Compute the non-zero basis functions and their derivatives at another parametric location
			double* localBasisFunctionsAndDerivatives = new double[(derivDegree
					+ 1) * (derivDegree + 1) * noBasisFunctions];

			theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
					localBasisFunctionsAndDerivatives, derivDegree, u,
					uKnotSpan, v, vKnotSpan);

			// Clear the heap from the pointer
			delete[] localBasisFunctionsAndDerivatives;
		}
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the Coordinates of a point on the IGA patch for leakage
	 ***********/
	void testIGAPatchSurfacePointOnSurface4Leakage() {

		// The parameters on the IGA surface and their knot span indices
		double u = 8.000000000021;
		double v = -2.11111111231;
		int uKnotSpan =
				theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
						u);
		int vKnotSpan =
				theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
						v);

		// Initialize the Cartesian Coordinates
		double cartesianCoordinatesPointOnPatch[3];

		for (int i = 0; i < 1e9; i++) {
			// Compute the Cartesian Coordinates of the point with the above surface parameters
			theIGAPatchSurface->computeCartesianCoordinates(
					cartesianCoordinatesPointOnPatch, u, uKnotSpan, v,
					vKnotSpan);
		}
	}

	/***********************************************************************************************
	 * \brief Test case: Test the computation of the projection of a point on the IGA patch for leakage
	 ***********/
	void testIGAPatchSurfacePointProjection4Leakage() {
		// Modify the Control Point net for getting some curvature into the structure
		theIGAPatchSurface->getControlPointNet()[theIGAPatchSurface->getVNoControlPoints()
				+ 4]->setZ(102.78054 * 4.500000000000000);
		theIGAPatchSurface->getControlPointNet()[4
				* theIGAPatchSurface->getVNoControlPoints()]->setZ(
				102.78054 * 4.500000000000000 / 1.5672);

		// Initial guesses for the Newton-Rapson iteration
		double u = 10.0002;
		double v = -6.2365;

		// The vertex to be projected onto the NURBS patch
		double vertex[] = { .5, -1.72, 3.8 };

		// Flag on the convergence of the Newton-Rapson iterations
		bool flag = 1;

		// Compute the orthogonal projection of the point on the NURBS patch iteratively
		int noIterations = 1e9;
		for (int i = 0; i < noIterations; i++) {
			flag = theIGAPatchSurface->computePointProjectionOnPatch(u, v,
					vertex);
		}
	}

	void testIGAPatchSurfaceFindNearestKnotIntersection() {
		const int numNodes = 77;
		double nodes[numNodes][3] = {
				{ 0, -1.31111111111111, 9.80000000000000 }, { 0.777380192706971,
						-1.31111111111111, 9.80000000000000 }, {
						1.68974345802307, -1.31111111111111, 9.80000000000000 },
				{ 2.71945888878606, -1.31111111111111, 9.80000000000000 }, {
						3.83641272531922, -1.31111111111111, 9.80000000000000 },
				{ 5, -1.31111111111111, 9.80000000000000 }, { 6.16358727468078,
						-1.31111111111111, 9.80000000000000 }, {
						7.28054111121394, -1.31111111111111, 9.80000000000000 },
				{ 8.31025654197693, -1.31111111111111, 9.80000000000000 }, {
						9.22261980729303, -1.31111111111111, 9.80000000000000 },
				{ 10.0000000000000, -1.31111111111111, 9.80000000000000 }, { 0,
						-1.45432098765432, 8.66111111111111 },
				{ 0.777380192706971, -1.45432098765432, 8.66111111111111 }, {
						1.68974345802307, -1.45432098765432, 8.66111111111111 },
				{ 2.71945888878606, -1.45432098765432, 8.66111111111111 }, {
						3.83641272531922, -1.45432098765432, 8.66111111111111 },
				{ 5.00000000000000, -1.45432098765432, 8.66111111111111 }, {
						6.16358727468078, -1.45432098765432, 8.66111111111111 },
				{ 7.28054111121394, -1.45432098765432, 8.66111111111111 }, {
						8.31025654197693, -1.45432098765432, 8.66111111111111 },
				{ 9.22261980729303, -1.45432098765432, 8.66111111111111 }, {
						10.0000000000000, -1.45432098765432, 8.66111111111112 },
				{ 0, -0.901234567901235, 7.84444444444444 },
				{ 0.777380192706971, -0.901234567901235, 7.84444444444444 },
				{ 1.68974345802307, -0.901234567901235, 7.84444444444444 },
				{ 2.71945888878606, -0.901234567901235, 7.84444444444444 },
				{ 3.83641272531922, -0.901234567901235, 7.84444444444444 },
				{ 5.00000000000000, -0.901234567901235, 7.84444444444445 },
				{ 6.16358727468078, -0.901234567901235, 7.84444444444444 },
				{ 7.28054111121394, -0.901234567901235, 7.84444444444445 },
				{ 8.31025654197693, -0.901234567901234, 7.84444444444444 },
				{ 9.22261980729303, -0.901234567901235, 7.84444444444444 }, {
						10, -0.901234567901235, 7.84444444444445 }, { 0,
						-4.16333634234434e-17, 7.55000000000000 }, {
						0.777380192706971, 5.46437894932694e-17,
						7.55000000000000 }, { 1.68974345802307,
						-8.67361737988404e-19, 7.55000000000000 },
				{ 2.71945888878606, 3.81639164714898e-17, 7.55000000000000 }, {
						3.83641272531922, -3.98986399474666e-17,
						7.55000000000000 }, { 5, -6.93889390390723e-18,
						7.55000000000000 }, { 6.16358727468078,
						-6.93889390390723e-18, 7.55000000000000 },
				{ 7.28054111121394, 6.93889390390723e-18, 7.55000000000000 }, {
						8.31025654197693, -8.32667268468867e-17,
						7.55000000000000 }, { 9.22261980729303, 0,
						7.55000000000000 }, { 10.0000000000000,
						6.93889390390723e-17, 7.55000000000000 }, { 0,
						0.901234567901235, 7.84444444444445 },
				{ 0.777380192706971, 0.901234567901235, 7.84444444444445 }, {
						1.68974345802307, 0.901234567901234, 7.84444444444444 },
				{ 2.71945888878606, 0.901234567901235, 7.84444444444444 }, {
						3.83641272531922, 0.901234567901235, 7.84444444444445 },
				{ 5, 0.901234567901235, 7.84444444444445 }, { 6.16358727468078,
						0.901234567901234, 7.84444444444444 }, {
						7.28054111121394, 0.901234567901235, 7.84444444444444 },
				{ 8.31025654197693, 0.901234567901235, 7.84444444444444 }, {
						9.22261980729303, 0.901234567901235, 7.84444444444445 },
				{ 10.0000000000000, 0.901234567901234, 7.84444444444444 }, { 0,
						1.45432098765432, 8.66111111111111 }, {
						0.777380192706971, 1.45432098765432, 8.66111111111111 },
				{ 1.68974345802307, 1.45432098765432, 8.66111111111111 }, {
						2.71945888878606, 1.45432098765432, 8.66111111111111 },
				{ 3.83641272531922, 1.45432098765432, 8.66111111111111 }, { 5,
						1.45432098765432, 8.66111111111111 }, {
						6.16358727468078, 1.45432098765432, 8.66111111111111 },
				{ 7.28054111121394, 1.45432098765432, 8.66111111111111 }, {
						8.31025654197693, 1.45432098765432, 8.66111111111111 },
				{ 9.22261980729303, 1.45432098765432, 8.66111111111111 }, {
						10.0000000000000, 1.45432098765432, 8.66111111111111 },
				{ 0, 1.31111111111111, 9.80000000000000 }, { 0.777380192706971,
						1.31111111111111, 9.80000000000000 }, {
						1.68974345802307, 1.31111111111111, 9.80000000000000 },
				{ 2.71945888878606, 1.31111111111111, 9.80000000000000 }, {
						3.83641272531922, 1.31111111111111, 9.80000000000000 },
				{ 5, 1.31111111111111, 9.80000000000000 }, { 6.16358727468078,
						1.31111111111111, 9.80000000000000 }, {
						7.28054111121394, 1.31111111111111, 9.80000000000000 },
				{ 8.31025654197693, 1.31111111111111, 9.80000000000000 }, {
						9.22261980729303, 1.31111111111111, 9.80000000000000 },
				{ 10, 1.31111111111111, 9.80000000000000 } };

		double correctU[numNodes] = { 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1, 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1, 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1, 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1, 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1, 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1, 0, 0.0833333333333333, 0.166666666666667,
				0.333333333333333, 0.416666666666667, 0.500000000000000,
				0.583333333333333, 0.666666666666667, 0.833333333333333,
				0.916666666666667, 1 };
		double correctV[numNodes] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0.200000000000000, 0.200000000000000, 0.200000000000000,
				0.200000000000000, 0.200000000000000, 0.200000000000000,
				0.200000000000000, 0.200000000000000, 0.200000000000000,
				0.200000000000000, 0.200000000000000, 0.400000000000000,
				0.400000000000000, 0.400000000000000, 0.400000000000000,
				0.400000000000000, 0.400000000000000, 0.400000000000000,
				0.400000000000000, 0.400000000000000, 0.400000000000000,
				0.400000000000000, 0.600000000000000, 0.600000000000000,
				0.600000000000000, 0.400000000000000, 0.600000000000000,
				0.600000000000000, 0.600000000000000, 0.600000000000000,
				0.400000000000000, 0.400000000000000, 0.600000000000000,
				0.600000000000000, 0.600000000000000, 0.600000000000000,
				0.600000000000000, 0.600000000000000, 0.600000000000000,
				0.600000000000000, 0.600000000000000, 0.600000000000000,
				0.600000000000000, 0.600000000000000, 0.800000000000000,
				0.800000000000000, 0.800000000000000, 0.800000000000000,
				0.800000000000000, 0.800000000000000, 0.800000000000000,
				0.800000000000000, 0.800000000000000, 0.800000000000000,
				0.800000000000000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		double u, v;
		for (int i = 0; i < numNodes; i++) {
			theIGAPatchSurface->findInitialGuess4PointProjection(u, v, nodes[i]);
			cout << i << ": " << "(" << u << "," << v << "), (" << correctU[i]
					<< "," << correctV[i] << ")" << endl;
			CPPUNIT_ASSERT(fabs(u - correctU[i])<=Tol);
			CPPUNIT_ASSERT(fabs(v - correctV[i])<=Tol);
		}

	}

// Make the tests
CPPUNIT_TEST_SUITE(TestIGAPatchSurface);
	CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testIGAPatchSurfaceKnotSpan);
	CPPUNIT_TEST(testIGAPatchSurfaceNurbsBasisFunctions);
	CPPUNIT_TEST(testIGAPatchSurfaceBasisFunctionsAndDerivativesTest1);
	CPPUNIT_TEST(testIGAPatchSurfaceBasisFunctionsAndDerivativesTest2);
	CPPUNIT_TEST(testIGAPatchSurfacePointOnSurface);
	CPPUNIT_TEST(testIGAPatchSurfacePointOnSurfaceMethod2);
	CPPUNIT_TEST(testIGAPatchSurfaceBaseVectors);
	CPPUNIT_TEST(testIGAPatchSurfaceBaseVectorsAndDerivatives);
	CPPUNIT_TEST(testProjectionOnIGAPatch);
	//CPPUNIT_TEST(testIGAPatchSurfaceFindNearestKnotIntersection);  can't find the correct solution from MATLAB yet, do it later

// Make the tests for leakage
	// CPPUNIT_TEST(testIGAPatchSurfaceCopyConstructor4Leakage);
	// CPPUNIT_TEST(testIGAPatchSurfaceBasisFunctions4Leakage);
	// CPPUNIT_TEST(testIGAPatchSurfaceBasisFunctionsAndDerivativesTest14Leakage);
	// CPPUNIT_TEST(testIGAPatchSurfaceBasisFunctionsAndDerivativesTest24Leakage);
	// CPPUNIT_TEST(testIGAPatchSurfacePointOnSurface4Leakage);
	// CPPUNIT_TEST(testIGAPatchSurfacePointProjection4Leakage);

	CPPUNIT_TEST_SUITE_END()
	;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatchSurface);
