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

class TestIGAPatchCurve: public CppUnit::TestFixture {

private:
    IGAPatchCurve* theIGAPatchCurve;
	double Tol;
    double TolRel;
    double TolRel10000;
    double TolRel1000000;

public:
	void setUp() {
        // Define the tolerances of the test case
        Tol = 1e-15;
        TolRel = Tol*10;
        TolRel10000 = Tol*1e4;
        TolRel1000000 = Tol*1e6;

        // Polynomial order
        int p = 3;

        // number of knots
        int noKnots = 11;

        // Knot vector
        double* knotVector = new double[noKnots];
        knotVector[0] = 0.0;
        knotVector[1] = 0.0;
        knotVector[2] = 0.0;
        knotVector[3] = 0.0;
        knotVector[4] = 0.30971;
        knotVector[5] = 0.46851;
        knotVector[6] = 0.81417;
        knotVector[7] = 1.0;
        knotVector[8] = 1.0;
        knotVector[9] = 1.0;
        knotVector[10] = 1.0;

        // Number of Control Points
        int noCPs = noKnots - p - 1;

        // Control Point coordinates and weights
        double* controlPointCoords = new double[4*noCPs];

        // 1st
        controlPointCoords[0*4 + 0] = -6.3606;
        controlPointCoords[0*4 + 1] = -0.0055991;
        controlPointCoords[0*4 + 2] = 0.0;
        controlPointCoords[0*4 + 3] = 1.0;

        // 2nd
        controlPointCoords[1*4 + 0] = -5.405;
        controlPointCoords[1*4 + 1] = 0.40215;
        controlPointCoords[1*4 + 2] = 0.0;
        controlPointCoords[1*4 + 3] = 1.0;

        // 3rd
        controlPointCoords[2*4 + 0] = -4.1291;
        controlPointCoords[2*4 + 1] = 0.72015;
        controlPointCoords[2*4 + 2] = 0.0;
        controlPointCoords[2*4 + 3] = 1.0;

        // 4th
        controlPointCoords[3*4 + 0] = -2.1779;
        controlPointCoords[3*4 + 1] = -0.68471;
        controlPointCoords[3*4 + 2] = 0.0;
        controlPointCoords[3*4 + 3] = 1.0;

        // 5th
        controlPointCoords[4*4 + 0] = 0.78705;
        controlPointCoords[4*4 + 1] = 1.5032;
        controlPointCoords[4*4 + 2] = 0.0;
        controlPointCoords[4*4 + 3] = 1.0;

        // 6th
        controlPointCoords[5*4 + 0] = 0.32737;
        controlPointCoords[5*4 + 1] = -0.21804;
        controlPointCoords[5*4 + 2] = 0.0;
        controlPointCoords[5*4 + 3] = 1.0;

        // 7th
        controlPointCoords[6*4 + 0] = 0.11198;
        controlPointCoords[6*4 + 1] = -1.103;
        controlPointCoords[6*4 + 2] = 0.0;
        controlPointCoords[6*4 + 3] = 1.0;

        theIGAPatchCurve = new IGAPatchCurve(0,p,noKnots,knotVector,noCPs,controlPointCoords);

        // Delete pointers
        delete[] knotVector;
        delete[] controlPointCoords;
	}

	void tearDown() {
        delete theIGAPatchCurve;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the constructor
	 ***********/
	void testConstructor() {
//        CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getId()==0);
        CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getPolynomialDegree() == 3);
        CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getNoKnots() == 11);
        for(int i = 0;i < 4; i++)
            CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getKnotVector()[i] == 0.0);
        CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getKnotVector()[4] == 0.30971);
        CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getKnotVector()[5] == 0.46851);
        CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getKnotVector()[6] == 0.81417);
        for(int i = 7;i < 11; i++)
            CPPUNIT_ASSERT(theIGAPatchCurve->getIGABasis()->getKnotVector()[i] == 1.0);
        CPPUNIT_ASSERT(theIGAPatchCurve->getNoControlPoints() == 7);
	}

	/***********************************************************************************************
     * \brief Test case: Test the computation of the basis functions and their derivatives
	 ***********/
    void testIGAPatchCurveBasisFunctionsAndDerivatives() {
        // Define a parametric coordinate
        double u = 0.1557;

        // Check the knot span index
        int knotSpanIndex = theIGAPatchCurve->getIGABasis()->findKnotSpan(u);
        CPPUNIT_ASSERT(knotSpanIndex == 3);

        // Get the polynomial order of the basis
        int p = theIGAPatchCurve->getIGABasis()->getPolynomialDegree();

        // Number of basis functions
        int noLocBasisFunctions = p + 1;

        // Get the number of derivatives for the basis functions
        int noDeriv = 1;

        // Define the expected solution in terms of the basis functions and their derivatives
        double expSolLocBasisFunctionsAndDerivatives[8]  = {0.122964876512325, 0.515334468590218, 0.329750229214093, 0.031950425683364,
                                                            -2.395264135685829, -1.354810285989855, 3.134459283267898, 0.615615138407786};

        // Compute the basis functions and their first derivatives
        double* locBasisFunctionsAndDerivatives = new double[noLocBasisFunctions*(noDeriv + 1)];
        theIGAPatchCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(locBasisFunctionsAndDerivatives, noDeriv, u, knotSpanIndex);

        // Verify the solution against the values from Matlab
        for(int i = 0; i < noLocBasisFunctions*(noDeriv + 1); i++)
            CPPUNIT_ASSERT(fabs(locBasisFunctionsAndDerivatives[i] - expSolLocBasisFunctionsAndDerivatives[i]) <= TolRel);

        // Delete pointers
        delete[] locBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Cartesian coordinates of a given parametric location
     ***********/
    void testComputeCartesianCoordinates() {
        // Number of Cartesian coordinates
        int noCoord = 3;

        // Initialize Cartesian coordinates
        double cartesianCoordinates[3];

        // Define a parametric coordinate
        double u = 0.358809;

        // Find the knot span where the parameter lies in
        int knotSpanIndex = theIGAPatchCurve->getIGABasis()->findKnotSpan(u);

        // Get the polynomial order of the basis
        int p = theIGAPatchCurve->getIGABasis()->getPolynomialDegree();

        // Get the number of basis functions
        int noLocBasisFunctions = p + 1;

        // Get the number of derivatives for the basis functions
        int noDeriv = 2;

        // Compute the basis functions and their first derivatives
        double* locBasisFunctionsAndDerivatives = new double[noLocBasisFunctions*(noDeriv + 1)];
        theIGAPatchCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(locBasisFunctionsAndDerivatives, noDeriv, u, knotSpanIndex);

        // Compute the Cartesian coordinates of the given parametric location
        theIGAPatchCurve->computeCartesianCoordinates(cartesianCoordinates, knotSpanIndex, locBasisFunctionsAndDerivatives);

        // Define the expected solution in terms of the Cartesian coordinates of the parametric location
        double expSolCartesianCoordinates[3] = {-3.417424547383805, 0.170160712217972, 0};

        // Verify the solution against the values from Matlab
        for(int iCoord = 0; iCoord < noCoord; iCoord++)
            CPPUNIT_ASSERT(fabs(cartesianCoordinates[iCoord] - expSolCartesianCoordinates[iCoord]) <= TolRel);

        // Delete pointers
        delete[] locBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vector and its derivatives
     ***********/
    void testComputeBaseVectorAndDerivatives() {
        // Number of Cartesian coordinates
        int noCoord = 3;

        // Define a parametric coordinate
        double u = 0.358809;

        // Find the knot span where the parameter lies in
        int knotSpanIndex = theIGAPatchCurve->getIGABasis()->findKnotSpan(u);

        // Get the polynomial order of the basis
        int p = theIGAPatchCurve->getIGABasis()->getPolynomialDegree();

        // Get the number of basis functions
        int noLocBasisFunctions = p + 1;

        // Get the number of derivatives for the basis functions
        int noDeriv = 2;

        // Compute the basis functions and their first derivatives
        double* locBasisFunctionsAndDerivatives = new double[noLocBasisFunctions*(noDeriv + 1)];
        theIGAPatchCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(locBasisFunctionsAndDerivatives, noDeriv, u, knotSpanIndex);

        // Compute the base vector and its derivatives
        double* baseVectorAndDerivatives = new double[noDeriv*noCoord];
        theIGAPatchCurve->computeBaseVectorAndDerivatives(baseVectorAndDerivatives, knotSpanIndex, locBasisFunctionsAndDerivatives, noDeriv - 1);

        // Define the expected solution in terms of the base vector and its derivatives
        double expSolBaseVectorAndDerivatives[6] = {7.519629163661582e+00,
                                                    -3.567924456192999e+00,
                                                    0.000000000000000e+00,
                                                    4.091421147205891e+00,
                                                    -3.268961096450674e+00,
                                                    0.000000000000000e+00};

        // Verify the solution against the values from Matlab
        for(int i = 0; i < noDeriv*noCoord; i++)
            CPPUNIT_ASSERT(fabs(baseVectorAndDerivatives[i] - expSolBaseVectorAndDerivatives[i]) <= TolRel);

        // Delete pointers
        delete[] locBasisFunctionsAndDerivatives;
        delete[] baseVectorAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the projection on a 2D curve for the first initial guess
     ***********/
    void testComputePointProjectionOn2DCurve1() {
        // Number of Cartesian coordinates
        int noCoordPrm = 2;

        // Define the point to be projected
        double P[3] = {-1.417, 0.28035};

        // Define the initial guess
        double u = 0.5;

        // Define Newton-Raphson convergence properties
        int maxIt = 30;
        double tol = 1e-5;

        // Compute the closest point projection on the curve
        bool isConvergent = theIGAPatchCurve->computePointProjectionOn2DCurve(u, P, noCoordPrm, maxIt, tol);

        // Verify the solution in terms of the convergence flag
        CPPUNIT_ASSERT(isConvergent == 1);

        // Define the expected solution in terms of the Cartesian coordinates of the projected point
        double expSolP[2] = {-1.359647350723798, 0.125633733514165};

        // Verify the result against the values from Matlab
        for(int iCoord = 0; iCoord < noCoordPrm; iCoord++)
            CPPUNIT_ASSERT(fabs(P[iCoord] - expSolP[iCoord]) < Tol);

        // Define the expected solution in terms of the parametric coordinate of the projected point
        double expSolu = 0.588309166996443;

        // Verify the result against the values from Matlab
        CPPUNIT_ASSERT(fabs(u - expSolu) < TolRel10000);
    }

    /***********************************************************************************************
     * \brief Test case: Test the projection on a 2D curve for the second initial guess
     ***********/
    void testComputePointProjectionOn2DCurve2() {
        // Number of Cartesian coordinates
        int noCoordPrm = 2;

        // Define the point to be projected
        double P[3] = {-1.417, 0.28035};

        // Define the initial guess
        double u = 0.95;

        // Define Newton-Raphson convergence properties
        int maxIt = 30;
        double tol = 1e-5;

        // Compute the closest point projection on the curve
        bool isConvergent = theIGAPatchCurve->computePointProjectionOn2DCurve(u, P, noCoordPrm, maxIt, tol);

        // Verify the solution in terms of the convergence flag
        CPPUNIT_ASSERT(isConvergent == 1);

        // Define the expected solution in terms of the Cartesian coordinates of the projected point
        double expSolP[2] = {0.359603628496996, -0.000698427792315};

        // Verify the result against the values from Matlab
        for(int iCoord = 0; iCoord < noCoordPrm; iCoord++)
            CPPUNIT_ASSERT(fabs(P[iCoord] - expSolP[iCoord]) < Tol);

        // Define the expected solution in terms of the parametric coordinate of the projected point
        double expSolu = 0.906662493347622;

        // Verify the result against the values from Matlab
        CPPUNIT_ASSERT(fabs(u - expSolu) < TolRel1000000);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vector and its derivatives for memory leakage
     ***********/
    void testComputeBaseVectorAndDerivatives4Leakage() {
        // Number of Cartesian coordinates
        int noCoord = 3;

        // Define a parametric coordinate
        double u = 0.358809;

        // Find the knot span where the parameter lies in
        int knotSpanIndex = theIGAPatchCurve->getIGABasis()->findKnotSpan(u);

        // Get the polynomial order of the basis
        int p = theIGAPatchCurve->getIGABasis()->getPolynomialDegree();

        // Get the number of basis functions
        int noLocBasisFunctions = p + 1;

        // Get the number of derivatives for the basis functions
        int noDeriv = 2;

        int noIterations = 1e9;
        for(int i = 0; i < noIterations; i++) {
            // Compute the basis functions and their first derivatives
            double* locBasisFunctionsAndDerivatives = new double[noLocBasisFunctions*(noDeriv + 1)];
            theIGAPatchCurve->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(locBasisFunctionsAndDerivatives, noDeriv, u, knotSpanIndex);

            // Compute the base vector and its derivatives
            double* baseVectorAndDerivatives = new double[noDeriv*noCoord];
            theIGAPatchCurve->computeBaseVectorAndDerivatives(baseVectorAndDerivatives, knotSpanIndex, locBasisFunctionsAndDerivatives, noDeriv - 1);

            // Delete pointers
            delete[] locBasisFunctionsAndDerivatives;
            delete[] baseVectorAndDerivatives;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the projection on the 2D curve for leakage
     ***********/
    void testComputePointProjectionOn2DCurve4Leakage() {
        // Initialize auxiliary variables
        bool isConvergent;

        // Number of Cartesian coordinates
        int noCoordPrm = 2;

        // Define the point to be projected
        double P[3] = {-1.417, 0.28035};

        // Define the initial guess
        double u = 0.5;

        // Define Newton-Raphson convergence properties
        int maxIt = 30;
        double tol = 1e-5;

        // Compute the closest point projection on the curve
        int noIterations = 1e9;
        for(int i = 0; i < noIterations; i++)
            isConvergent = theIGAPatchCurve->computePointProjectionOn2DCurve(u, P, noCoordPrm, maxIt, tol);
    }

// Make the tests
CPPUNIT_TEST_SUITE(TestIGAPatchCurve);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testComputeBaseVectorAndDerivatives);
    CPPUNIT_TEST(testComputeCartesianCoordinates);
    CPPUNIT_TEST(testComputeBaseVectorAndDerivatives);
    CPPUNIT_TEST(testComputePointProjectionOn2DCurve1);
    CPPUNIT_TEST(testComputePointProjectionOn2DCurve2);

// Make the tests for leakage
    // CPPUNIT_TEST(testComputeBaseVectorAndDerivatives4Leakage);
    // CPPUNIT_TEST(testComputePointProjectionOn2DCurve4Leakage);

	CPPUNIT_TEST_SUITE_END()
	;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatchCurve);

