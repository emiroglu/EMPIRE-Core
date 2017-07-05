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

// Inclusion of user-defined libraries
#include "IGAPatchSurface.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class IGAPatchSurface
 ***********/

class TestIGAPatchCurve: public CppUnit::TestFixture {

private:
  
	IGAPatchSurfaceTrimming* theIGAPatchSurfaceTrimming;
	
	double Tol;
	double relTol;
	double TolDeriv;

public:
	void setUp() {
		// Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
		Tol = 1e-15;

		// Assign a relaxed tolerance value (corresponding to maximum accuracy provided by MATLAB) for the Newton-Rapson iteration error (accumulative error appears here)
		relTol = 1e-14;

		// Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the derivative functional values
		TolDeriv = 1e-13;
		
		// Curve direction
		int direction = 0;

		// Provide an id for the basis
		int IDBasis = 1;

		// The polynomial degree
		int pDegree = 3;

		// Number of knots
		const int uNoKnots = 50;

		// The knot vectors
		double uKnotVector[uNoKnots] = { -23.20663072597694, -23.20663072597694, -23.20663072597694, -23.20663072597694, -21.659522010911811, -21.659522010911811, -21.659522010911811, -20.112413295846682, -20.112413295846682, -20.112413295846682, -18.56530458078155, -18.56530458078155, -18.56530458078155, -17.018195865716422, -17.018195865716422, -17.018195865716422, -15.471087150651293, -15.471087150651293, -15.471087150651293, -13.923978435586164, -13.923978435586164, -13.923978435586164, -12.376869720521034, -12.376869720521034, -12.376869720521034, -10.829761005455905, -10.829761005455905, -10.829761005455905, -9.2826522903907751, -9.2826522903907751, -9.2826522903907751, -7.7355435753256465, -7.7355435753256465, -7.7355435753256465, -6.188434860260517, -6.188434860260517, -6.188434860260517, -4.6413261451953876, -4.6413261451953876, -4.6413261451953876, -3.0942174301302585, -3.0942174301302585, -3.0942174301302585, -1.5471087150651293, -1.5471087150651293, -1.5471087150651293, 0, 0, 0, 0};
		
		// The Control Point net
		const int uNoControlPoints = uNoKnots - pDegree - 1;
		
		// Control Points for a NURBS
		double controlPointNet[46*4] = { 10.0, 3.0, 0.0, 1.0,
						 9.0, 4.0, 0.0, 1.0,
						 8.5, 3.5, 0.0, 1.0,
						 8.0833333333333339, 3.0833333333333335, 0.0, 1.0,
						 7.666666666666667, 2.666666666666667, 0.0, 1.0,
						 7.3333333333333339, 2.3333333333333339, 0.0, 1.0,
						 7.0000000000000018, 2.0000000000000009, 0.0, 1.0,
						 6.666666666666667, 1.666666666666667, 0.0, 1.0,
						 6.333333333333333, 1.333333333333333, 0.0, 1.0,
						 5.9999999999999991, 1.3333333333333335, 0.0, 1.0,
						 5.6666666666666661, 1.3333333333333339, 0.0, 1.0,
						 5.333333333333333, 1.666666666666667, 0.0, 1.0,
						 5.3333333333333339, 1.8333333333333335, 0.0, 1.0,
						 5.3333333333333339, 2.0, 0.0, 1.0,
						 5.666666666666667, 2.0, 0.0, 1.0,
						 6.0, 2.1666666666666665, 0.0, 1.0,
						 6.3333333333333339, 2.333333333333333, 0.0, 1.0,
						 6.666666666666667, 2.6666666666666661, 0.0, 1.0,
						 6.6666666666666661, 2.9999999999999991, 0.0, 1.0,
						 6.6666666666666661, 3.333333333333333, 0.0, 1.0,
						 6.333333333333333, 3.666666666666667, 0.0, 1.0,
						 5.9999999999999991, 3.6666666666666665, 0.0, 1.0,
						 5.6666666666666661, 3.6666666666666661, 0.0, 1.0,
						 5.3333333333333339, 3.3333333333333335, 0.0, 1.0,
						 5.0, 3.333333333333333, 0.0, 1.0,
						 4.666666666666667, 3.333333333333333, 0.0, 1.0,
						 4.333333333333333, 3.666666666666667, 0.0, 1.0,
						 4.0, 3.666666666666667, 0.0, 1.0,
						 3.6666666666666661, 3.6666666666666661, 0.0, 1.0,
						 3.3333333333333335, 3.3333333333333335, 0.0, 1.0,
						 3.166666666666667, 3.0, 0.0, 1.0,
						 3.0, 2.666666666666667, 0.0, 1.0,
						 3.0, 2.333333333333333, 0.0, 1.0,
						 2.8333333333333335, 1.9999999999999998, 0.0, 1.0,
						 2.6666666666666665, 1.6666666666666665, 0.0, 1.0,
						 2.333333333333333, 1.3333333333333333, 0.0, 1.0,
						 1.9999999999999998, 1.3333333333333335, 0.0, 1.0,
						 1.6666666666666665, 1.3333333333333333, 0.0, 1.0,
						 1.3333333333333333, 1.6666666666666667, 0.0, 1.0,
						 1.3333333333333335, 1.8333333333333335, 0.0, 1.0,
						 1.333333333333333, 2.0, 0.0, 1.0,
						 1.6666666666666665, 2.0, 0.0, 1.0,
						 1.8333333333333333, 2.25, 0.0, 1.0,
						 2.0, 2.5, 0.0, 1.0,
						 2.0, 3.0, 0.0, 1.0,
						 0.0, 4.0, 0.0, 1.0};
		

		// Test just one object of the class (that works pretty also)
		theIGAPatchSurfaceTrimming = new IGAPatchSurfaceTrimming();
		theIGAPatchSurfaceTrimming->addTrimLoop(0,1);
		theIGAPatchSurfaceTrimming->addTrimCurve(direction, IDBasis, pDegree, uNoKnots, uKnotVector, uNoControlPoints, controlPointNet);
		theIGAPatchSurfaceTrimming->linearizeLoops();

	}

	void tearDown() {
		delete theIGAPatchSurfaceTrimming;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the constructor
	 ***********/
	void testComputeIntersectionsWithKnotBisection() {
		const IGAPatchSurfaceTrimmingLoop* theIGAPatchSurfaceTrimmingLoop = &theIGAPatchSurfaceTrimming->getFirstLoop();
		const IGAPatchCurve* theIGAPatchCurve = &theIGAPatchSurfaceTrimmingLoop->getIGACurve(0);
		
		std::vector<double> uvSurface;
		std::vector<double> uTilde;
		unsigned int dir=0;
		double knot = 2.5;
		
		theIGAPatchCurve->computeIntersectionsWithKnotBisection(uvSurface,uTilde,dir,knot);
		
	}

	// Make the tests
	CPPUNIT_TEST_SUITE(TestIGAPatchCurve);
	CPPUNIT_TEST(testComputeIntersectionsWithKnotBisection);
	CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatchCurve);

