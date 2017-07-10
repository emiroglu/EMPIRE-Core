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

class TestIGAPatchSurfaceTrimming: public CppUnit::TestFixture {

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
		const int direction = 0;

		// Provide an id for the basis
		const int IDBasis = 1;

		// The polynomial degree
		const int pDegree = 3;

		// Number of knots
		const int uNoKnots = 27;

		// The knot vectors
		double uKnotVector[uNoKnots] = { 0.0, 0.0, 0.0, 0.0, 0.008673617703141023, 0.022387811450882932, 0.039735046857164974, 0.06184851981804492, 0.08786937292746798, 0.1185352424040215, 0.15322971321658557, 0.19250118763803997, 0.23586927615374506, 0.2837708955818746, 0.33581260180072076, 0.3923576711928459, 0.45307299511483307, 0.5182695279168437, 0.5876584695419718, 0.6615116630705355, 0.7395742223988048, 0.8220876382016415, 0.9088238152330518, 1.0, 1.0, 1.0, 1.0};
		
		// The Control Point net
		const int uNoControlPoints = uNoKnots - pDegree - 1;
		
		// Control Points for a NURBS
		double controlPointNet[uNoControlPoints*4] = 
					       { 0.0,			 0.0,			0.0,	1.0,
						 0.03711725733156428,	 0.042688513773407946,	0.0,	1.0,
						 0.18128490012859777,	 0.2015404446885209,	0.0,	1.0,
						-0.14018240692059364,	-0.29428358253262943,	0.0,	1.0,
						 0.4025638661429545,	 0.35024387491125203,	0.0,	1.0,
						-0.2506977318483484,	-0.5087560176292382,	0.0,	1.0,
						 0.6053600837771562,	 0.45482419877466773,	0.0,	1.0,
						-0.352029338412692,	-0.7078157439904523,	0.0,	1.0,
						 0.8048537252505352,	 0.5539426255559453,	0.0,	1.0,
						-0.4519405440215286,	-0.9060985482545199,	0.0,	1.0,
						 1.0041413411962206,	 0.6530513405768528,	0.0,	1.0,
						-0.5517187135874752,	-1.1049102552107424,	0.0,	1.0,
						 1.2035698332094393,	 0.7524572591630635,	0.0,	1.0,
						-0.6515148911699897,	-1.3040934520772662,	0.0,	1.0,
						 1.4030950325255178,	 0.8520183481137549,	0.0,	1.0,
						-0.7512083064268196,	-1.5033630672439289,	0.0,	1.0,
						 1.602148870014657,	 0.9511105698078708,	0.0,	1.0,
						-0.8486228325121628,	-1.7004019069582168,	0.0,	1.0,
						 1.7919217828929315,	 1.0406544585008441,	0.0,	1.0,
						-0.9082676483904522,	-1.8586489184805504,	0.0,	1.0,
						 1.8298779316639826,	 0.9739756886167465,	0.0,	1.0,
						 0.40526878743025746,	-0.5648480610926698,	0.0,	1.0,
						 0.0,			-1.0,			0.0,	1.0};

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
	 * \brief Test case: Test computeIntersectionsWithKnotBisection
	 ***********/
	void testComputeTrimmingCurveIntersectionsWithKnotBisection() {
	  
		double referenceIntersectionsUV[80] = {9.652790268736409e-01,  5.000047881316947e-02,
						       9.479975322156601e-01,  5.000076591968322e-02,
						       8.596546618057613e-01,  5.000049536997528e-02,
						       8.469506619609048e-01,  5.000018496834335e-02,
						       7.600118680222836e-01,  5.000082334771844e-02,
						       7.469608648486293e-01,  5.000078017353132e-02,
						       6.607759534905656e-01,  4.999911088078902e-02,
						       6.471011257716630e-01,  5.000036357416815e-02,
						       5.617451773630358e-01,  4.999993121513078e-02,
						       5.473497446603758e-01,  4.999904871752374e-02,
						       4.629757435588678e-01,  5.000056029441495e-02,
						       4.477903419093880e-01,  5.000017489796048e-02,
						       3.645895395110711e-01,  4.999909533839611e-02,
						       3.485973643384774e-01,  5.000029178781272e-02,
						       2.668358812678921e-01,  4.999971910309909e-02,
						       2.501173631924725e-01,  4.999904824846019e-02,
						       1.705447915704064e-01,  4.999977577547785e-02,
						       1.523160000903682e-01,  5.000004050918205e-02,
						       7.899320413070741e-02,  5.000058986896855e-02,
						       4.469512884086848e-02,  5.000051570989300e-02,
						       4.999986807714255e-02, -9.462900837797903e-01,
						       5.000085403248532e-02, -8.563078941173277e-01,
						       5.000062340654288e-02, -8.408832547793020e-01,
						       4.999972599376677e-02, -7.558567576972531e-01,
						       4.999982544229256e-02, -7.393594376767862e-01,
						       5.000061592275257e-02, -6.561675201075120e-01,
						       5.000078626225002e-02, -6.382799760380714e-01,
						       5.000002847395595e-02, -5.566100071160699e-01,
						       5.000061508460731e-02, -5.369462634207405e-01,
						       4.999915481936426e-02, -4.571880067779903e-01,
						       4.999987024413353e-02, -4.351396268954141e-01,
						       5.000062191785445e-02, -3.579816766484307e-01,
						       5.000006303389427e-02, -3.325115462953434e-01,
						       4.999946100788041e-02, -2.591602154025976e-01,
						       5.000016541867213e-02, -2.282837044713583e-01,
						       5.000029686498839e-02, -1.608609332963082e-01,
						       4.999975019801586e-02, -1.204287905907703e-01,
						       4.999984170848370e-02, -6.109982764394972e-02,
						       5.000070119537275e-02, -1.625974383110812e-03,
						       4.999978403725074e-02,  5.570488484748966e-02};
						       
		double referenceIntersectionsUTilde[40] = {9.305734073414522e-01,
							   9.010774311085488e-01,
							   7.520563462201286e-01,
							   7.294755823471967e-01,
							   5.991635260452104e-01,
							   5.779676998362822e-01,
							   4.639320373535157e-01,
							   4.440250883914517e-01,
							   3.460063934326173e-01,
							   3.274699798241215e-01,
							   2.453495074983381e-01,
							   2.282979628619026e-01,
							   1.619341233197381e-01,
							   1.465156779569738e-01,
							   9.570462884087756e-02,
							   8.211651970358454e-02,
							   4.645695405847887e-02,
							   3.492563991156171e-02,
							   1.325076431740871e-02,
							   3.155384991247042e-03,
							   9.963611153995289e-01,
							   8.364172542796416e-01,
							   8.119467720421896e-01,
							   6.735737220204276e-01,
							   6.510873682358684e-01,
							   5.296314463895910e-01,
							   5.083840314079735e-01,
							   4.030711230109720e-01,
							   3.830946751881789e-01,
							   2.937931068689920e-01,
							   2.751726262709673e-01,
							   2.017704988026441e-01,
							   1.846120497759651e-01,
							   1.269637286219885e-01,
							   1.113963407628676e-01,
							   6.928007608513843e-02,
							   5.542687808766085e-02,
							   2.840849692633041e-02,
							   1.589647503089962e-02,
							   3.524780236337090e-03};
	  
		const IGAPatchSurfaceTrimmingLoop* theIGAPatchSurfaceTrimmingLoop = &theIGAPatchSurfaceTrimming->getFirstLoop();
		const IGAPatchCurve* theIGAPatchCurve = &theIGAPatchSurfaceTrimmingLoop->getIGACurve(0);
		
		std::vector<double> uvSurface;
		std::vector<double> uTilde;
		unsigned int dir_u=0;
		unsigned int dir_v=1;
		
		double knot_u1 = 0;
		double knot_u2 = 0.05;
		double knot_v1 = 0.1;
		double knot_v2 = 0.05;
		
		// Intersect the trimming curve with the parameter line v=0.05 in u direction
		theIGAPatchCurve->computeIntersectionsWithKnotBisection(uvSurface,uTilde,dir_u,knot_v2);
		// Intersect the trimming curve with the parameter line v=0.05 in v direction
		theIGAPatchCurve->computeIntersectionsWithKnotBisection(uvSurface,uTilde,dir_v,knot_u2);
		
		// Assert number of intersections (u and v)
		CPPUNIT_ASSERT(uvSurface.size()==80);
		
		// Assert number of intersections (uTilde)
		CPPUNIT_ASSERT(uTilde.size()==40);
		
		// Assert the intersection parameters in the patch parameter space
		for (int i=0; i<uvSurface.size(); i++)
		  CPPUNIT_ASSERT(fabs(uvSurface[i]-referenceIntersectionsUV[i])<Tol);
		
		for (int i=0; i<uTilde.size(); i++)
		  CPPUNIT_ASSERT(fabs(uTilde[i]-referenceIntersectionsUTilde[i])<Tol);
		
	}

	// Make the tests
	CPPUNIT_TEST_SUITE(TestIGAPatchSurfaceTrimming);
	CPPUNIT_TEST(testComputeTrimmingCurveIntersectionsWithKnotBisection);
	CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatchSurfaceTrimming);

