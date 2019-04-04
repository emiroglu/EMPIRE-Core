
///*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
//*  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
//*
//*  All rights reserved.
//*
//*  This file is part of EMPIRE.
//*
//*  EMPIRE is free software: you can redistribute it and/or modify
//*  it under the terms of the GNU General Public License as published by
//*  the Free Software Foundation, either version 3 of the License, or
//*  (at your option) any later version.
//*
//*  EMPIRE is distributed in the hope that it will be useful,
//*  but WITHOUT ANY WARRANTY; without even the implied warranty of
//*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//*  GNU General Public License for more details.
//*
//*  You should have received a copy of the GNU General Public License
//*  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
//*/
//// inclusion of standard libraries   (only if really necessary here in *.h)
//#include "cppunit/TestFixture.h"
//#include "cppunit/TestAssert.h"
//#include "cppunit/extensions/HelperMacros.h"
//#include <iostream>
//#include <string>
//#include <math.h>
//#include "MathLibrary.h"

//// Inclusion of user-defined libraries
//#include "IGAPatchSurface.h"
//#include "MatrixVectorMath.h"

//using namespace std;

//namespace EMPIRE {

///********//**
//* \brief Test the class IGAPatchSurface
//***********/

//class TestIGAPatchSurfaceCurvedPlate: public CppUnit::TestFixture {

//private:
//   IGAPatchSurface* theIGAPatchSurface;
//   double Tol;
//   double relTol;
//   double TolDeriv;

//public:
//   void setUp() {
//       // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
//       Tol = 1e-15;

//       // Assign a relaxed tolerance value (corresponding to maximum accuracy provided by MATLAB) for the Newton-Rapson iteration error (accumulative error appears here)
//       relTol = 1e-14;

//       // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the derivative functional values
//       TolDeriv = 1e-13;

//       // Provide an id for the basis
//       int id_basis = 3;

//       // Provide an id for the patch itself
//       int id_patch = 1;


//       // The polynomial degrees
//       int p = 2;
//       int q = 2;

//       // Number of knots in both directions
//       int uNoKnots = 7;
//       int vNoKnots = 6;

//       double* uKnotVector = new double[uNoKnots];
//       uKnotVector[0]=0;
//       uKnotVector[1]=0;
//       uKnotVector[2]=0;
//       uKnotVector[3]=0.5;
//       uKnotVector[4]=1;
//       uKnotVector[5]=1;
//       uKnotVector[6]=1;

//       double* vKnotVector = new double[vNoKnots];
//       vKnotVector[0]=0;
//       vKnotVector[1]=0;
//       vKnotVector[2]=0;
//       vKnotVector[3]=1;
//       vKnotVector[4]=1;
//       vKnotVector[5]=1;


//       // The Control Point net
//       int uNoControlPoints = uNoKnots - p - 1;
//       int vNoControlPoints = vNoKnots - q - 1;
//       IGAControlPoint** controlPointNet = new IGAControlPoint*[uNoControlPoints
//               * vNoControlPoints];
//       double* controlPointWeights = new double[uNoControlPoints
//               * vNoControlPoints];


//       // the numbering of the gauss points is first in xi direction then eta
//       // First row
//       controlPointNet[0] = new IGAControlPoint(0, 0.0, 0.0, 0.0,
//                       1.0);
//       controlPointNet[1] = new IGAControlPoint(3, 1.0, 0.0, 0.0,
//                       2.0);
//       controlPointNet[2] = new IGAControlPoint(6, 2.0, 0.0, 0.0,
//                       1.0);
//       controlPointNet[3] = new IGAControlPoint(9, 3.0, 0.0, 0.0,
//                       2.0);

//       // Second row
//       controlPointNet[4] = new IGAControlPoint(1, 0.0, 1.5, -0.25,
//                       2.0);
//       controlPointNet[5] = new IGAControlPoint(4, 1.0, 1.0, 1.0,
//                       1.0);
//       controlPointNet[6] = new IGAControlPoint(7, 1.75, 1.0, 0.0,
//                       2.0);
//       controlPointNet[7] = new IGAControlPoint(10, 3.0, 1.0, 0.0,
//                       1.0);

//       // Third row
//       controlPointNet[8] = new IGAControlPoint(2, 0.0, 2.0, 0.0,
//                       1.0);
//       controlPointNet[9] = new IGAControlPoint(5, 1.0, 2.0, 0.5,
//                       2.0);
//       controlPointNet[10] = new IGAControlPoint(8, 2.0, 2.0, 0.4,
//                       1.0);
//       controlPointNet[11] = new IGAControlPoint(11, 3.0, 2.0, 0.0,
//                       2.0);


//       theIGAPatchSurface = new IGAPatchSurface(id_basis, p, uNoKnots, uKnotVector, q, vNoKnots, vKnotVector,
//               uNoControlPoints, vNoControlPoints, controlPointNet);
//   }

//   void tearDown() {
//       delete theIGAPatchSurface;
//   }

//   /***********************************************************************************************
//    * \brief Test case: Test the constructor
//    ***********/
//   void testConstructor() {
//       CPPUNIT_ASSERT(theIGAPatchSurface->getIGABasis()->getId()==3);
//       CPPUNIT_ASSERT(
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()==2);
//       CPPUNIT_ASSERT(
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()==2);
//       CPPUNIT_ASSERT(
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getNoKnots()==7);
//       CPPUNIT_ASSERT(
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getNoKnots()==6);
//       for (int i = 0; i < 3; i++)
//           CPPUNIT_ASSERT(
//                   theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 0);
//       for (int i = 3; i < 4; i++)
//           CPPUNIT_ASSERT(
//                   theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 0.5);
//       for (int i = 4; i < 7; i++)
//           CPPUNIT_ASSERT(
//                   theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 1);

//       for (int i = 0; i < 3; i++)
//           CPPUNIT_ASSERT(
//                   theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == 0);
//       for (int i = 3; i < 6; i++)
//           CPPUNIT_ASSERT(
//                   theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == 1);

//       CPPUNIT_ASSERT(theIGAPatchSurface->getUNoControlPoints()==4);
//       CPPUNIT_ASSERT(theIGAPatchSurface->getVNoControlPoints()==3);

//   }

//   /***********************************************************************************************
//    * \brief Test case: Test the knot span
//    ***********/
//   void testIGAPatchSurfaceKnotSpan() {

//       // The parametric coordinates on the Spline parameter space
//       double u = 0.4;
//       double v = 0.65;
//       int uCorrectknotSpan = 2;
//       int vCorrectknotSpan = 2;

//       // Check the find knot span function
//       CPPUNIT_ASSERT(
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u)==uCorrectknotSpan);
//       CPPUNIT_ASSERT(
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v)==vCorrectknotSpan);
//   }

//   void testIGAPatchSurfaceComputeSurfaceNormalAndDerivatives() {

//       // The parameters on the IGA surface and their knot span indices
//       double u = 0.4;
//       double v = 0.65;

//       int uKnotSpan =
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
//                       u);
//       int vKnotSpan =
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
//                       v);
//       // Get the polynomial degrees of the patch
//       int pDegree =
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
//       int qDegree =
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

//       // Decide on absolute degree of the partial derivatives for the base vectors
//       int derivDegreeBaseVec = 1;

//       // Initialize the index of the base vector
//       int indexBaseVct = 0;

//       // Initialize the indices related to the derivatives of each component
//       int uDeriv = 1;
//       int vDeriv = 1;

//       // Indices related to the base vectors
//       int indexUBaseVec = 0;
//       int indexVBaseVec = 1;

//       // Compute the local basis functions and their derivatives
//       int noBasisFcts = (pDegree + 1) * (qDegree + 1);
//       int derivDegree = derivDegreeBaseVec + 1;
//       double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
//               * (derivDegree + 2) * noBasisFcts / 2];
//       theIGAPatchSurface->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
//               localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v,
//               vKnotSpan);

//       // Compute the base vectors and their derivatives
//       int noBaseVec = 2;
//       int noCoord = 3;
//       double* baseVctsAndDerivs = new double[(derivDegreeBaseVec + 1)
//               * (derivDegreeBaseVec + 2) * noCoord * noBaseVec / 2];
//       theIGAPatchSurface->computeBaseVectorsAndDerivatives(baseVctsAndDerivs,
//               localBasisFunctionsAndDerivatives, derivDegreeBaseVec,
//               uKnotSpan, vKnotSpan);

//       vector<double> normalAndDerivatives;

//       theIGAPatchSurface->computeSurfaceNormalAndDerivatives(normalAndDerivatives, baseVctsAndDerivs);

//       // normal and derivatives calculated with matlab
//       double correctNormalVctAndDerivs[] = {
//           0.037639568219582482 , -0.170756724002525034 , 0.984594030101832995 ,
//           3.607507521909729231 , -0.014399785093197726 , -0.140406991486770050 ,
//           -0.087369807255650897 , 0.393120217785823900 , 0.071518290986910094 };

//       for (int i = 0; i < 9; i++)
//           CPPUNIT_ASSERT(
//                   fabs(normalAndDerivatives[i]-correctNormalVctAndDerivs[i])<=Tol);

//       delete[] localBasisFunctionsAndDerivatives;
//       delete[] baseVctsAndDerivs;

//   }

//   void testIGAPatchSurfaceComputeCurvatureAtPoint() {

//       double u = 0.8;
//       double v = 0;

//       int uKnotSpan =
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
//       int vKnotSpan =
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

//       int p = theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
//       int q = theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
//       int noLocalBasisFunctions = (p+1) * (q+1);

//       double* B_disp = new double[noLocalBasisFunctions];

//       double* tang = new double[3];
//       tang[0] = 1.0;
//       tang[1] = 0.0;
//       tang[2] = 0.0;



//       vector<vector<double> > B_rot;
//       vector<double> normalAndDerivs;
//       theIGAPatchSurface->computeCurvatureAtPoint(B_rot, B_disp, normalAndDerivs, u,uKnotSpan,v,vKnotSpan, tang);

//        double correct_Brot[] = {
//            0.000000000000000000, 0.011803278688524601,
//            -0.230163934426229733, 0.000000000000000000,
//            0.027540983606557396, -0.537049180327869302,
//            0.000000000000000000, 0.011803278688524576,
//            -0.230163934426229122, 0.000000000000000000,
//            -0.002622950819672129, 0.051147540983606528,
//            0.000000000000000000, -0.036721311475409815,
//            0.716065573770491626, 0.000000000000000000,
//            -0.011803278688524587, 0.230163934426229511,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.015758385222510281, -0.307288511838950451,
//            0.000000000000000000, 0.019697981528137867,
//            -0.384110639798688314, 0.000000000000000000,
//            -0.035456366750648088, 0.691399151637638321,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000,
//            0.000000000000000000, 0.000000000000000000};

//        int counter = 0;
//        for(int i = 0 ; i < B_rot.size() ; i++) {
//            for(int j = 0 ; j < B_rot[i].size() ; j++) {
//                CPPUNIT_ASSERT(
//                        fabs(B_rot[i][j]-correct_Brot[counter++])<=Tol);
//            }
//        }

//       delete[] B_disp;
//       delete[] tang;

//   }

//   void testIGAPatchSurfaceComputeCurvatureAtPoint4Leakage() {

//       double u = 0.8;
//       double v = 0;

//       int uKnotSpan =
//               theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
//       int vKnotSpan =
//               theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

//       int p = theIGAPatchSurface->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
//       int q = theIGAPatchSurface->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
//       int noLocalBasisFunctions = (p+1) * (q+1);

//       double* B_disp = new double[noLocalBasisFunctions];

//       double* tang = new double[3];
//       tang[0]=1;
//       tang[1]=0;
//       tang[2]=0;

//       int derivDegree = 2;

//       int noCoord=3;
//       int noBaseVec=2;


//       for(int i = 0 ; i < 10E6 ; i++) {
//           vector<vector<double> > B_rot;
//           vector<double> normalAndDerivs;
//           theIGAPatchSurface->computeCurvatureAtPoint(B_rot, B_disp, normalAndDerivs, u, uKnotSpan, v, vKnotSpan, tang);

//       }
//       delete[] B_disp;
//       delete[] tang;
//   }

//// Make the tests
//CPPUNIT_TEST_SUITE(TestIGAPatchSurfaceCurvedPlate);
//   CPPUNIT_TEST(testConstructor);
//   CPPUNIT_TEST(testIGAPatchSurfaceKnotSpan);
//   CPPUNIT_TEST(testIGAPatchSurfaceComputeSurfaceNormalAndDerivatives);
//   CPPUNIT_TEST(testIGAPatchSurfaceComputeCurvatureAtPoint);
////   CPPUNIT_TEST(testIGAPatchSurfaceComputeCurvatureAtPoint4Leakage);

//   CPPUNIT_TEST_SUITE_END()
//   ;
//}
//;

//} /* namespace EMPIRE */

//CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatchSurfaceCurvedPlate);



