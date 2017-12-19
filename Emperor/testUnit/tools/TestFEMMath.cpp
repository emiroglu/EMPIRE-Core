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
#include <math.h>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include "AuxiliaryParameters.h"
#include "MathLibrary.h"

namespace EMPIRE {
using namespace std;
using namespace MathLibrary;

/********//**
 * \brief This class manages tests the MathLibrary of EMPIRE
 **************************************************************************************************/
class TestFEMMath: public CppUnit::TestFixture {
private:
    // Number of Gauss points for the integration in the bi-unit interval
    int noGPsBiunitInterval;

    // Basic tolerance
    double tolBasic;

    // Absolute tolerance
    double tolAbs;
public:
    /***********************************************************************************************
     * \brief Set up the Gauss quadrature
     * \author Andreas Apostolatos
     ***********/
    void setUp() {
        noGPsBiunitInterval = 50;
        tolBasic = 1e-1;
        tolAbs = 1e-15;
    }

    /***********************************************************************************************
     * \brief Delete test vectors and matrices
     * \author Andreas Apostolatos
     ***********/
    void tearDown() {
    }

    /***********************************************************************************************
     * \brief Solve the integral of function xlnx in [1 5] interval numerically for all quadratures
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnBiunitInterval() {
        // Initialize auxiliary values
        double integral;
        double gaussPointValueOnParentSpace;
        double gaussPointValue;
        double gaussPointWeight;
        double TolWeighted;
        double integralAnalytical = 25.0*log(5.0)/2.0 - 6.0;

        // Loop over all quadratures
        for (int i = 0; i < noGPsBiunitInterval; i++){
            // Create a Gauss quadrature rule
            IGAGaussQuadratureOnBiunitInterval* theIGAGaussQuadratureOnBiunitInterval = new IGAGaussQuadratureOnBiunitInterval(i + 1);

            // Initialize the integral value
            integral = 0.0;

            // Loop over all Gauss points and compute the integral value numerically
            for (int j = 0; j < i + 1; j++){
                gaussPointValueOnParentSpace = *(theIGAGaussQuadratureOnBiunitInterval->getGaussPoint(j));
                gaussPointWeight = theIGAGaussQuadratureOnBiunitInterval->getGaussWeight(j);
                gaussPointValue = 2*gaussPointValueOnParentSpace + 3;
                integral += gaussPointValue*log(gaussPointValue)*2*gaussPointWeight;
            }

            // Verify the solution against the analytical one
            if (i <= 15)
                TolWeighted = tolBasic*pow(10.0, - i + 1);
            else
                TolWeighted = 1e-14;
            CPPUNIT_ASSERT(fabs(integralAnalytical - integral) <= TolWeighted);

            // Delete pointers
            delete theIGAGaussQuadratureOnBiunitInterval;
        }
    }

    /***********************************************************************************************
     * \brief Solve the integral of function sin(pi*xi/2)*cos(2*pi*eta) in [-1 1]x[-1 1] interval numerically for all quadratures
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnBiunitQuadrilateral() {
        // Initialize auxiliary values
        int noGPs;
        double integral;
        double gaussPointXiValue;
        double gaussPointEtaValue;
        double gaussPointWeight;
        double integralAnalytical = 8*(1.0/3.0)*(1.0/M_PI);

        // Loop over all quadratures
        for (int i = 0; i < 10; i++) {
            // Number of Gauss points for each quadrature
            noGPs = pow(i + 1,2.0);

            // Create a Gauss quadrature rule
            IGAGaussQuadratureOnBiunitQuadrilateral* theIGAGaussQuadratureOnBiunitQuadrilateral = new IGAGaussQuadratureOnBiunitQuadrilateral(noGPs);

            // Initialize the integral value
            integral = 0.0;

            // Loop over all Gauss points and compute the integral value numerically
            for (int j = 0; j < noGPs; j++){
                gaussPointXiValue = *theIGAGaussQuadratureOnBiunitQuadrilateral->getGaussPoint(j);
                gaussPointEtaValue = *(theIGAGaussQuadratureOnBiunitQuadrilateral->getGaussPoint(j) + 1);
                gaussPointWeight = theIGAGaussQuadratureOnBiunitQuadrilateral->getGaussWeight(j);
                integral += pow(gaussPointXiValue,2.0)*cos(M_PI*gaussPointEtaValue/2.0)*gaussPointWeight;
            }

            // Verify the solution against the analytical one
            CPPUNIT_ASSERT(fabs(integralAnalytical - integral) <= 7.0 * pow(10.0,-2.0*i));

            delete theIGAGaussQuadratureOnBiunitQuadrilateral;
        }
    }

    /***********************************************************************************************
     * \brief Solve the integral of function ln(zeta1 + 1)*sin(pi*zeta2/2) in the canonical triangle using the symmetric rule
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnCanonicalTriangleUsingTheSymmetricRule() {
        // Initialize auxiliary values
        int noGPs;
        double integral;
        double gaussPointZeta1Value;
        double gaussPointZeta2Value;
        double gaussPointWeight;
        double integralAnalytical = 0.050909793169365;
        double tolWeighted = 1e-1;
        int noGPsPerPolOrder[8] = {1, 3, 4, 6, 7, 12, 13, 16};

        // Loop over all polynomial orders
        for (int i = 0; i < 8; i++) {
            // Number of Gauss points for each quadrature
            noGPs = noGPsPerPolOrder[i];

            // Create a Gauss quadrature rule
            IGAGaussQuadratureOnTriangle* theIGAGaussQuadratureOnTriangle = new IGAGaussQuadratureOnTriangle(noGPs);

            // Initialize the integral value
            integral = 0.0;

            // Loop over all Gauss points and compute the integral value numerically
            for (int j = 0; j < noGPs; j++){
                gaussPointZeta1Value = *theIGAGaussQuadratureOnTriangle->getGaussPoint(j);
                gaussPointZeta2Value = *(theIGAGaussQuadratureOnTriangle->getGaussPoint(j) + 1);
                gaussPointWeight = theIGAGaussQuadratureOnTriangle->getGaussWeight(j);
                integral += log(gaussPointZeta1Value + 1)*sin(M_PI*gaussPointZeta2Value/2.0)*gaussPointWeight;
            }

            // Verify the solution against the analytical one
            CPPUNIT_ASSERT(fabs(integralAnalytical - integral) <= tolWeighted*pow(10,-i + 2));

            // Delete pointers
            delete theIGAGaussQuadratureOnTriangle;
        }
    }

    /***********************************************************************************************
     * \brief Solve the integral of function ln(zeta1 + 1)*sin(pi*zeta2/2) in the canonical triangle using the degenerated quadrilateral
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnCanonicalTriangleUsingTheDegeneratedQuadrilateral() {
        // Initialize auxiliary values
        int noGPs;
        double integral;
        double gaussPointZeta1Value;
        double gaussPointZeta2Value;
        double gaussPointWeight;
        double integralAnalytical = 0.050909793169365;
        double tolStrickened;

        // Loop over all quadratures
        for (int i = 0; i < 10; i++) {
            // Number of Gauss points for each quadrature
            noGPs = pow(i + 1,2.0);

            // Create a Gauss quadrature rule
            IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral* theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral =
                    new IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral(noGPs);

            // Initialize the integral value
            integral = 0.0;

            // Loop over all Gauss points and compute the integral value numerically
            for (int j = 0; j < noGPs; j++){
                gaussPointZeta1Value = *theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral->getGaussPoint(j);
                gaussPointZeta2Value = *(theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral->getGaussPoint(j) + 1);
                gaussPointWeight = theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral->getGaussWeight(j);
                integral += log(gaussPointZeta1Value + 1)*sin(M_PI*gaussPointZeta2Value/2.0)*gaussPointWeight;
            }

            // Verify the solution against the analytical one
            tolStrickened = tolBasic*pow(10.0,- 3.0*(i + 1)/2.0 + (i == 0)*1.0);
            CPPUNIT_ASSERT(fabs(integralAnalytical - integral) <= tolStrickened*(tolStrickened > tolAbs) + tolAbs*(tolStrickened <= tolAbs));

            // Delete pointers
            delete theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral;
        }
    }

    /***********************************************************************************************
     * \brief Tests the class IGAGaussQuadratureOnBiunitInterval for memory leakage
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnBiunitInterval4Leakage() {
        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            for(int i = 0; i < noGPsBiunitInterval; i++){
                // Create a Gauss quadrature rule
                IGAGaussQuadratureOnBiunitInterval* theIGAGaussQuadratureOnBiunitInterval = new IGAGaussQuadratureOnBiunitInterval(i + 1);

                // Delete pointers
                delete theIGAGaussQuadratureOnBiunitInterval;
            }
    }

    /***********************************************************************************************
     * \brief Tests the class IGAGaussQuadratureOnBiunitInterval for memory leakage
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnBiunitQuadrilateral4Leakage() {
        // Initialize auxiliary values
        int noGPs;

        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            for(int i = 0; i < 10; i++){
                // Create a Gauss quadrature rule
                noGPs = pow(i + 1,2.0);
                IGAGaussQuadratureOnBiunitQuadrilateral* theIGAGaussQuadratureOnBiunitQuadrilateral = new IGAGaussQuadratureOnBiunitQuadrilateral(noGPs);

                // Delete pointers
                delete theIGAGaussQuadratureOnBiunitQuadrilateral;
            }
    }

    /***********************************************************************************************
     * \brief Tests the class IGAGaussQuadratureOnTriangle for memory leakage
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnTriangle4Leakage() {
        // Initialize auxiliary arrays
        int noGPs;
        int noGPsPerPolOrder[8] = {1, 3, 4, 6, 7, 12, 13, 16};

        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            for (int i = 0; i < 8; i++) {
                // Number of Gauss points for each quadrature
                noGPs = noGPsPerPolOrder[i];

                // Create a Gauss quadrature rule
                IGAGaussQuadratureOnTriangle* theIGAGaussQuadratureOnTriangle = new IGAGaussQuadratureOnTriangle(noGPs);

                // Delete pointers
                delete theIGAGaussQuadratureOnTriangle;
            }
    }

    /***********************************************************************************************
     * \brief Tests the class IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral for memory leakage
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral4Leakage() {
        // Initialize auxiliary values
        int noGPs;

        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            for (int i = 0; i < 10; i++) {
                // Number of Gauss points for each quadrature
                noGPs = pow(i + 1,2.0);

                // Create a Gauss quadrature rule
                IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral* theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral =
                        new IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral(noGPs);

                // Delete pointers
                delete theIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral;
            }
    }

    CPPUNIT_TEST_SUITE( TestFEMMath );
    // Make the tests
    CPPUNIT_TEST(testIGAGaussQuadratureOnBiunitInterval);
    CPPUNIT_TEST(testIGAGaussQuadratureOnBiunitQuadrilateral);
    CPPUNIT_TEST(testIGAGaussQuadratureOnCanonicalTriangleUsingTheSymmetricRule);
    CPPUNIT_TEST(testIGAGaussQuadratureOnCanonicalTriangleUsingTheDegeneratedQuadrilateral);

    // Make the tests for memory leakage
    // CPPUNIT_TEST(testIGAGaussQuadratureOnBiunitInterval4Leakage);
    // CPPUNIT_TEST(testIGAGaussQuadratureOnBiunitQuadrilateral4Leakage);
    // CPPUNIT_TEST(testIGAGaussQuadratureOnTriangle4Leakage);
    // CPPUNIT_TEST(testIGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral4Leakage);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestFEMMath);
