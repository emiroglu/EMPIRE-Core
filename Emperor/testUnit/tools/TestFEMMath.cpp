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
    int noGPs;
    double Tol;
public:
    /***********************************************************************************************
     * \brief Set up the Gauss quadrature
     * \author Andreas Apostolatos
     ***********/
    void setUp() {
        // Number of Gauss points that can be withdrawn from the function
        noGPs = 50;

        // Basic tolerance
        Tol = 1e-1;
    }
    /***********************************************************************************************
     * \brief Delete test vectors and matrices
     * \author Andreas Apostolatos
     ***********/
    void tearDown() {
    }
    /***********************************************************************************************
     * \brief Solve the integral of function xlnx in [1 5] interval analytically for all Gauss Point quadratures
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

        for(int i = 0; i < noGPs; i++){
            // Create a Gauss quadrature rule
            IGAGaussQuadratureOnBiunitInterval* theIGAGaussQuadratureOnBiunitInterval = new IGAGaussQuadratureOnBiunitInterval(i + 1);

            // Initialize the integral value
            integral = 0.0;

            // Loop over all Gauss points and compute the integral value numerically
            for(int j = 0; j < i + 1; j++){
                gaussPointValueOnParentSpace = theIGAGaussQuadratureOnBiunitInterval->gaussPoints[j];
                gaussPointWeight = theIGAGaussQuadratureOnBiunitInterval->weights[j];
                gaussPointValue = 2*gaussPointValueOnParentSpace + 3;
                integral += gaussPointValue*log(gaussPointValue)*2*gaussPointWeight;
            }

            // Verify the solution against the analytical one
            if(i <= 15)
                TolWeighted = Tol*pow(10.0, - i + 1);
            else
                TolWeighted = 1e-14;
            CPPUNIT_ASSERT(fabs(integralAnalytical - integral) <= TolWeighted);

            // Delete pointers
            delete theIGAGaussQuadratureOnBiunitInterval;
        }
    }

    /***********************************************************************************************
     * \brief Tests the class IGAGaussQuadratureOnBiunitInterval for memory leakage
     * \author Andreas Apostolatos
     ***********/
    void testIGAGaussQuadratureOnBiunitInterval4Leakage() {
        // Create and destroy the same objects iteratively
        for(int i = 0; i < 1e10; i++)
            for(int i = 0; i < noGPs; i++){
                // Create a Gauss quadrature rule
                IGAGaussQuadratureOnBiunitInterval* theIGAGaussQuadratureOnBiunitInterval = new IGAGaussQuadratureOnBiunitInterval(i + 1);

                // Delete pointers
                delete theIGAGaussQuadratureOnBiunitInterval;
            }
    }

    CPPUNIT_TEST_SUITE( TestFEMMath );
    // Make the tests
    CPPUNIT_TEST(testIGAGaussQuadratureOnBiunitInterval);

    // Make the tests for memory leakage
    // CPPUNIT_TEST(testIGAGaussQuadratureOnBiunitInterval4Leakage);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestFEMMath);
