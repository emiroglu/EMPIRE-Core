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

#include <fstream>
#include <vector>
#include <cstdlib>
#include <map>
#include <vector>
#include <assert.h>
#include <typeinfo>
#include <iostream>
#include <cmath>
#include "AuxiliaryParameters.h"

namespace EMPIRE {
namespace MathLibrary {

/***********************************************************************************************
 * \brief Compute mass matrix of a triangle element
 * \param[in] triangle the triangle
 * \param[in] numGaussPoints number of Gauss points used in the Gauss quadrature
 * \param[in] dual whether dual or not
 * \param[out] mass matrix (3x3)
 * \author Tianyang Wang
 ***********/
void computeMassMatrixOfTrianlge(const double *triangle, int numGaussPoints, bool dual,
        double *massMatrix);

/***********************************************************************************************
 * \brief Compute mass matrix of a quad element
 * \param[in] quad the quad
 * \param[in] numGaussPoints number of Gauss points used in the Gauss quadrature
 * \param[in] dual whether dual or not
 * \param[out] mass matrix (4x4)
 * \author Tianyang Wang
 ***********/
void computeMassMatrixOfQuad(const double *quad, int numGaussPoints, bool dual, double *massMatrix);

/***********************************************************************************************
 * \brief Compute the shape function value by local coordinates in a quadrilateral
 * \param[in] xi_eta local coordinates
 * \param[out] shapeFuncValues shape function values of xi_eta
 * \author Tianyang Wang
 ***********/
void computeShapeFuncOfQuad(const double *xi_eta, double *shapeFuncValues);

/***********************************************************************************************
 * \brief Compute the determinant of Jocobian matrix by local coordinates in a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] xi_eta local coordinates
 * \return determinant of Jocobian matrix
 * \author Tianyang Wang
 ***********/
double computeDetJOfQuad(const double *quad, const double *xi_eta);

/***********************************************************************************************
 * \brief Compute global coordinates of a point in a triangle
 * \param[in] triangle the triangle
 * \param[in] localCoor local coordinates of the point
 * \param[out] globalCoor global coordinates of the point
 * \author Tianyang Wang
 ***********/
void computeGlobalCoorInTriangle(const double *triangle, const double *localCoor,
        double *globalCoor);

/***********************************************************************************************
 * \brief Compute global coordinates of a point in a quad
 * \param[in] quad the quad
 * \param[in] localCoor local coordinates of the point
 * \param[out] globalCoor global coordinates of the point
 * \author Tianyang Wang
 ***********/
void computeGlobalCoorInQuad(const double *quad, const double *localCoor, double *globalCoor);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle
 * \param[in] triangle the triangle
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[in] point the point
 * \param[out] localCoor local coordinates of the point
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Tianyang Wang
 ***********/
bool computeLocalCoorInTriangle(const double *triangle, int planeToProject, const double *point,
        double *localCoor);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[in] point the point
 * \param[out] localCoor local coordinates of the point
 * \return a boolean saying whether the point is inside the quadrilateral or not (true of false)
 * \author Tianyang Wang
 ***********/
bool computeLocalCoorInQuad(const double *quad, int planeToProject, const double *point,
        double *localCoor);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle in a 2D space
 * \param[in] _coordsTriangle, coordinates of the triangle. double[6].
 * \param[in] _coordsNode, coordinates of the point. double[2].
 * \param[out] _localCoords, local coordinates of the point. double[3]
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInTriangle(const double *_coordsTri, const double *_coordsNode,
        double* _localCoords);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadriliteral in a 2D space by solving a nonlinear
 *        system using the Newton-Raphson scheme
 * \param[in] _coordsQuad Coordinates of the quadriliteral. double[8].
 * \param[in] _coordsNode Coordinates of the point. double[2].
 * \param[out] _localCoords local coordinates of the point. double[2]
 * \return a boolean saying whether the point is inside the quadriliteral or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInQuad(const double *_coordsQuad, const double *_coordsNode,
        double* _localCoords);

/***********************************************************************************************
 * \brief Computes the values of the low order shape functions (linear for triangle and bilinear
 *        for the quadrilateral)
 * \param[in] _nNodes The number of nodes in the element level
 * \param[in] _coords The coordinates of the point where to evaluate the shape functions
 * \param[in/out] _shapeFuncs The evaluated shape functions
 * \author Chenshen Wu
 ***********/
void computeLowOrderShapeFunc(int _nNodes, const double* _coords, double* _shapeFuncs);

/***********************************************************************************************
 * \brief Computes the value of a data field in the interior of an element
 * \param[in] _nNodes The number of nodes of the element
 * \param[in] _nValue Takes value 1 for a scalar, or the values 2-3 for a vector
 * \param[in] _values The values on the nodes of the element
 * \param[in] _shapeFuncs The values of the shape functions at the interior of the element
 * \param[in/out] _returnValue The resulting linear combination of the nodal values
 * \author Chenshen Wu
 ***********/
void computeLinearCombination(int _nNodes, int _nValue, const double * _values,
        const double *_shapeFuncs, double* _returnValue);

/********//**
 * \brief Class IntegrandFunction is the mother class of all integrand functions
 * \author Tianyang Wang
 ***********/
class IntegrandFunction {
public:
    IntegrandFunction() {
    }
    virtual ~IntegrandFunction() {
    }
    /***********************************************************************************************
     * \brief Compute the function value on the Gauss point
     * \param[in] gaussPoint x,y,z coordinates of the Gauss point
     * \return the function value
     ***********/
    virtual double operator()(double *gaussPoint) = 0;
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on triangle
 * \author Tianyang Wang
 ***********/
class GaussQuadratureOnTriangle {
public:
    /***********************************************************************************************
     * \brief constructor of the class GaussQuadratureOnTriangle
     * \param[in] _triangle triangle to integrate on
     * \param[in] _nGaussPoints number of gauss points
     * \author Tianyang Wang
     ***********/
    GaussQuadratureOnTriangle(double *_triangle, int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor of the class.
     * \author Tianyang Wang
     ***********/
    virtual ~GaussQuadratureOnTriangle();

    /***********************************************************************************************
     * \brief To define the integrand function
     * \param[in] _integrandFunc Object of the class integrand functions
     * \author Tianyang Wang
     ***********/
    void setIntegrandFunc(IntegrandFunction *_integrandFunc);

    /***********************************************************************************************
     * \brief To Compute integral on the triangle
     * \author Tianyang Wang
     ***********/
    double computeIntegral();

    const int numGaussPoints;
    double *gaussPointsGlobal;
private:
    const double *triangle;
    const double *gaussPointsLocal;
    const double *weights;
    IntegrandFunction *integrandFunc;
    double area;
};

/********//**
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Tianyang Wang
 ***********/
class GaussQuadratureOnQuad {
public:
    /***********************************************************************************************
     * \brief constructor of the class GaussQuadratureOnQuad
     * \param[in] _quad quad to integrate on
     * \param[in] _nGaussPoints number of gauss points
     * \author Tianyang Wang
     ***********/
    GaussQuadratureOnQuad(double *_quad, int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~GaussQuadratureOnQuad();

    /***********************************************************************************************
     * \brief To define the integrand function
     * \param[in] _integrandFunc Object of the class integrand functions
     * \author Tianyang Wang
     ***********/
    void setIntegrandFunc(IntegrandFunction *_integrandFunc);

    /***********************************************************************************************
     * \brief To Compute integral on the quad
     * \author Tianyang Wang
     ***********/
    double computeIntegral();

    const int numGaussPoints;
    double *gaussPointsGlobal;
private:
    const double *quad;
    const double *gaussPointsLocal;
    const double *weights;
    double *detJ;
    IntegrandFunction *integrandFunc;
};

// IGA Classes
// %%%%%%%%%%%%%%%%%%%%%%%%%

/********//**
 * \brief Class GaussQuadrature base class for general quadrature
 * \author Andreas Apostolatos
 ***********/
class IGAGaussQuadrature {

private:
    // The number of Gauss Points
    int numGaussPoints;

    // Number of dimensions
    int numDimensions;

    // Array containing the Gauss Point locations in the quadrature space
    const double *gaussPoints;

    // Array containing the quadrature weights
    const double *weights;

    /// Constructor, destructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    IGAGaussQuadrature(int _numGaussPoints, int _numDimensions) :
            numGaussPoints(_numGaussPoints), numDimensions(_numDimensions) {
    }

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~IGAGaussQuadrature() {
    }

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Returns the number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    double getNumGaussPoints() {
        return numGaussPoints;
    }

    /***********************************************************************************************
     * \brief Returns the coordinates of the Gauss point
     * \param[in] _index Index to the Gauss point
     * \author Andreas Apostolatos
     ***********/
    const double* getGaussPoint(int _index) {
        return &gaussPoints[_index * numDimensions];
    }

    /***********************************************************************************************
     * \brief Returns the weight of the Gauss point
     * \param[in] _index Index to the Gauss weight
     * \author Andreas Apostolatos
     ***********/
    double getGaussWeight(int _index) {
        return weights[_index];
    }

    /***********************************************************************************************
     * \brief Sets the coordinates of the Gauss points
     * \param[in] _values Pointer to the array with the values
     * \author Andreas Apostolatos
     ***********/
    void setGaussPoints(const double* _values) {
        gaussPoints = _values;
    }

    /***********************************************************************************************
     * \brief Sets the Gauss weights
     * \param[in] _values Pointer to the array with the values
     * \author Andreas Apostolatos
     ***********/
    void setGaussWeights(const double* _values) {
        weights = _values;
    }
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on triangle using a symmetric rule
 * \author Andreas Apostolatos
 ***********/
class IGAGaussQuadratureOnTriangle: public IGAGaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    IGAGaussQuadratureOnTriangle(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~IGAGaussQuadratureOnTriangle() {
    }
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on canonical triangle using the degenerated quadrilateral
 * \author Andreas Apostolatos
 ***********/
class IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral: public IGAGaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~IGAGaussQuadratureOnTriangleUsingDegeneratedQuadrilateral() {
    }
};

/**********
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Andreas Apostolatos, Altug Emiroglu
 ***********/
class IGAGaussQuadratureOnBiunitInterval: public IGAGaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    IGAGaussQuadratureOnBiunitInterval(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~IGAGaussQuadratureOnBiunitInterval() {
    }
};

/**********
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Andreas Apostolatos
 ***********/
class IGAGaussQuadratureOnBiunitQuadrilateral: public IGAGaussQuadrature {
public:
    /***********************************************************************************************
     * \brief Constructor
     * param[in] _numGaussPoints, number of Gauss points
     * \author Andreas Apostolatos
     ***********/
    IGAGaussQuadratureOnBiunitQuadrilateral(int _numGaussPoints);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    virtual ~IGAGaussQuadratureOnBiunitQuadrilateral() {
    }
};

}
}
