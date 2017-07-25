/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Munich
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
 * \file IGAPatchCurve.h
 * This file holds the class IGAPatchCurve.h
 * \date 20/11/2014
 **************************************************************************************************/
#ifndef IGAPATCHCURVE_H_
#define IGAPATCHCURVE_H_

// Inclusion of user defined libraries
#include "Message.h"
#include "AbstractMesh.h"
#include "BSplineBasis1D.h"
#include "NurbsBasis1D.h"
#include "IGAControlPoint.h"

namespace EMPIRE {

class IGAPatchCurve {
protected:
    /// The basis functions of the 1D NURBS patch
    BSplineBasis1D* IGABasis;

    /// Number of Control Points in u-direction
    int uNoControlPoints;

    /// The set of the Control Points of the patch
    std::vector<IGAControlPoint> ControlPointNet;

    // List of points making up the linearized version of the curve (u,v)
    std::vector<double> polyline;

    // List of points making up the linearized version of the curve (uTilde)
    std::vector<double> polylineKnots;

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _IDBasis The id of the underlying basis to the IGA 1D patch
     * \param[in] _pDegree The polynomial degree of the IGA 1D patch in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 1D patch in the u-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 2D NURBS patch
     * \author Fabien Pean
     ***********/
    IGAPatchCurve(int, int, int, double*, int, double*);

    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean
     ***********/
	virtual ~IGAPatchCurve();

    /// Basis related functions
public:
    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates in 3D of a point on a B-Spline surface given the basis functions and the knot span index at the parametric location
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose basis functions are _locBasisFunctions
     * \param[in] _knotSpanIndex The knot span index of the knot span where the point lies
     * \param[in] _locBasisFunctions The non-zero basis functions at the parametric location where the point is
     ***********/
    void computeCartesianCoordinates(double* _cartesianCoordinates, int _knotSpanIndex, double* _locBasisFunctions);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates in 2D of a point on a B-Spline surface given the basis functions and the knot span index at the parametric location
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose basis functions are _locBasisFunctions
     * \param[in] _knotSpanIndex The knot span index of the knot span where the point lies
     * \param[in] _locBasisFunctions The non-zero basis functions at the parametric location where the point is
     ***********/
    void computeCartesianCoordinates(double* _cartesianCoordinates, double* _locBasisFunctions, int _knotSpanIndex);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters are known
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose surface parameters are _uPrm
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \author Fabien Pean, Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double, int) const;

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \author Fabien Pean, Chenshen Wu
     ***********/
    void computeCartesianCoordinates(double*, double) const;

    /***********************************************************************************************
     * \brief Returns the base vector components and their derivatives
     * \param[in/out] _baseVectorAndDerivatives The base vector and its derivatives
     * \param[in] _knotSpanIndex The knot span index where the given parameter lies in
     * \param[in] _basisFunctionsAndDerivatives The basis functions and their derivatives
     * \param[in] _derivOrder The derivative order of the base vector
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectorAndDerivatives(double* _baseVectorAndDerivatives, int _knotSpanIndex, double* _basisFunctionsAndDerivatives, int _derivOrder) const;

    /// Geometric operation functions
public:
    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates and the Curve Parameters of intersections with a parameter line in an ascending order with respect to their curve parameters
     * \param[in/out] _uvSurface Coordinates of the intersection in the patch parameter space
     * \param[in/out] _uTilde Coordinate of the intersection in the curve parameter space
     * \param[in] _dir Flag to indicate if knot is along u(True) or v(False)
     * \param[in] _knot Knot to intersect
     * \return The flag on whether or not the iterations have converged
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    bool computeIntersectionsWithKnotBisection(std::vector<double>& _uTilde, std::vector<double>& _uvSurface, unsigned int _dir, double _knot) const;

    /***********************************************************************************************
     * \brief Returns the curve parameters of the intersections with a parameter line in an ascending order with respect to their curve parameters
     * This functions expects either a sorted or a zero sized _uTilde
     * \param[in/out] _uTilde Coordinate of the intersection in the curve parameter space
     * \param[in] _dir Flag to indicate if knot is along u(True) or v(False)
     * \param[in] _knot Knot to intersect
     * \return The flag on whether or not the iterations have converged
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    bool computeIntersectionsWithKnotBisection(std::vector<double>& _uTilde, unsigned int _dir, double _knot) const;

    /***********************************************************************************************
     * \brief Solves the bisection problem for computing the intersection of a parameter line with a trimming curve
     * \param[in/out] _uvP Coordinates of the intersection in the patch parameter space
     * \param[in/out] _uTildeP Coordinate of the intersection in the curve parameter space
     * \param[in] _uTildeP1 Curve parameter of the first point to use for bracketing
     * \param[in] _uTildeP2 Curve parameter of the second point to use for bracketing
     * \param[in] _dir Direction of the parameter line
     * \param[in] _knot Knot to intersect
     * \return The flag on whether or not the iterations have converged
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    bool solveIntersectionWithKnotBisection(double* _uvP, double& _uTildeP, double _uTildeP1, double _uTildeP2, unsigned int _dir, double _knot) const;

    /***********************************************************************************************
     * \brief Computes the point projection on a 2D B-Spline curve. The point and the curve are assumed to be 2D meaning that the 3rd coordinate is not taken in to consideration
     * \param[in/out] _uPrm The initial guess for the Newton-Raphson algorithm which is overwritten by the parametric coordinates of the projected Point after the function call
     * \param[in/out] _P The point to be projected, after the function the point is overwritten by the found projection
     * \param[in] _noCoord The number of coordinates of the given point _p
     * \param[in] _maxIt Maximum number of iterations for the Newton-Raphson algorithm
     * \param[in] _tol Convergence tolerance for the Newton-Raphson iterations
     * \param[out] Boolean on the convergence of the Newton-Raphson algorithm
     * \author Andreas Apostolatos, Altug Emiroglu
     */
    bool computePointProjectionOn2DCurve(double& _uPrm, double* _P, int _noCoord, int _maxIt=MAX_NUM_ITERATIONS_NEWTONRAPHSON, double _tol=TOL_CONVERGENCE_NEWTONRAPHSON);

    /***********************************************************************************************
     * \brief Create a linear approximation for the curve (Implementation migrated from IGAPatchSurfaceTrimming object)
     * \param[in] _type (0: NCPxP, 1: Greville, 2: Combined) Type of linearization
     * \param[in] _dir Direction of the linearization. True starts from beginning, False starts from the end
     * \author Altug Emiroglu
     ***********/
    void linearize(int _type = 2, bool _dir = true);

    void linearizeUsingGreville(bool _dir = true);

    void linearizeUsingNCPxP(bool _dir = true);

    void linearizeCombined(bool _dir = true);

    /// get functions
public:
    /***********************************************************************************************
     * \brief Get the underlying IsoGeometric basis
     * \author Fabien Pean
     ***********/
    inline const BSplineBasis1D* getIGABasis() const {
        return IGABasis;
    }
    /***********************************************************************************************
     * \brief Get the number of the Control Points
     * \author Fabien Pean
     ***********/
    inline int getNoControlPoints() const {
        return uNoControlPoints;
    }
    /***********************************************************************************************
     * \brief Get the Control Points
     * \author Fabien Pean
     ***********/
    inline const std::vector<IGAControlPoint>& getControlPointNet() const {
        return ControlPointNet;
    }
    /***********************************************************************************************
     * \brief Get a specific Control Point
     * \author Fabien Pean
     ***********/
    inline const IGAControlPoint& getControlPoint(const unsigned int i) const {
        return ControlPointNet.at(i);
    }
    inline const IGAControlPoint& operator[](const unsigned int i) {
    	return ControlPointNet.at(i);
    }
    /***********************************************************************************************
     * \brief Find knot span on u direction
     * \author Fabien Pean
     ***********/
    inline int findKnotSpan(double _u) const {
        return getIGABasis()->findKnotSpan(_u);
    }

    /***********************************************************************************************
     * \brief Get the patch parameters of the linearization of the trimming curve
     * \author Altug Emiroglu
     ***********/
    inline const std::vector<double>* getPolyline() const {
        return &polyline;
    }

    /***********************************************************************************************
     * \brief Get the curve parameters of the linearization of the trimming curve
     * \author Altug Emiroglu
     ***********/
    inline const std::vector<double>* getPolylineKnots() const {
        return &polylineKnots;
    }

    /// set functions
public:

    /***********************************************************************************************
     * \brief Add a vertex to the linearization of the trimming curve
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    inline void addPolylineVertex(double _knot, double _u, double _v){
        polylineKnots.push_back(_knot);
        polyline.push_back(_u);
        polyline.push_back(_v);
    }

    /// The maximum number of iterations for the computation of the intersection with a line
    static int MAX_NUM_ITERATIONS;
    /// The tolerance for the iterations for the computation of the intersection with a line
    static double TOL_CONVERGENCE;
    /// The maximum number of iterations for the projection of a point on a curve using the Newton-Raphson algorithm
    static int MAX_NUM_ITERATIONS_NEWTONRAPHSON;
    /// The tolerance for the iterations for the projection of a point on a curve using the Newton-Raphson algorithm
    static double TOL_CONVERGENCE_NEWTONRAPHSON;

    /// The tolerance in the curve parameter space for the linearization of the trimming loops
    static double TOL_LINEARIZATION;

};

} /* namespace EMPIRE */
#endif /* IGAPATCHCURVE_H_ */
