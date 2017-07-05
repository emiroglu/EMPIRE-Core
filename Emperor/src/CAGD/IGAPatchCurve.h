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

    // List of points making up the linearized version of the curve (u,v,uTilde)
    std::vector<double> polyline;

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

    /// Geometric operation functions
public:
    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
     * \param[in/out] _uvSurface Coordinates of the intersection in the patch parameter space
     * \param[in/out] _uTilde Coordinate of the intersection in the curve parameter space
     * \param[in] _dir Flag to indicate if knot is along u(True) or v(False)
     * \param[in] _knot Knot to intersect
     * \return The flag on whether or not the iterations have converged
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    bool computeIntersectionsWithKnotBisection(std::vector<double>& _uvSurface , std::vector<double>& _uTilde, unsigned int _dir, double _knot) const;

    bool solveIntersectionWithKnotBisection(double* _uvP, double& _uTildeP, double _uTildeP1, double _uTildeP2, unsigned int _dir, double _knot) const;

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

    /// set functions
public:

    /***********************************************************************************************
     * \brief Add a vertex to the linearization of the trimming curve
     * \author Altug Emiroglu, Andreas Apostolatos
     ***********/
    inline void addPolylineVertex(double _u, double _v, double _knot){
        polyline.push_back(_u);
        polyline.push_back(_v);
        polyline.push_back(_knot);
    }

    /// The maximum number of iterations for the computation of the intersection with a line
    static int MAX_NUM_ITERATIONS;
    /// The tolerance for the iterations for the computation of the intersection with a line
    static double TOL_CONVERGENCE;

};

} /* namespace EMPIRE */
#endif /* IGAPATCHCURVE_H_ */
