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
 * \date 20/11/2013
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
    /// The basis functions of the 2D NURBS patch
    BSplineBasis1D* IGABasis;

    /// Number of Control Points in u-direction
    int uNoControlPoints;

    /// The set of the Control Points of the patch
    std::vector<IGAControlPoint> ControlPointNet;


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
    inline const IGAControlPoint& getControlPoint(int i) const {
        return ControlPointNet.at(i);
    }
};

} /* namespace EMPIRE */
#endif /* IGAPATCHCURVE_H_ */
