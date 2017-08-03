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
#ifndef WEAKIGADIRICHLETCURVECONDITION_H
#define WEAKIGADIRICHLETCURVECONDITION_H

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"
#include "AbstractCondition.h"

namespace EMPIRE {

class Message;
class IGAPatchCurve;
class IGAPatchSurface;

/********//**
 * \brief Class WeakIGADirichletCurveCondition is the class that holds and applies the weak Dirichlet conditions on a curve for IGA
 ***********/
class WeakIGADirichletCurveCondition: public AbstractCondition {

protected:

    /// The curve defined on a patch parameter space to apply Dirichlet conditions
    IGAPatchCurve* dirichletCurve;

    /// Flag on whether the curve is a trimming curve or an arbitrary curve
    bool isTrimmingCurve;

    /// If the GP data for the Dirichlet condition is initialized
    bool isGPDataInitialized;

    /// Number of GPs on the trimming curve
    int curveNumGP;

    /// GP weights on the trimming curve
    double* curveGPWeights;

    /// Jacobian products at the GPs on the trimming curve
    double* curveGPJacobianProducts;

    /// The index of the patch in the sending order
    int patchIndex;

    /// The index of the patch boundary loop in the sending order
    int patchBLIndex;

    /// The index of the patch boundary loop trimming curve in the sending order
    int patchBLTrCurveIndex;

    /// Curve GP coordinates in the parametric space of the patch
    double* curveGPs;

    /// Curve tangents in the parametric space of the patch
    double* curveGPTangents;

public:

    /***********************************************************************************************
     * brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _p The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGADirichletCurveCondition(int _ID,
                                   int _patchIndex, int _p, int _uNoKnots, double* _uKnotVector, int _uNoControlPoints, double* _controlPointNet);

    /***********************************************************************************************
     * \brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _dirichletCurve The curve to apply Dirichlet conditions
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGADirichletCurveCondition(int _ID,
                                   int _patchIndex, IGAPatchCurve* _dirichletCurve);

    /***********************************************************************************************
     * \brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _patchBLIndex The index of the patch boundary loop in the EMPIRE data structure
     * \param[in] _patchBLTrCurveIndex The index of the patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGADirichletCurveCondition(int _ID,
                                   int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex);

    /***********************************************************************************************
     * \brief Set the GP data for the Dirichlet condition
     * \param[in] _curveNumGP The total number of GPs on the trimming curve
     * \param[in] _curveGPs The parametric coordinates of the GPs on the trimming curve in the patch parameter space
     * \param[in] _curveGPWeights The GP weights
     * \param[in] _curveGPTangents The tangent to the trimming curve vector in the physical space
     * \param[in] _curveGPJacobianProducts The Jacobian products for the transformation of the integrals
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void addWeakDirichletCurveConditionGPData(int _curveNumGP,
                                              double* _curveGPs, double* _curveGPWeights, double* _curveGPTangents,
                                              double* _curveGPJacobianProducts);

    /***********************************************************************************************
     * \brief Create the GP data for the Dirichlet condition
     * \param[in] _patch The patch on which the trimming curve exists
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void createGPData(IGAPatchSurface* _patch);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    ~WeakIGADirichletCurveCondition();

    /***********************************************************************************************
     * \brief get patchIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getPatchIndex() const {
        return patchIndex;
    }

    /***********************************************************************************************
     * \brief get curveNumGP
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getCurveNumGP() const {
        return curveNumGP;
    }

    /***********************************************************************************************
     * \brief get curveGPs
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getCurveGPs() const {
        return curveGPs;
    }

    /***********************************************************************************************
     * \brief get curveGPWeights
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getCurveGPWeights() const {
        return curveGPWeights;
    }

    /***********************************************************************************************
     * \brief get curveGPTangents
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getCurveGPTangents() const {
        return curveGPTangents;
    }

    /***********************************************************************************************
     * \brief get curveGPJacobianProducts
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getCurveGPJacobianProducts() const {
        return curveGPJacobianProducts;
    }

    /***********************************************************************************************
     * \brief get getCurveGPData
     * \param[out] _curveGP         The parametric coordinates of the GP on the curve in the patch parameter space
     * \param[out] _curveGPWeight   The GP weight
     * \param[out] _curveGPTangent  The tangent to the trimming curve vector in the physical space
     * \param[out] _curveGPJacobianProduct The Jacobian product for the transformation of the integrals
     * \param[in]  _iGP The requested GP counter
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void getCurveGPData(double* _curveGP, double& _curveGPWeight, double* _curveGPTangent, double& _curveGPJacobianProduct, int _iGP);

    /***********************************************************************************************
     * \brief get getCurveAllGPData
     * \param[out] _curveGPs The parametric coordinates of the GPs on the trimming curve in the patch parameter space
     * \param[out] _curveGPWeights The GP weights
     * \param[out] _curveGPTangents The tangent to the trimming curve vector in the physical space
     * \param[out] _curveGPJacobianProducts The Jacobian products for the transformation of the integrals
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void getCurveAllGPData(double* _curveGPs, double* _curveGPWeights, double* _curveGPTangents, double* _curveGPJacobianProducts);

};

} /* namespace EMPIRE */

#endif /* WEAKIGADIRICHLETCURVECONDITION_H */
