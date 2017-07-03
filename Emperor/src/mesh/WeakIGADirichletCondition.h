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
#ifndef WEAKIGADIRICHLETCONDITION_H_
#define WEAKIGADIRICHLETCONDITION_H_

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"
#include "AbstractCondition.h"

namespace EMPIRE {

class Message;
/********//**
 * \brief Class WeakIGADirichletCondition is the class that holds and applies the weak continuity conditions for IGA
 ***********/
class WeakIGADirichletCondition: public AbstractCondition {

protected:

    /// ID of the condition
    int ID;

    /// If the GP for coupling are provided
    bool isGPProvided;

    /// Number of GPs on the coupled trimming curve
    int trCurveNumGP;

    /// GP weights on the coupled trimming curve
    double* trCurveGPWeights;

    /// Jacobian products at the GPs on the coupled trimming curve
    double* trCurveGPJacobianProducts;

    /// The index of the patch in the sending order
    int patchIndex;

    /// The index of the patch boundary loop in the sending order
    int patchBLIndex;

    /// The index of the patch boundary loop trimming curve in the sending order
    int patchBLTrCurveIndex;

    /// GP parametric coordinates on the trimming curve of the patch
    double* trCurveGPs;

    /// Trimming curve tangents in the parametric space of the patch
    double* trCurveGPTangents;

public:
    /***********************************************************************************************
     * \brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _patchBLIndex The index of the patch boundary loop in the EMPIRE data structure
     * \param[in] _patchBLTrCurveIndex The index of the patch trimming curve in the current boundary loop in the EMPIRE data structure
     * \param[in] _isGPprovided Flag if the GP data is provided
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGADirichletCondition(int _ID,
                    int _patchIndex, int _patchBLIndex, int _patchBLTrCurveIndex,
                    int _isGPProvided);

    /***********************************************************************************************
     * \brief Set the GP data for the coupling condition
     * \param[in] _trCurveNumGP The total number of GPs on the trimming curve
     * \param[in] _trCurveGPs The parametric coordinates of the GPs on the trimming curve in the patch parameter space
     * \param[in] _trCurveGPWeights The GP weights
     * \param[in] _trCurveGPTangents The tangent to the trimming curve vector in the parameter space of the patch
     * \param[in] _trCurveGPJacobianProducts The Jacobian products for the transformation of the integrals
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void addWeakDirichletConditionGPData(int _trCurveNumGP,
               double* _trCurveGPs, double* _trCurveGPWeights, double* _trCurveGPTangents,
               double* _trCurveGPJacobianProducts);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    ~WeakIGADirichletCondition();

    /***********************************************************************************************
     * \brief get patchIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getPatchIndex() const {
        return patchIndex;
    }

    /***********************************************************************************************
     * \brief get trCurveNumGP
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getTrCurveNumGP() const {
        return trCurveNumGP;
    }

    /***********************************************************************************************
     * \brief get trCurveGPs
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveGPs() const {
        return trCurveGPs;
    }

    /***********************************************************************************************
     * \brief get trCurveGPWeights
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveGPWeights() const {
        return trCurveGPWeights;
    }

    /***********************************************************************************************
     * \brief get trCurveGPTangents
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveGPTangents() const {
        return trCurveGPTangents;
    }

    /***********************************************************************************************
     * \brief get trCurveGPJacobianProducts
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getTrCurveGPJacobianProducts() const {
        return trCurveGPJacobianProducts;
    }

};

} /* namespace EMPIRE */

#endif /* WEAKIGADIRICHLETCONDITION_H_ */
