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
#ifndef WEAKIGADIRICHLETSURFACECONDITION_H
#define WEAKIGADIRICHLETSURFACECONDITION_H

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"
#include "AbstractCondition.h"

namespace EMPIRE {

class Message;
class IGAPatchSurfaceTrimmingLoop;
class IGAPatchSurface;

/********//**
 * \brief Class WeakIGADirichletSurfaceCondition is the class that holds and applies the weak Dirichlet conditions on a Surface for IGA
 ***********/
class WeakIGADirichletSurfaceCondition: public AbstractCondition {
private:
    /// Type definitions
    typedef std::pair<double,double> Point2D;
    typedef std::vector<Point2D> Polygon2D;
    typedef std::vector<Polygon2D> ListPolygon2D;
protected:

    /// The Surface defined on a patch parameter space to apply Dirichlet conditions
    IGAPatchSurfaceTrimmingLoop* conditionBoundaryLoop;

    /// Flag on whether the Surface is a trimming Surface or an arbitrary Surface
    bool isBoundaryLoop;

    /// If the GP data for the Dirichlet condition is initialized
    bool isGPDataInitialized;

    /// The surface polygons to apply Dirichlet condition
    ListPolygon2D conditionPolygons;

    /// Number of GPs on the trimming Surface
    int surfaceNumGP;

    /// Surface GP coordinates in the parametric space of the patch
    double* surfaceGPs;

    /// GP weights on the trimming Surface
    double* surfaceGPWeights;

    /// Jacobian products at the GPs on the trimming Surface
    double* surfaceGPJacobians;

    /// The index of the patch in the sending order
    int patchIndex;

    /// The index of the patch boundary loop in the sending order
    int patchBLIndex;

public:
    /***********************************************************************************************
     * \brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGADirichletSurfaceCondition(int _ID,
                                     int _patchIndex, IGAPatchSurfaceTrimmingLoop* _conditionBoundaryLoop);

    /***********************************************************************************************
     * \brief Constructor, initializing the class
     * \param[in] _ID ID of the condition
     * \param[in] _patchIndex The index of the patch in the EMPIRE data structure
     * \param[in] _patchBLIndex The index of the patch boundary loop in the EMPIRE data structure
     * \param[in] _patchBLTrSurfaceIndex The index of the patch trimming Surface in the current boundary loop in the EMPIRE data structure
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    WeakIGADirichletSurfaceCondition(int _ID,
                                     int _patchIndex, int _patchBLIndex);

    /***********************************************************************************************
     * \brief Create the GP data for the Dirichlet condition
     * \param[in] _patch The patch on which the boundary loop exists
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    void createGPData(const std::vector<IGAPatchSurface*>& _surfacePatches);

    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming curves of the given patch
     * This function does the same as in IGAMortarMapper
     * \param[in] _thePatch 	The patch for which trimming curves are used
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of trimming polygon
     * \author Altug Emiroglu
     ***********/
    void clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV);

    /***********************************************************************************************
     * \brief Clip the input polygon by the trimming window of the given trimming loop
     * This function does the same as in IGAMortarMapper
     * \param[in] _theTrimmingLoop 	The trimming loop for which trimming curves are used
     * \param[in] _polygonUV 	An input polygon defined in parametric (i.e. 2D) space
     * \param[out] _listPolygon	A set of polygons after application of trimming polygon
     * \author Altug Emiroglu
     ***********/
    void clipByCondition(const IGAPatchSurfaceTrimmingLoop* _theTrimmingLoop, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV);

    ListPolygon2D triangulatePolygon(const Polygon2D& _polygonUV);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    ~WeakIGADirichletSurfaceCondition();

    /***********************************************************************************************
     * \brief get patchIndex
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getPatchIndex() const {
        return patchIndex;
    }

    /***********************************************************************************************
     * \brief get surfaceNumGP
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    int getSurfaceNumGP() const {
        return surfaceNumGP;
    }

    /***********************************************************************************************
     * \brief get surfaceGPs
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getSurfaceGPs() const {
        return surfaceGPs;
    }

    /***********************************************************************************************
     * \brief get surfaceGPWeights
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getSurfaceGPWeights() const {
        return surfaceGPWeights;
    }

    /***********************************************************************************************
     * \brief get surfaceGPJacobians
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    double* getSurfaceGPJacobians() const {
        return surfaceGPJacobians;
    }

    /***********************************************************************************************
     * \brief get getConditionBoundaryLoop
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    IGAPatchSurfaceTrimmingLoop* getConditionBoundaryLoop() const {
        return conditionBoundaryLoop;
    }

    /***********************************************************************************************
     * \brief get getIsBoundaryLoop
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    bool getIsBoundaryLoop() const {
        return isBoundaryLoop;
    }

//    /***********************************************************************************************
//     * \brief get getSurfaceGPData
//     * \param[out] _surfaceGP         The parametric coordinates of the GP inside the boundary loop in the patch parameter space
//     * \param[out] _surfaceGPWeight   The GP weight
//     * \param[out] _surfaceGPJacobianProduct The Jacobian product for the transformation of the integrals
//     * \param[in]  _iGP The requested GP counter
//     * \author Andreas Apostolatos, Altug Emiroglu
//     ***********/
//    void getSurfaceGPData(double* _surfaceGP, double& _surfaceGPWeight, double& _surfaceGPJacobianProduct, int _iGP);

//    /***********************************************************************************************
//     * \brief get getSurfaceAllGPData
//     * \param[out] _SurfaceGPs The parametric coordinates of the GPs on the trimming Surface in the patch parameter space
//     * \param[out] _SurfaceGPWeights The GP weights
//     * \param[out] _surfaceGPJacobians The Jacobian products for the transformation of the integrals
//     * \author Andreas Apostolatos, Altug Emiroglu
//     ***********/
//    void getSurfaceAllGPData(double* _surfaceGPs, double* _surfaceGPWeights, double* _surfaceGPTangents, double* _surfaceGPJacobians);

};

} /* namespace EMPIRE */

#endif /* WEAKIGADIRICHLETSURFACECONDITION_H */
