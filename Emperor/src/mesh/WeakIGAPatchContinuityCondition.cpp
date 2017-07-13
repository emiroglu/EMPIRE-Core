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

// Inclusion of standard libraries
#include <iostream>
#include <assert.h>
#include <math.h>

// Inclusion of user defined libraries
#include "IGAPatchSurface.h"
#include "WeakIGAPatchContinuityCondition.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

WeakIGAPatchContinuityCondition::WeakIGAPatchContinuityCondition(int _ID, 
								 int _masterPatchIndex, int _masterPatchBLIndex, int _masterPatchBLTrCurveIndex,
								 int _slavePatchIndex, int _slavePatchBLIndex, int _slavePatchBLTrCurveIndex,
                                 int _isGPProvided) :
        AbstractCondition(_ID),
        masterPatchIndex(_masterPatchIndex), masterPatchBLIndex(_masterPatchBLIndex), masterPatchBLTrCurveIndex(_masterPatchBLTrCurveIndex),
        slavePatchIndex(_slavePatchIndex), slavePatchBLIndex(_slavePatchBLIndex), slavePatchBLTrCurveIndex(_slavePatchBLTrCurveIndex),
        isGPProvided(_isGPProvided)
        {

    type = EMPIRE_WeakIGAPatchContinuityCondition;

    isGPDataInitialized = false;
}

void WeakIGAPatchContinuityCondition::addWeakContinuityConditionGPData(int _trCurveNumGP,
							double* _trCurveMasterGPs, double* _trCurveSlaveGPs, double* _trCurveGPWeights,
							double* _trCurveMasterGPTangents, double* _trCurveSlaveGPTangents,
							double* _trCurveGPJacobianProducts){
  
    trCurveNumGP = _trCurveNumGP;

    trCurveMasterGPs = new double[2*_trCurveNumGP];
    trCurveSlaveGPs = new double[2*_trCurveNumGP];
    trCurveGPWeights = new double[_trCurveNumGP];
    trCurveMasterGPTangents = new double[3*_trCurveNumGP];
    trCurveSlaveGPTangents = new double[3*_trCurveNumGP];
    trCurveGPJacobianProducts = new double[_trCurveNumGP];

    for (int i = 0; i < trCurveNumGP; i++){
        for (int j = 0; j < 2; j++){
            trCurveMasterGPs[2*i+j]=_trCurveMasterGPs[2*i+j];
            trCurveSlaveGPs[2*i+j]=_trCurveSlaveGPs[2*i+j];
        }
        trCurveGPWeights[i] = _trCurveGPWeights[i];
        for (int j = 0; j < 3; j++){
            trCurveMasterGPTangents[3*i+j]=_trCurveMasterGPTangents[3*i+j];
            trCurveSlaveGPTangents[3*i+j]=_trCurveSlaveGPTangents[3*i+j];
        }
        trCurveGPJacobianProducts[i] = _trCurveGPJacobianProducts[i];
    }

    isGPDataInitialized = true;
}

void WeakIGAPatchContinuityCondition::createGPData(IGAPatchSurface* _masterPatch, IGAPatchSurface* _slavePatch) {

        int noCoordParam = 2;
        int noCoord = 3;

        // Parameter and coordinate sets to store the intersections and transfer them
        std::vector<double> masterUTilde;
        std::vector<double> masterUV;
        std::vector<double> masterXYZ;
        std::vector<double> slaveUTilde;
        std::vector<double> slaveUV;
        std::vector<double> slaveXYZ;
        std::vector<double> slaveUTildeFromMaster;
        std::vector<double> slaveUVFromMaster;
        std::vector<double> slaveXYZFromMaster;

        // Compute the knot intersections of the master trimming curve and their respective parameters/coordinates on the master patch
        _masterPatch->computeKnotIntersectionsWithTrimmingCurve(masterUTilde,masterUV,masterXYZ,
                                                               masterPatchBLIndex,masterPatchBLTrCurveIndex);

        // Compute the knot intersections of slave the trimming curve and their respective parameters/coordinates on the slave patch
        _slavePatch->computeKnotIntersectionsWithTrimmingCurve(slaveUTilde,slaveUV,slaveXYZ,
                                                              slavePatchBLIndex,slavePatchBLTrCurveIndex);

        // Compute the coordinates of the knot intersections of the master trimming curve on the slave trimming curve and slave patch
        _slavePatch->computePointProjectionOnTrimmingCurve(slaveUTildeFromMaster,slaveUVFromMaster,slaveXYZFromMaster,
                                                          masterXYZ,slavePatchBLIndex,slavePatchBLTrCurveIndex);

//        std::cout << "initialSlaveUTilde" << std::endl;
//        for (int i=0; i< slaveUTilde.size();i++)
//            std::cout << slaveUTilde[i] << std::endl;

//        std::cout << "initialslaveUV" << std::endl;
//        for (int i=0; i< slaveUV.size()/2;i++)
//            std::cout << slaveUV[i*2] << " " << slaveUV[i*2+1] << std::endl;

//        std::cout << "initialXYZCoords" << std::endl;
//        for (int i=0; i< slaveXYZ.size()/3;i++)
//            std::cout << slaveXYZ[i*3] << " " << slaveXYZ[i*3+1] << " " << slaveXYZ[i*3+2] << std::endl;

//        std::cout << "slaveUTildeFromMaster" << std::endl;
//        for (int i=0; i< slaveUTildeFromMaster.size();i++)
//            std::cout << slaveUTildeFromMaster[i] << std::endl;

//        std::cout << "slaveUVFromMaster" << std::endl;
//        for (int i=0; i< slaveUVFromMaster.size()/2;i++)
//            std::cout << slaveUVFromMaster[i*2] << " " << slaveUVFromMaster[i*2+1] << std::endl;

//        std::cout << "slaveXYZFromMaster" << std::endl;
//        for (int i=0; i< slaveXYZFromMaster.size()/3;i++)
//            std::cout << slaveXYZFromMaster[i*3] << " " << slaveXYZFromMaster[i*3+1] << " " << slaveXYZFromMaster[i*3+2] << std::endl;

        // Sort the uTilde, uvParams and xyzCoords of both intersections in ascending order with respect to UTilde

        // Initiate iterators
        std::vector<double>::iterator iSlaveUTilde = slaveUTilde.begin();
        std::vector<double>::iterator iSlaveUTildeFromMaster = slaveUTildeFromMaster.begin();

        // Loop over each uTilde of master curve knot intersections on the slave curve
        while (iSlaveUTildeFromMaster != slaveUTildeFromMaster.end()){

            // Loop over each uTilde of slave curve knot intersections
            while (iSlaveUTilde != slaveUTilde.end()) {

                // If the slaveUTildeFromMaster is greater, advance the iterator on the slaveUTilde
                if (*iSlaveUTildeFromMaster > *iSlaveUTilde)  iSlaveUTilde++;

                // If the uTildes are equal break loop
                if (*iSlaveUTildeFromMaster == *iSlaveUTilde)   break;

                // If the slaveUTildeFromMaster is smaller than slaveUTilde, insert into place
                else if (*iSlaveUTildeFromMaster < *iSlaveUTilde && fabs(*iSlaveUTildeFromMaster - *iSlaveUTilde) > 1e-6){

                    // Compute index to insert uTilde as an integer
                    int indexUTilde = std::distance(slaveUTilde.begin(),iSlaveUTilde);

                    // Compute index to insert uvParam as an integer
                    int indexUV = indexUTilde*noCoordParam;

                    // Compute index to insert xyzCoords as an integer
                    int indexXYZCoords = indexUTilde*noCoord;

                    // Compute index to retrieve from slaveUVFromMaster as an integer
                    int indexUVFromMaster = std::distance(slaveUTildeFromMaster.begin(),iSlaveUTildeFromMaster)*noCoordParam;

                    // Compute index to retrieve from slaveXYZFromMaster as an integer
                    int indexXYZFromMaster = std::distance(slaveUTildeFromMaster.begin(),iSlaveUTildeFromMaster)*noCoord;

                    // Insert UTilde
                    slaveUTilde.insert(iSlaveUTilde,*iSlaveUTildeFromMaster);

                    // Insert UV
                    slaveUV.insert(slaveUV.begin()+indexUV,*(slaveUVFromMaster.begin()+indexUVFromMaster+1));
                    slaveUV.insert(slaveUV.begin()+indexUV,*(slaveUVFromMaster.begin()+indexUVFromMaster));

                    // Insert XYZ
                    slaveXYZ.insert(slaveXYZ.begin()+indexXYZCoords,*(slaveXYZFromMaster.begin()+indexXYZFromMaster+2));
                    slaveXYZ.insert(slaveXYZ.begin()+indexXYZCoords,*(slaveXYZFromMaster.begin()+indexXYZFromMaster+1));
                    slaveXYZ.insert(slaveXYZ.begin()+indexXYZCoords,*(slaveXYZFromMaster.begin()+indexXYZFromMaster));

                    // Reset the iterator after each insert
                    iSlaveUTilde = slaveUTilde.begin();
                }
                iSlaveUTilde++;
            }
            iSlaveUTildeFromMaster++;
        }

//        std::cout << "sortedSlaveUTilde" << std::endl;
//        for (int i=0; i< slaveUTilde.size();i++)
//            std::cout << slaveUTilde[i] << std::endl;

//        std::cout << "sortedslaveUV" << std::endl;
//        for (int i=0; i< slaveUV.size()/2;i++)
//            std::cout << slaveUV[i*2] << " " << slaveUV[i*2+1] << std::endl;

//        std::cout << "sortedslaveXYZ" << std::endl;
//        for (int i=0; i< slaveXYZ.size()/3;i++)
//            std::cout << slaveXYZ[i*3] << " " << slaveXYZ[i*3+1] << " " << slaveXYZ[i*3+2] << std::endl;

}

WeakIGAPatchContinuityCondition::~WeakIGAPatchContinuityCondition() {
    if (isGPDataInitialized) {
        delete trCurveMasterGPs;
        delete trCurveSlaveGPs;
        delete trCurveGPWeights;
        delete trCurveMasterGPTangents;
        delete trCurveSlaveGPTangents;
        delete trCurveGPJacobianProducts;
    }
}


} /* namespace EMPIRE */
