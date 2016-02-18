/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Ragnar Björnsson, Munich
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
 * \file IGAPatchCouplingCaratData.h
 * This file holds the class IGAPatchCouplingCaratData.h
 * \date 6/8/2013
 **************************************************************************************************/

#ifndef IGAPATCHCOUPLINGCARATDATA_H
#define IGAPATCHCOUPLINGCARATDATA_H

#include <string>
#include <cfloat>
// Inclusion of user defined libraries
//#include "AbstractMesh.h"

namespace EMPIRE {

/********//**
 * \brief class IGAMesh is a specialization of the class AbstractMesh used for IGA Mesh containing number of IGA surface patches
 ***********/

class IGAPatchCouplingCaratData
{
private:
    std::vector<std::vector<double*> > AllGP_master;
    std::vector<std::vector<double*> > AllGP_slave;
    std::vector<std::vector<double*> > AllGP_weight;
    std::vector<std::vector<double*> > AllTangents_master;
    std::vector<std::vector<double*> > AllTangents_slave;
    std::vector<std::vector<double*> > AllMappings;
    std::vector<std::vector<int> > AllIDs_slave;
    std::vector<std::vector<int> > AllNumElemsPerBRep;
    std::vector<std::vector<int> > AllNumGPsPerElem;

    int numPatches;
    int* numBrepsPerPatch;


public:

    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numPatches number of patches
     * \param[in] _numBrepsPerPatch the number of BRePs for each patch
     * \author Ragnar Björnsson
     ***********/
    IGAPatchCouplingCaratData(int _numPatches, int* _numBrepsPerPatch) {
        numPatches=_numPatches;
        numBrepsPerPatch=_numBrepsPerPatch;

        AllGP_master.resize(numPatches);
        AllGP_slave.resize(numPatches);
        AllGP_weight.resize(numPatches);
        AllTangents_master.resize(numPatches);
        AllTangents_slave.resize(numPatches);
        AllMappings.resize(numPatches);
        AllIDs_slave.resize(numPatches);
        AllNumElemsPerBRep.resize(numPatches);
        AllNumGPsPerElem.resize(numPatches);

        for(int i = 0 ; i < numPatches ; i++) {
            AllGP_master[i].resize(numBrepsPerPatch[i]);
            AllGP_slave[i].resize(numBrepsPerPatch[i]);
            AllGP_weight[i].resize(numBrepsPerPatch[i]);
            AllTangents_master[i].resize(numBrepsPerPatch[i]);
            AllTangents_slave[i].resize(numBrepsPerPatch[i]);
            AllMappings[i].resize(numBrepsPerPatch[i]);
            AllIDs_slave[i].resize(numBrepsPerPatch[i]);
            AllNumElemsPerBRep[i].resize(numBrepsPerPatch[i]);
            AllNumGPsPerElem[i].resize(numBrepsPerPatch[i]);
        }
    }

    /***********************************************************************************************
     * \brief Destructor
     * \author Ragnar Björnsson
     ***********/
    ~IGAPatchCouplingCaratData() {
        for(int i = 0 ; i < AllGP_master.size() ; i++) {
            for(int j = 0 ; j < AllGP_master[i].size() ; j++) {
                delete[] AllGP_master[i][j];
                delete[] AllGP_slave[i][j];
                delete[] AllGP_weight[i][j];
                delete[] AllTangents_master[i][j];
                delete[] AllTangents_slave[i][j];
                delete[] AllMappings[i][j];
            }
        }
        delete[] numBrepsPerPatch;
    }

    /***********************************************************************************************
     * \brief add coupling data of BReP
     * \param[in] patchCounter patch counter
     * \param[in] BRepCounter BReP counter
     * \param[in] GP_m gauss points on the master side
     * \param[in] GP_s gauss points on the slave side
     * \param[in] GP_w gauss point weights
     * \param[in] tang_m tangents on the master side
     * \param[in] tang_s tangents on the slave side
     * \param[in] map mapping of gauss points from parent element space to the physical space
     * \param[in] ID_s Ids of the slave patch
     * \param[in] NumElemsOfBRep number of linear elements on the BReP
     * \param[in] NumGPsOfElem number of gauss points each elements (they are of uniform length)
     * \author Ragnar Björnsson
     ***********/
    void addBRePCouplingData(int patchCounter, int BRepCounter, double* GP_m, double* GP_s, double* GP_w,
                             double* tang_m, double* tang_s, double* map,
                             int ID_s, int NumElemsOfBRep, int NumGPsOfElem) {
        AllGP_master[patchCounter][BRepCounter] = GP_m;
        AllGP_slave[patchCounter][BRepCounter] = GP_s;
        AllGP_weight[patchCounter][BRepCounter] = GP_w;
        AllTangents_master[patchCounter][BRepCounter] = tang_m;
        AllTangents_slave[patchCounter][BRepCounter] = tang_s;
        AllMappings[patchCounter][BRepCounter] = map;
        AllIDs_slave[patchCounter][BRepCounter] = ID_s;
        AllNumElemsPerBRep[patchCounter][BRepCounter] = NumElemsOfBRep;
        AllNumGPsPerElem[patchCounter][BRepCounter] = NumGPsOfElem;
    }

    /***********************************************************************************************
     * \brief get number of patches
     * \author Ragnar Björnsson
     ***********/
    int getNumPatches() {
        return numPatches;
    }

    /***********************************************************************************************
     * \brief get number of BRePs for each patch
     * \author Ragnar Björnsson
     ***********/
    int* getNumBrepsPerPatch() {
        return numBrepsPerPatch;
    }

    /***********************************************************************************************
     * \brief get gauss points on the master side
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    double* getGPs_master(int patchNum, int BRepNum) {
        return AllGP_master[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get gauss points on the slave side
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    double* getGPs_slave(int patchNum, int BRepNum) {
        return AllGP_slave[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get gauss points weights
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    double* getGPs_weight(int patchNum, int BRepNum) {
        return AllGP_weight[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get tangents on the master side
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    double* getTangents_master(int patchNum, int BRepNum) {
        return AllTangents_master[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get tangents on the slave side
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    double* getTangents_slave(int patchNum, int BRepNum) {
        return AllTangents_slave[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get mapping of gauss points from parent element space to the physical space
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    double* getMappings(int patchNum, int BRepNum) {
        return AllMappings[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get slave patch Id
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    int getSlaveID(int patchNum, int BRepNum) {
        return AllIDs_slave[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get number of elements on a BReP
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    int getNumElemsOfBRep(int patchNum, int BRepNum) {
        return AllNumElemsPerBRep[patchNum][BRepNum];
    }

    /***********************************************************************************************
     * \brief get number of gauss points on each element
     * \param[in] patchNum number of patch
     * \param[in] BRepNum number of BReP
     * \author Ragnar Björnsson
     ***********/
    int getNumGPsOfElem(int patchNum, int BRepNum) {
        return AllNumGPsPerElem[patchNum][BRepNum];
    }
};

} /* namespace EMPIRE */

#endif // IGAPATCHCOUPLINGCARATDATA_H
