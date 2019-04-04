/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Munich
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
 * \file DataFieldIntegrationAdapter.h
 * This file holds the class DataFieldIntegrationAdapter
 * \date 1/4/2015
 **************************************************************************************************/

#ifndef DATAFIELDINTEGRATIONADAPTER_H_
#define DATAFIELDINTEGRATIONADAPTER_H_

namespace EMPIRE {

class AbstractMesh;
class FEMesh;
class IGAMesh;
class AbstractDataFieldIntegration;

/********//**
 * \brief Class DataFieldIntegrationAdapter is an adapter for different kind of DataFieldIntegration depending of input mesh type
 * \author Fabien Pean
 ***********/
class DataFieldIntegrationAdapter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name name of the mapper
     * \param[in] _mesh mesh
     * \author Fabien Pean
     ***********/
	DataFieldIntegrationAdapter(AbstractMesh *_mesh);
    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean
     ***********/
	virtual ~DataFieldIntegrationAdapter();
    /***********************************************************************************************
     * \brief Initialize initDataFieldIntegrationMesh
     * \author Fabien Pean
     ***********/
    void initDataFieldIntegrationMesh(FEMesh *_mesh);
    /***********************************************************************************************
     * \brief Initialize initDataFieldIntegrationNURBS
     * \author Fabien Pean
     ***********/
    void initDataFieldIntegrationNURBS(IGAMesh *_mesh);
    /***********************************************************************************************
     * \brief Do integration (massMatrix*tractions=forces).
     * \param[in] tractions tractions (in one direction)
     * \param[out] forces forces (in one direction)
     * \author Tianyang Wang
     ***********/
    void integrate(const double *tractions, double *forces);
    /***********************************************************************************************
     * \brief Do deintegration (massMatrix^(-1)*forces=tractions).
     * \param[in] forces forces (in one direction)
     * \param[out] tractions tractions (in one direction)
     * \author Tianyang Wang
     ***********/
    void deIntegrate(const double *forces, double *tractions);
    /***********************************************************************************************
     * \brief is it mesh or not
     * \param[in] _mesh mesh
     * \return true if it is mesh
     * \author Tianyang Wang
     ***********/
    inline bool isMesh(AbstractMesh *_mesh) {
        return mesh == _mesh;
    };
private:
    /// the adapted mapper
    AbstractDataFieldIntegration *dataFieldIntegrationImpl;
    /// mesh
    const AbstractMesh *mesh;
};

} /* namespace EMPIRE */
#endif /* DATAFIELDINTEGRATIONADAPTER_H_ */
