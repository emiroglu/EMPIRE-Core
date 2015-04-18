/*
 * DataFieldIntegrationAdapter.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: fabien
 */

#include "DataFieldIntegrationAdapter.h"
#include "AbstractDataFieldIntegration.h"
#include "DataFieldIntegration.h"
#include "DataFieldIntegrationNURBS.h"
#include "FEMesh.h"
#include "IGAMesh.h"
#include "EMPEROR_Enum.h"
#include "Message.h"

namespace EMPIRE {

DataFieldIntegrationAdapter::DataFieldIntegrationAdapter(AbstractMesh *_mesh) :
	mesh(_mesh) {
	if(mesh->type == EMPIRE_Mesh_FEMesh) {
	    FEMesh *feMesh = dynamic_cast<FEMesh*>(_mesh);
	    initDataFieldIntegrationMesh(feMesh);
   } else if(mesh->type == EMPIRE_Mesh_IGAMesh) {
	    IGAMesh *igaMesh = dynamic_cast<IGAMesh*>(_mesh);
	    initDataFieldIntegrationNURBS(igaMesh);
   } else
		ERROR_BLOCK_OUT("DataFieldIntegrationAdapter","DataFieldIntegrationAdapter","DataFieldIntegration not implemented for mesh other than standard FE mesh or NURBS surfaces");
}

DataFieldIntegrationAdapter::~DataFieldIntegrationAdapter() {
	delete dataFieldIntegrationImpl;
}

void DataFieldIntegrationAdapter::initDataFieldIntegrationMesh(FEMesh *_mesh) {
	dataFieldIntegrationImpl=new DataFieldIntegration(_mesh);
}

void DataFieldIntegrationAdapter::initDataFieldIntegrationNURBS(IGAMesh *_mesh) {
	dataFieldIntegrationImpl=new DataFieldIntegrationNURBS(_mesh);
}

void DataFieldIntegrationAdapter::integrate(const double *tractions, double *forces) {
	dataFieldIntegrationImpl->integrate(tractions,forces);
}

void DataFieldIntegrationAdapter::deIntegrate(const double *forces, double *tractions) {
	dataFieldIntegrationImpl->deIntegrate(forces, tractions);
}

} /* namespace EMPIRE */
