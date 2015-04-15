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

#include "DataFieldIntegrationNURBS.h"
#include "IGAMesh.h"
#include "IGAPatchSurface.h"
#include "MathLibrary.h"
#include "ClipperAdapter.h"
#include "TriangulatorAdaptor.h"
#include "Message.h"

namespace EMPIRE {

using std::make_pair;

const int DataFieldIntegrationNURBS::numGPsMassMatrixTri =  16;
const int DataFieldIntegrationNURBS::numGPsMassMatrixQuad = 25;

DataFieldIntegrationNURBS::DataFieldIntegrationNURBS(IGAMesh* _mesh) {
	numNodes=_mesh->getNumNodes();
    massMatrix = new EMPIRE::MathLibrary::SparseMatrix<double>(numNodes,false);

    /// 1. Loop over every patch to compute the matrix
	for(unsigned int i=0; i< _mesh->getNumPatches(); i++) {
		IGAPatchSurface* patch = _mesh->getSurfacePatch(i);
		Polygon2D polygonUV(4);
		{
			const double u0 = patch->getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
			const double v0 = patch->getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
			const double u1 = patch->getIGABasis()->getUBSplineBasis1D()->getLastKnot();
			const double v1 = patch->getIGABasis()->getVBSplineBasis1D()->getLastKnot();
			polygonUV[0]=make_pair(u0,v0);
			polygonUV[1]=make_pair(u1,v0);
			polygonUV[2]=make_pair(u1,v1);
			polygonUV[3]=make_pair(u0,v1);
		}
		/**************************/
		/// WARNING change order, clip by knot span first and then by trimming ///
		/**************************/

		/// 1.1 Clip by knot span
		Polygon2D listSpan;
		ListPolygon2D listKnotPolygonUV;
		clipByKnotSpan(patch,polygonUV,listKnotPolygonUV,listSpan);
		/// 1.3 For each subelement output of the trimmed polygon, clip by knot span
		for(int index=0;index<listSpan.size();index++) {
			ClipperAdapter::cleanPolygon(listKnotPolygonUV[index]);
			if(listKnotPolygonUV[index].size()<3)
				continue;
			/// 1.3.1 Init list of trimmed polyggons in case patch is not trimmed
			ListPolygon2D listTrimmedPolygonUV(1, listKnotPolygonUV[index]);
			/// 1.3.2 Apply trimming
			if(patch->isTrimmed())
				clipByTrimming(patch,listKnotPolygonUV[index],listTrimmedPolygonUV);
			/// 1.3.2 For each subelement clipped by knot span, compute canonical element and integrate
			for(int trimmedPolygonIndex=0;trimmedPolygonIndex<listTrimmedPolygonUV.size();trimmedPolygonIndex++) {
				ClipperAdapter::cleanPolygon(listTrimmedPolygonUV[trimmedPolygonIndex]);
				if(listTrimmedPolygonUV[trimmedPolygonIndex].size()<3)
					continue;
				ListPolygon2D triangulatedPolygons = triangulatePolygon(listTrimmedPolygonUV[trimmedPolygonIndex]);
				for(ListPolygon2D::iterator triangulatedPolygon=triangulatedPolygons.begin();
						triangulatedPolygon != triangulatedPolygons.end(); triangulatedPolygon++) {
					/// WARNING hard coded tolerance. Cleaning of triangle. Avoid heavily distorted triangle to go further.
					ClipperAdapter::cleanPolygon(*triangulatedPolygon,1e-6);
					if(triangulatedPolygon->size()<3)
						continue;
					assemble(patch,*triangulatedPolygon,listSpan[index].first,listSpan[index].second);
				}
			}
		}
	}
}

DataFieldIntegrationNURBS::~DataFieldIntegrationNURBS() {
}

void DataFieldIntegrationNURBS::deIntegrate(const double *forces, double *tractions) {
	int numNodesReduced=tableCnn.size();
	if(numNodesReduced) {
		double* tmpVecReduced = new double[numNodesReduced];
		for(int i=0;i<numNodesReduced;i++)
			tmpVecReduced[i]=forces[tableCnn[i]];
		double* tmpVecSolReduced = new double[numNodesReduced];
		// 2. solve C_NN * x_master = x_tmp
		massMatrix->solve(tmpVecSolReduced, tmpVecReduced);
		for(int i=0;i<numNodesReduced;i++)
			tractions[tableCnn[i]]=tmpVecSolReduced[i];
		delete[] tmpVecReduced;
		delete[] tmpVecSolReduced;
    } else {
        massMatrix->solve(tractions,const_cast<double*>(forces));
    }
}

void DataFieldIntegrationNURBS::clipByTrimming(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygonUV) {
	ClipperAdapter c;
	for(int loop=0;loop<_thePatch->getTrimming().getNumOfLoops();loop++) {
		const std::vector<double> clippingWindow=_thePatch->getTrimming().getLoop(loop).getPolylines();
		c.addPathClipper(clippingWindow);
	}
	// Setup filling rule to have for sure clockwise loop as hole and counterclockwise as boundaries
	c.setFilling(ClipperAdapter::POSITIVE, 0);
	c.addPathSubject(_polygonUV);
	c.clip();
	c.getSolution(_listPolygonUV);
}

bool DataFieldIntegrationNURBS::computeKnotSpanOfProjElement(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, int* _span) {
	int minSpanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(
            _polygonUV[0].first);
    int minSpanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(
            _polygonUV[0].second);
    int maxSpanU = minSpanU;
    int maxSpanV = minSpanV;

    for (int nodeCount = 1; nodeCount < _polygonUV.size(); nodeCount++) {
        int spanU = _thePatch->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_polygonUV[nodeCount].first);
        int spanV = _thePatch->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_polygonUV[nodeCount].second);

        if (spanU < minSpanU)
            minSpanU = spanU;
        if (spanU > maxSpanU)
            maxSpanU = spanU;
        if (spanV < minSpanV)
            minSpanV = spanV;
        if (spanV > maxSpanV)
            maxSpanV = spanV;
    }

    bool OnSameKnotSpan=(minSpanU == maxSpanU && minSpanV == maxSpanV);
    if(_span!=NULL) {
		_span[0]=minSpanU;
		_span[1]=maxSpanU;
		_span[2]=minSpanV;
		_span[3]=maxSpanV;
    }
    return OnSameKnotSpan;
}

void DataFieldIntegrationNURBS::clipByKnotSpan(const IGAPatchSurface* _thePatch, const Polygon2D& _polygonUV, ListPolygon2D& _listPolygon, Polygon2D& _listSpan) {
	/// 1.find the knot span which the current element located in.
	//      from minSpanu to maxSpanu in U-direction, and from minSpanV to max SpanV in V-direction
    const double *knotVectorU = _thePatch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();
    const double *knotVectorV = _thePatch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

	int span[4];
	int isOnSameKnotSpan = computeKnotSpanOfProjElement(_thePatch, _polygonUV,span);
	int minSpanU = span[0];
	int maxSpanU = span[1];
	int minSpanV = span[2];
	int maxSpanV = span[3];
	/// If on same knot span then returned the same polygon as input
	if (isOnSameKnotSpan) {
		_listPolygon.push_back(_polygonUV);
		_listSpan.push_back(make_pair(minSpanU,minSpanV));
	/// Else clip the polygon for every knot span window it is crossing
	} else {
		for (int spanU = minSpanU; spanU <= maxSpanU; spanU++) {
			for (int spanV = minSpanV; spanV <= maxSpanV; spanV++) {
				ClipperAdapter c(1e-9);
				if (knotVectorU[spanU] != knotVectorU[spanU + 1]
						&& knotVectorV[spanV] != knotVectorV[spanV + 1]) {
					Polygon2D knotSpanWindow(4);
					knotSpanWindow[0]=make_pair(knotVectorU[spanU],knotVectorV[spanV]);
					knotSpanWindow[1]=make_pair(knotVectorU[spanU+1],knotVectorV[spanV]);
					knotSpanWindow[2]=make_pair(knotVectorU[spanU+1],knotVectorV[spanV+1]);
					knotSpanWindow[3]=make_pair(knotVectorU[spanU],knotVectorV[spanV+1]);
					// WARNING : here we assume to get only a single output polygon from the clipping !
					Polygon2D solution = c.clip(_polygonUV,knotSpanWindow);
					/// Store polygon and its knot span for integration
					_listPolygon.push_back(solution);
					_listSpan.push_back(make_pair(spanU,spanV));
				}
			}
		}
	}
}

DataFieldIntegrationNURBS::ListPolygon2D DataFieldIntegrationNURBS::triangulatePolygon(const Polygon2D& _polygonUV) {
	// If already easily integrable by quadrature rule, do nothing
	if(_polygonUV.size()<4)
		return ListPolygon2D(1,_polygonUV);
	// Otherwise triangulate polygon
	TriangulatorAdaptor triangulator;
	// Fill adapter
	for(Polygon2D::const_iterator it=_polygonUV.begin();it!=_polygonUV.end();it++)
		triangulator.addPoint(it->first,it->second,0);
	// Triangulate
	int numTriangles = _polygonUV.size() - 2;
	int triangleIndexes[3 * numTriangles];
	bool triangulated=triangulator.triangulate(triangleIndexes);
	if(!triangulated)
		return ListPolygon2D();
	// Fill output structure
	ListPolygon2D out(numTriangles, Polygon2D(3));
	for(int i = 0; i < numTriangles; i++)
		for(int j=0;j<3;j++)
			out[i][j]=_polygonUV[triangleIndexes[3*i + j]];
	return out;
}

void DataFieldIntegrationNURBS::assemble(IGAPatchSurface* _thePatch, Polygon2D _polygonUV,
        int _spanU, int _spanV) {
    assert(!_polygonUV.empty());
    int numNodesUV=_polygonUV.size();
    assert(numNodesUV>2);
    assert(numNodesUV<5);

    // Definitions
    int numNodesElMaster = 0;

    int pDegree = _thePatch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = _thePatch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
    int nShapeFuncsIGA = (pDegree + 1) * (qDegree + 1);

    numNodesElMaster = nShapeFuncsIGA;

    // Initialize element matrix
    double elementCouplingMatrixNN[numNodesElMaster * (numNodesElMaster + 1) / 2];

    for (int arrayIndex = 0; arrayIndex < numNodesElMaster * (numNodesElMaster + 1) / 2;
            arrayIndex++)
        elementCouplingMatrixNN[arrayIndex] = 0;

/// 1. Copy input polygon into contiguous C format
	double nodesUV[8];
	for (int i = 0; i < numNodesUV; i++) {
		nodesUV[i*2]=_polygonUV[i].first;
		nodesUV[i*2+1]=_polygonUV[i].second;
	}

/// 2. Loop through each quadrature
/// 2.1 Choose gauss triangle or gauss quadriliteral
	MathLibrary::IGAGaussQuadrature *theGaussQuadrature;
	int nNodesQuadrature=numNodesUV;
	if (nNodesQuadrature == 3)
		theGaussQuadrature = new MathLibrary::IGAGaussQuadratureOnTriangle(numGPsMassMatrixTri);
	else
		theGaussQuadrature = new MathLibrary::IGAGaussQuadratureOnQuad(numGPsMassMatrixQuad);

	double *quadratureUV = nodesUV;

/// 2.2 Loop throught each Gauss point
	for (int GPCount = 0; GPCount < theGaussQuadrature->numGaussPoints; GPCount++) {

		/// 2.2.1 compute shape functions from Gauss points(in the quadrature).
		const double *GP = theGaussQuadrature->getGaussPoint(GPCount);

		double shapeFuncs[nNodesQuadrature];
		MathLibrary::computeLowOrderShapeFunc(nNodesQuadrature, GP, shapeFuncs);

		/// 2.2.2 evaluate the coordinates in IGA patch from shape functions
		double GPIGA[2];
		MathLibrary::computeLinearCombination(nNodesQuadrature, 2, quadratureUV, shapeFuncs,
				GPIGA);

		int derivDegree = 1;

		/// 2.2.5 Compute the local basis functions(shape functions of IGA) and their derivatives(for Jacobian)
		double localBasisFunctionsAndDerivatives[(derivDegree + 1) * (derivDegree + 2)
				* nShapeFuncsIGA / 2];

		_thePatch->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
				localBasisFunctionsAndDerivatives, derivDegree, GPIGA[0], _spanU, GPIGA[1],
				_spanV);

		/// 2.2.6 Compute the Jacobian from parameter space on IGA patch to physical
		double baseVectors[6];
		_thePatch->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, _spanU,
				_spanV);

		double JacobianUVToPhysical = MathLibrary::computeAreaTriangle(baseVectors[0],
				baseVectors[1], baseVectors[2], baseVectors[3], baseVectors[4], baseVectors[5])
				* 2;

		/// 2.2.7 Compute the Jacobian from the canonical space to the parameter space of IGA patch
		double JacobianCanonicalToUV;
		if (nNodesQuadrature == 3) {
			JacobianCanonicalToUV = MathLibrary::computeAreaTriangle(
					quadratureUV[2] - quadratureUV[0], quadratureUV[3] - quadratureUV[1], 0,
					quadratureUV[4] - quadratureUV[0], quadratureUV[5] - quadratureUV[1], 0);
		} else {
			double dudx = .25
					* (-(1 - GP[2]) * quadratureUV[0] + (1 - GP[2]) * quadratureUV[2]
							+ (1 + GP[2]) * quadratureUV[4] - (1 + GP[2]) * quadratureUV[6]);
			double dudy = .25
					* (-(1 - GP[1]) * quadratureUV[0] - (1 + GP[1]) * quadratureUV[2]
							+ (1 + GP[1]) * quadratureUV[4] + (1 - GP[1]) * quadratureUV[6]);
			double dvdx = .25
					* (-(1 - GP[2]) * quadratureUV[1] + (1 - GP[2]) * quadratureUV[3]
							+ (1 + GP[2]) * quadratureUV[5] - (1 + GP[2]) * quadratureUV[7]);
			double dvdy = .25
					* (-(1 - GP[1]) * quadratureUV[1] - (1 + GP[1]) * quadratureUV[3]
							+ (1 + GP[1]) * quadratureUV[5] + (1 - GP[1]) * quadratureUV[7]);
			JacobianCanonicalToUV = fabs(dudx * dvdy - dudy * dvdx);
		}
		double Jacobian = JacobianUVToPhysical * JacobianCanonicalToUV;
		/// 2.2.8 integrate the shape function product for C_NN(Linear shape function multiply linear shape function)
		int count = 0;
		for (int i = 0; i < numNodesElMaster; i++) {
			for (int j = i; j < numNodesElMaster; j++) {
					double IGABasisFctsI =
							localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
									1, 0, 0, i)];
					double IGABasisFctsJ =
							localBasisFunctionsAndDerivatives[_thePatch->getIGABasis()->indexDerivativeBasisFunction(
									1, 0, 0, j)];
					elementCouplingMatrixNN[count++] += IGABasisFctsI * IGABasisFctsJ * Jacobian
							* theGaussQuadrature->weights[GPCount];
			}
		}
	}
/// 3.Assemble the element coupling matrix to the global coupling matrix.
    int dofIGA[nShapeFuncsIGA];
    _thePatch->getIGABasis()->getBasisFunctionsIndex(_spanU, _spanV, dofIGA);

    for (int i = 0; i < nShapeFuncsIGA; i++)
        dofIGA[i] = _thePatch->getControlPointNet()[dofIGA[i]]->getDofIndex();

    int count = 0;
    int dof1, dof2;
    for (int i = 0; i < numNodesElMaster; i++) {
        for (int j = i; j < numNodesElMaster; j++) {
			dof1 = dofIGA[i];
			dof2 = dofIGA[j];
			//Because matrix not instantiated as symmetric !!!
            //if (dof1 < dof2)
                (*massMatrix)(dof1, dof2) += elementCouplingMatrixNN[count];
            //else
            if(dof1 != dof2) // if different from main diagonal add symmetric value
                (*massMatrix)(dof2, dof1) += elementCouplingMatrixNN[count];
            count++;
        }
    }
}

void DataFieldIntegrationNURBS::reduceCnn() {
	tableCnn.reserve(numNodes);
	/// Build the index table
	for(int i=0;i<numNodes;i++) {
		if(!massMatrix->isRowEmpty(i))
			tableCnn.push_back(i);
	}
	int numNodesReduced=tableCnn.size();
	/// If nothing missing then do not recreate a matrix and clear table
	if(numNodesReduced == numNodes) {
		tableCnn.clear();
		return;
	}
	/// Create and fill the reduced matrix
	MathLibrary::SparseMatrix<double>* tmp_Cnn = new MathLibrary::SparseMatrix<double>(numNodesReduced, true);
	for(int i=0;i<numNodesReduced;i++)
		for(int j=i;j<numNodesReduced;j++)
			(*tmp_Cnn)(i,j)=(*massMatrix)(tableCnn[i],tableCnn[j]);
	/// Replace the matrix
	delete massMatrix;
	massMatrix=tmp_Cnn;
}

} /* namespace EMPIRE */
