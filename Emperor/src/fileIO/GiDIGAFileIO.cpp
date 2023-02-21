/*
 * GiDNurbsFileIO.cpp
 *
 *  Created on: May 13, 2015
 *      Author: Philipp Bucher
 */

#include "GiDIGAFileIO.h"

#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "BSplineBasis2D.h"
#include "BSplineBasis1D.h"
#include "IGAControlPoint.h"
#include "DataField.h"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <assert.h>
#include <vector>
#include "stdlib.h"
#include <math.h>

#include <sys/stat.h>

using namespace std;
namespace EMPIRE {
namespace GiDIGAFileIO {

void writeIGAMesh(const std::string _fileName, const EMPIRE::IGAMesh* igaMesh) {
	// set up names and path for files
	string folderName = _fileName + ".gid";
	string geometryFileName = _fileName + ".geo";
	string geometryFileAndPathName = folderName + "/" + geometryFileName;

	// create .gid-folder for files .geo and .post.res
	/*
	 To make it plattform independent, the flags must be implemented in CMAKELISTS
	 if(-DUSE_LINUX==true){// create folder in linux
	 mkdir(folderName.c_str(),0777);
	 }
	 if(-DUSE_WINDOWS==true){// create folder in windows
	 mkdir(folderName.c_str());
	 }*/

	mkdir(folderName.c_str(), 0777); // create folder in linux => delete once the flags work

	ofstream dotGeoFile(geometryFileAndPathName.c_str()); // open geometry file in created folder
	assert(!dotGeoFile.fail());
	dotGeoFile.precision(12);
	dotGeoFile << std::dec;

	// start writing data to file
	dotGeoFile << headerDotGeo;	// write header

	int pointCounter = 1;
	int curveCounter = 1;
	int surfaceCounter = 1;
	double cartCoord[3], surfCoord[2];
	int noSampPoints;
	int startIndex;
	int endIndex;
	int firstLastPoint;
	std::vector<double> polylines;// vector to store parametric coordinates of linearized NURBS curve

	const int numPatches = igaMesh->getNumPatches();
	
	for (int indexPatch = 0; indexPatch < numPatches; indexPatch++) { // loop over patches in IGA-mesh
		const EMPIRE::IGAPatchSurface* patch = igaMesh->getSurfacePatch(
				indexPatch);
		const int numCPSs = patch->getNoControlPoints();
		int curveCounterBeginning = curveCounter;
		bool patchIsTrimmed = patch->getTrimming().isTrimmed();
		
		// ********** trimmed patch **********
		if (patchIsTrimmed) {
			
			const int numLoops = patch->getTrimming().getNumOfLoops();
			for (int indexLoop = 0; indexLoop < numLoops; indexLoop++) { // loop over trimming loops
				const EMPIRE::IGAPatchSurfaceTrimmingLoop* loop =
						&patch->getTrimming().getLoop(indexLoop);
				const int numCurves = loop->getNoCurves();
				firstLastPoint = pointCounter;

				for (int indexCurve = 0; indexCurve < numCurves; indexCurve++) { // loop over curves in trimming loop
					const EMPIRE::IGAPatchCurve* curve = &loop->getIGACurve(
							indexCurve);

					linearizeTrimmingCurve(polylines, curve,
							loop->getDirection(indexCurve)); // compute linearization of NURBS curve

					noSampPoints = polylines.size() / 2;

					startIndex = 1;
					endIndex = noSampPoints;

					if (indexCurve == 0) // first curve in loop
						startIndex = 0;

					if (indexCurve == numCurves - 1) // last curve in loop
						endIndex = noSampPoints - 1;

					// write points
					for (int i = startIndex; i < endIndex; i++) {
						surfCoord[0] = polylines.at(2 * i);
						surfCoord[1] = polylines.at(2 * i + 1);
						patch->computeCartesianCoordinates(cartCoord,
								surfCoord); // compute coord in cartesian coordinates
						WritePoint(dotGeoFile, cartCoord, pointCounter); // write point
					}

					// write line segments
					for (int i = 0; i < noSampPoints - 2; i++)
						WriteLineSegment(dotGeoFile, curveCounter, curveCounter,
								curveCounter + 1);
					// write last segment
					if (indexCurve != numCurves - 1) // check if curve is last curve in loop
						WriteLineSegment(dotGeoFile, curveCounter, curveCounter,
								curveCounter + 1);
					else
						WriteLineSegment(dotGeoFile, curveCounter, curveCounter,
								firstLastPoint);

					dotGeoFile << endl;
				} // loop curve
			} // loop trimming

			
		} // end trimming

		// ********** untrimmed patch **********
		else {
			
			/* 	   4	 |2|	 3
			 * 		x-----------x
			 * 		|			|
			 * 		|			|
			 * 	|3|	|	PATCH	| |1|
			 * 		|			|
			 * 		|			|
			 * 		x-----------x
			 * 	   1	 |0|	 2 
			 */
			
			int noSampPoints_U = getNoLinearizationPoints(patch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree(), 
					patch->getUNoControlPoints());
			int noSampPoints_V = getNoLinearizationPoints(patch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree(), 
					patch->getVNoControlPoints());
			firstLastPoint = pointCounter;
			
			const double U_FirstKnot = patch->getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
			const double U_LastKnot = patch->getIGABasis()->getUBSplineBasis1D()->getLastKnot();

			const double V_FirstKnot = patch->getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
			const double V_LastKnot = patch->getIGABasis()->getVBSplineBasis1D()->getLastKnot();
			
			bool isOnU;
			double UV_Start;
			double UV_End;
			double UV_Fixed;
			double U;
			double V;
			double dUV_Running;
			
			// Loop over the 4 boundary curves of the patch
			for (int indexCurve = 0; indexCurve < 4; indexCurve++) {
				// Get fixed and running parameters on the boundary curve
				if (indexCurve == 0){
					isOnU = true;
					UV_Start = U_FirstKnot;
					UV_End = U_LastKnot;
					UV_Fixed = V_FirstKnot;
					U = UV_Start;
					V = UV_Fixed;
					noSampPoints = noSampPoints_U;
				}
				else if (indexCurve == 1) {
					isOnU = false;
					UV_Start = V_FirstKnot;
					UV_End = V_LastKnot;
					UV_Fixed = U_LastKnot;
					U = UV_Fixed;
					V = UV_Start;
					noSampPoints = noSampPoints_V;
				}
				else if (indexCurve == 2) {
					isOnU = true;
					UV_Start = U_LastKnot;
					UV_End = U_FirstKnot;
					UV_Fixed = V_LastKnot;
					U = UV_Start;
					V = UV_Fixed;
					noSampPoints = noSampPoints_U;
				}
				else {
					isOnU = false;
					UV_Start = V_LastKnot;
					UV_End = V_FirstKnot;
					UV_Fixed = U_FirstKnot;
					U = UV_Fixed;
					V = UV_Start;
					noSampPoints = noSampPoints_V;
				}
				
				// Initialize parameters and counters
				dUV_Running = (UV_End - UV_Start) / (noSampPoints - 1);
				startIndex = 1;
				endIndex = noSampPoints;
				if (indexCurve == 0) {
					startIndex = 0;}
				if (indexCurve == 3) {
					endIndex = noSampPoints - 1;}
				
				// write points
				for (int i = startIndex; i < endIndex; i++) {
					if (isOnU) {
						U = UV_Start + i * dUV_Running;
					}
					else {
						V = UV_Start + i * dUV_Running;
					}
					surfCoord[0] = U;
					surfCoord[1] = V;
					
					patch->computeCartesianCoordinates(cartCoord, surfCoord); // compute coord in cartesian coordinates
					WritePoint(dotGeoFile, cartCoord, pointCounter); // write point				
				}
				
				// write line segments
				for (int i = 0; i < noSampPoints - 2; i++)
					WriteLineSegment(dotGeoFile, curveCounter, curveCounter, curveCounter + 1);
				// write last segment
				if (indexCurve != 3) // check if curve is last curve in loop
					WriteLineSegment(dotGeoFile, curveCounter, curveCounter, curveCounter + 1);
				else
					WriteLineSegment(dotGeoFile, curveCounter, curveCounter, firstLastPoint);
				
				dotGeoFile << endl;
			}
		
		} // end untrimmed 
		
		// start of Patch declaration ------------------------------------------------------------------------------
		WriteSurfaceHeader(dotGeoFile, surfaceCounter); // write SurfaceHeader;

		// write number of curves per patch and the curve IDs
		int NoPatchCurves = curveCounter - curveCounterBeginning; // number of curves for patch
		dotGeoFile << NoPatchCurves << endl;

		for (int curveIndex = curveCounterBeginning;
				curveIndex < curveCounter; curveIndex++)
			dotGeoFile << curveIndex << SPACE; // write curve ID
		dotGeoFile << endl;

		// write orientations of curves
		for (int i = 0; i < NoPatchCurves; i++)
			dotGeoFile << "0 ";
		dotGeoFile << endl;

		// write and compute approx center and corresponding normal (where label is drawn)
		// compute middle of knotspan in U and V => used as approx center
		double U_FirstKnot =
				patch->getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
		double U_LastKnot =
				patch->getIGABasis()->getUBSplineBasis1D()->getLastKnot();

		double V_FirstKnot =
				patch->getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
		double V_LastKnot =
				patch->getIGABasis()->getVBSplineBasis1D()->getLastKnot();

		double U_MidKnotSpan = (U_FirstKnot + U_LastKnot) / 2;
		double V_MidKnotSpan = (V_FirstKnot + V_LastKnot) / 2;

		double center[3];
		double normal[3];

		patch->computeCartesianCoordinatesAndNormalVector(center, normal,
				U_MidKnotSpan, V_MidKnotSpan);

		WriteCoordinates(dotGeoFile, center); // write center
		WriteCoordinates(dotGeoFile, normal); // write normal

		// write further information (N ofControlPoints u/v and pDegree u/v)
		int noCPU = patch->getUNoControlPoints();
		int noCPV = patch->getVNoControlPoints();

		int pDegreeU =
				patch->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
		int pDegreeV =
				patch->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

		dotGeoFile << patchIsTrimmed << SPACE;		
		dotGeoFile << noCPU << SPACE << noCPV << SPACE << pDegreeU << SPACE	<< pDegreeV << endl;

		// compute and write Control Points of patch
		for (int indexCP = 0; indexCP < numCPSs; indexCP++) {
			// compute Point
			cartCoord[0] = patch->getControlPointNet()[indexCP]->getX();
			cartCoord[1] = patch->getControlPointNet()[indexCP]->getY();
			cartCoord[2] = patch->getControlPointNet()[indexCP]->getZ();

			WriteCoordinates(dotGeoFile, cartCoord); // write coordinates
		}

		// write knot vector
		// knot vector U
		int noKnots =
				patch->getIGABasis()->getUBSplineBasis1D()->getNoKnots();
		double* knotVector =
				patch->getIGABasis()->getUBSplineBasis1D()->getKnotVector();

		rescaleKnotVector(knotVector, noKnots);

		for (int indexKnot = 0; indexKnot < noKnots; indexKnot++)
			dotGeoFile << knotVector[indexKnot] << SPACE;
		dotGeoFile << endl;

		// knot vector V
		noKnots = patch->getIGABasis()->getVBSplineBasis1D()->getNoKnots();
		knotVector =
				patch->getIGABasis()->getVBSplineBasis1D()->getKnotVector();

		rescaleKnotVector(knotVector, noKnots);

		for (int indexKnot = 0; indexKnot < noKnots; indexKnot++)
			dotGeoFile << knotVector[indexKnot] << SPACE;
		dotGeoFile << endl;

		// write weights
		dotGeoFile << "1 "; // 1 stands for rational

		for (int indexCP = 0; indexCP < numCPSs; indexCP++)
			dotGeoFile << patch->getControlPointNet()[indexCP]->getW()
					<< SPACE;
		patch->getControlPointNet()[1]->getW();
		dotGeoFile << endl << endl;
		
		// end of Patch declaration ------------------------------------------------------------------------------
		
	} // loop patches

	dotGeoFile << "0"; // indicates the end of the definition of geometrical entities

	dotGeoFile.close();
}

void initDotPostRes(const std::string _fileName) {
	// set up names and path for files
	string folderName = _fileName + ".gid";
	string resultFileName = _fileName + ".post.res";
	string resultFileAndPathName = folderName + "/" + resultFileName;

	ofstream dotPostResFile(resultFileAndPathName.c_str(), ios_base::out); // open result file in previously created folder
	assert(!dotPostResFile.fail());
	dotPostResFile << headerDotPostRes;
	dotPostResFile.close();
}

void appendCPDataToDotRes(std::string _fileName, std::string _dataFieldName,
		std::string _analysisName, int _step, std::string _resultType,
		DataField* _dataField, const IGAMesh* const _mesh) {
	// build patch to dof id table
	// TODO should be done only once and reused, no recomputed
	int numPatches = _mesh->getNumPatches();
	std::map<int, std::vector<int> > patchToDof;
	for (int indexPatch = 0; indexPatch < numPatches; indexPatch++) {
		const EMPIRE::IGAPatchSurface* patch = _mesh->getSurfacePatch(
				indexPatch);
		for (int indexCP = 0; indexCP < patch->getNoControlPoints();
				indexCP++) {
			patchToDof[indexPatch].push_back((*patch)[indexCP]->getDofIndex());
		}
	}

	// set up names and path for files
	string folderName = _fileName + ".gid";
	string resultFileName = _fileName + ".post.res";
	string resultFileAndPathName = folderName + "/" + resultFileName;

	ofstream dotPostResFile(resultFileAndPathName.c_str(), ios_base::app); // open initialized result file in previously created folder
	dotPostResFile.precision(14);
	dotPostResFile << std::dec;

	// set dimension of problem
	int dimension = 0;
	if (_resultType == "vector")
		dimension = 3;
	else if (_resultType == "scalar")
		dimension = 1;
	else {
		WARNING_BLOCK_OUT("GidIFAFileIO", "appendCPDataToDotRes",
				"Dimension of datafield not available for output");
		return;
	}

	// start writing data to file ------------------------------------------------------------------------------
	dotPostResFile << "Result " << "\"" << _dataFieldName << "\"" << SPACE
			<< _analysisName << SPACE << _step << SPACE << _resultType << SPACE
			<< "OnNurbsSurface" << endl;
	dotPostResFile << "Values" << endl;

	// loop over patches
	for (map<int, vector<int> >::iterator itPatch = patchToDof.begin();
			itPatch != patchToDof.end(); itPatch++) {
		dotPostResFile << (itPatch->first) + 1 << endl;
		for (vector<int>::iterator itCP = itPatch->second.begin();
				itCP != itPatch->second.end(); itCP++) {
			for (int d = 0; d < dimension; d++)
				dotPostResFile << _dataField->data[*itCP * dimension + d]
						<< std::setprecision(20) << SPACE;
			dotPostResFile << endl;
		}
	}
	dotPostResFile << "End Values" << endl << endl;
	// end writing data to file ------------------------------------------------------------------------------

	dotPostResFile.close();
}

void linearizeTrimmingCurve(std::vector<double>& polylines,
		const EMPIRE::IGAPatchCurve* curve, bool curveDirection) {
	polylines.clear();
	int noCP = curve->getNoControlPoints();
	int p = curve->getIGABasis()->getPolynomialDegree();
	int noSampPoints = getNoLinearizationPoints(p, noCP);
	
	double u0 = curve->getIGABasis()->getFirstKnot();
	double u1 = curve->getIGABasis()->getLastKnot();
	double du = (u1 - u0) / (noSampPoints - 1);
	
	// Check direction to put points in the right sequence (counter clockwise for outter loop, clockwise for inner
	if (curveDirection) {
		for (int i = 0; i < noSampPoints; i++) {
			double knot = u0 + i * du;
			double parametricCoordinates[2] = { 0 };
			curve->computeCartesianCoordinates(parametricCoordinates, knot);
			polylines.push_back(parametricCoordinates[0]);
			polylines.push_back(parametricCoordinates[1]);
		}
	} else {
		for (int i = noSampPoints - 1; i >= 0; i--) {
			double knot = u0 + i * du;
			double parametricCoordinates[2] = { 0 };
			curve->computeCartesianCoordinates(parametricCoordinates, knot);
			polylines.push_back(parametricCoordinates[0]);
			polylines.push_back(parametricCoordinates[1]);
		}
	}
}

int getNoLinearizationPoints(const int p, const int noCPs) {
	int factor = fmax(0, 4 - p); // 4 is arbitrary degree upon linearization is fully NCPxP
	factor = 1 + factor * factor * factor;
	return noCPs * p * factor; // return number of linearization points
}

void rescaleKnotVector(double* _knotVector, int _noKnots, double _firstKnot,
		double _lastKnot) {
	if (_firstKnot != 0)
		WARNING_OUT("rescaling to non-zero start of span not yet implemented");

	double shift = -_knotVector[0];

	for (int i = 0; i < _noKnots; i++)
		_knotVector[i] += shift;

	double scaling = _lastKnot / _knotVector[_noKnots - 1];

	for (int i = 0; i < _noKnots; i++)
		_knotVector[i] *= scaling;

	/*test to rescale to arbirtary span => doesn't work because of multiplicity!
	 double spacing = fabs((firstKnot-lastKnot)/(noKnots-1));

	 knotVector[0] = firstKnot;

	 for(int i=1;i<noKnots;i++)
	 knotVector[i] = knotVector[0] + spacing * i;*/
}

void WriteCoordinates(std::ostream& myfile, double* _coords) {
	myfile << _coords[0] << SPACE << _coords[1] << SPACE << _coords[2] << endl;
}

void WritePoint(std::ostream& myfile, double* _coords,
		int &pointCounter) {
	myfile << "1 " << pointCounter << " 0 0 2 0 0 2 0" << endl;	// write point header
	pointCounter++;
	WriteCoordinates(myfile, _coords); // write point coordinates
	myfile << endl;
}

void WriteLineSegment(std::ostream& myfile, int &curveCounter,
		int startPoint, int endPoint) {
	myfile << "2 " << curveCounter << " 0 0 1 0 0 2 0" << endl; // write segment header
	myfile << startPoint << SPACE << endPoint << endl;// write points that define segment
	curveCounter++;
}

void WriteSurfaceHeader(std::ostream& myfile, int &surfaceCounter) {
	myfile << "14 " << surfaceCounter << " 0 0 0 0 0 2 0" << endl; // write surface header
	surfaceCounter++;
}
} /* namespace GiDIGAFileIO */
} /* namespace EMPIRE */
