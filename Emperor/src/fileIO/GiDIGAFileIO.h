/*  Copyright &copy; 2015, TU Muenchen, Chair of Structural Analysis,
 *  Philipp Bucher, Munich
 */
/***********************************************************************************************//**
 * \file GiDIGAFileIO.h
 * This file holds functions for GiD file formats for NURBS
 * \date 3/15/2013
 **************************************************************************************************/

#ifndef GIDNURBSFILEIO_H_
#define GIDNURBSFILEIO_H_

#include <string>
#include <fstream>
#include "IGAMesh.h"
#include "IGAPatchSurface.h"

namespace EMPIRE {
namespace GiDIGAFileIO {

const std::string SPACE = " ";

// string for header of geometry file
const std::string headerDotGeo = 	"RAMSAN-ASCII-gid-v7.6\n" 	// Header
									"UNKNOWN 0\n"				// Problem Type
									"0\n"						// Must Repair
									"1 Layer0 0 1 153 153 153\n"
									"0\n"
									"0\n\n";

// string for header of result file
const std::string headerDotPostRes = "GiD Post Results File 1.2\n\n"; // Header

/***********************************************************************************************
 * \brief Write a .geo file
 * \param[in] igaMesh Information on the IGA mesh (polynomial orders, knot vectors, Control Points etc.)
 * \author Philipp Bucher
 ***********/
void writeIGAMesh(const std::string _fileName, const EMPIRE::IGAMesh* igaMesh);
/***********************************************************************************************
 * \brief Initialize a .post.res file
 * \param[in] fileName name of the result file
 * \param[in] numNodesPerElem number of nodes per element, used to output Gauss points information
 * \author Philipp Bucher
 ***********/
void initDotPostRes(std::string _fileName);
/***********************************************************************************************
 * \brief Append data of a certain time step to a .post.res file
 * \param[in] fileName name of the result file
 * \param[in] _dataFieldName name data field
 * \param[in] _analysisName name of analysis
 * \param[in] _step step number
 * \param[in] _resultType type "Scalar" or "Vector"
 * \param[in] _dataField data field
 * \param[in] _mesh the nurbs geometry
 * \author Philipp Bucher
 ***********/
void appendCPDataToDotRes(std::string _fileName, std::string _dataFieldName, std::string _analysisName, int _step, std::string _resultType, DataField* _dataField, const IGAMesh* const _mesh);
/***********************************************************************************************
 * \brief Create a linear approximation for this  trimming curve using (#CPs times pDegree of the curve) nodes
 *  This allows a scalable linearization, quite refined.
 * \param[in] polylines vector to store parametric coordinates of linearized NURBS curve
 * \param[in] curve IGA-Curve to be linearized
 * \param[in] curveDirection orientation of curve
 * \author Philipp Bucher
 ***********/
void linearizeTrimmingCurve(std::vector<double>& polylines, const EMPIRE::IGAPatchCurve* curve, bool curveDirection);
/***********************************************************************************************
 * \brief compute number of sample points for linearization of a IGA curve
 * \param[in] p degree of curve to linearize
 * \param[in] noCPs number of control points of curve to linearize
 * \author Philipp Bucher
 ***********/
int getNoLinearizationPoints(const int p, const int noCPs);
/***********************************************************************************************
 * \brief rescale the span of the knot vector
 * This function is used to rescale the knotvector from the existing span (knotVector[0] to knotVector[noKnots-1])
 * to a new span (firstKnot to lastKnot). Currently the firstKnot is always 0 which results in a span from 0 to
 * lastKnot!
 * \param[in] knotVector pointer to knotvector
 * \param[in] noKnots number of knots
 * \param[in] firstKnot beginning of span of the knotvector
 * \param[in] lastKnot end of span of the knotvector
 * \author Philipp Bucher
 ***********/
void rescaleKnotVector(double* knotVector, int noKnots, double firstKnot=0, double lastKnot=1);
/***********************************************************************************************
 * \brief write coordinates (3) to file
 * \param[in] myfile file to write coordinates
 * \param[in] _coords coordinates to write (3)
 * \author Philipp Bucher
 ***********/
void WriteCoordinates(std::ostream& myfile, double* _coords);
/***********************************************************************************************
 * \brief write point to file
 * \param[in] myfile file to write point
 * \param[in] _coords coordinates to write (3)
 * \param[in] pointCounter counter of points
 * \author Philipp Bucher
 ***********/
void WritePoint(std::ostream& myfile, double* _coords, int &pointCounter);
/***********************************************************************************************
 * \brief write line segment to file
 * \param[in] myfile file to write line segment
 * \param[in] curveCounter counter of curves
 * \param[in] startPoint start-point ID of segment
 * \param[in] endPoint end-point ID of segment
 * \author Philipp Bucher
 ***********/
void WriteLineSegment(std::ostream& myfile, int &curveCounter, int startPoint, int endPoint);
/***********************************************************************************************
 * \brief write surface header to file
 * \param[in] myfile file to write header
 * \param[in] surfaceCounter counter of surfaces
 * \author Philipp Bucher
 ***********/
void WriteSurfaceHeader(std::ostream& myfile, int &surfaceCounter);
/***********************************************************************************************/

} /* namespace GiDFileIO */
} /* namespace EMPIRE */
#endif /* GIDNURBSFILEIO_H_ */
