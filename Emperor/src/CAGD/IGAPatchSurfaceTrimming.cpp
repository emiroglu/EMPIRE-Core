/*  Copyright &copy; 2014, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Andreas Apostolatos, Chenshen Wu, Munich
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
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits>

// Inclusion of user defined libraries
#include "IGAPatchSurfaceTrimming.h"
#include "ClipperAdapter.h"
#include "MathLibrary.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

IGAPatchSurfaceTrimming::IGAPatchSurfaceTrimming():outter() {
}

IGAPatchSurfaceTrimming::~IGAPatchSurfaceTrimming() {
	for(int i=0;i<loops.size();i++)
		delete loops[i];
}

void IGAPatchSurfaceTrimming::addTrimLoop(int _inner, int _numCurves) {
    if(_inner) {
        loops.push_back(new IGAPatchSurfaceTrimmingLoop(_numCurves));
    	return;
    } else {
		outter.push_back(loops.size());
        loops.push_back(new IGAPatchSurfaceTrimmingLoop(_numCurves));
	    return;
    }
}

void IGAPatchSurfaceTrimming::addTrimCurve(int _direction,int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
                                           int _uNoControlPoints, double* _controlPointNet) {
    // Read input
    bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
    
    if (ucondition) {
        ERROR_OUT() << " in IGAPatchSurfaceTrimming::IGAPatchSurfaceTrimming" << endl;
        ERROR_OUT()
        << "Number of Control Points, number of knots and polynomial degree do not match!"
        << endl;
        exit(-1);
    }
    // Get the loop currently worked on
    IGAPatchSurfaceTrimmingLoop& loop=*loops.back();
    // Check that size is not going over allocated space made during instantiation
    assert(loop.IGACurves.size()<loop.IGACurves.capacity());
    assert(loop.direction.size()<loop.direction.capacity());
    // Add direction of the curve
    loop.direction.push_back(_direction);
    // Create the NURBS or the B-Spline underlying basis
    loop.IGACurves.push_back(new IGAPatchCurve(_IDBasis, _pDegree, _uNoKnots, _uKnotVector,_uNoControlPoints,_controlPointNet));
}

void IGAPatchSurfaceTrimming::linearizeLoops() {
     for(int i=0;i<loops.size();i++) {
    	loops[i]->linearize();
    }
}

IGAPatchSurfaceTrimmingLoop::IGAPatchSurfaceTrimmingLoop(int _numCurves) {
	 // Reserve the place for the vectors
     IGACurves.reserve(_numCurves);
	 direction.reserve(_numCurves);
     assert(IGACurves.size()==0);
     assert(IGACurves.capacity()!=0);
}

IGAPatchSurfaceTrimmingLoop::~IGAPatchSurfaceTrimmingLoop() {
	for(int i=0;i<getNoCurves();i++)
        delete IGACurves[i];
}

void IGAPatchSurfaceTrimmingLoop::linearize(int _type) {

    // Linearize each trimming curve using the combined algorithm and
    // add the curve linearizations to the trimming loop linearization
    std::vector<double> curvePolyline;
    for(int j=0;j<IGACurves.size();j++) {
        IGACurves[j]->linearize(_type, direction[j]);
//        IGACurves[j]->linearizeUsingGreville(_type, direction[j]);
//        IGACurves[j]->linearizeUsingNCPxP(_type, direction[j]);
        curvePolyline = *(IGACurves[j]->getPolyline());
        for (int iVertex = 0; iVertex < curvePolyline.size()/2; iVertex++) {
            polylines.push_back(curvePolyline[iVertex*2]);
            polylines.push_back(curvePolyline[iVertex*2+1]);
        }
    }
    ClipperAdapter::cleanPolygon(polylines);
}

Message &operator<<(Message &message, const IGAPatchSurfaceTrimming &trim) {
    message << "\t" << "---------------------------------Start Trimming" << endl;
    message << "\t" << "Trimming Info"<< endl;
    message << "\t\t" << "Number of loops: "<<trim.getNumOfLoops()<< endl;
    message << "\t\t" << "Outer loop size: "<<trim.getOutterLoopIndex().size()<< endl;
	for(int i=0; i<trim.getOutterLoopIndex().size();i++) {
	    message << "\t" << "Outter Loop["<<trim.getOutterLoopIndex(i)<<"]"<< endl;
		message << trim.getLoop(i);
	}
	for(int i=0; i<trim.getNumOfLoops();i++) {
		if(find(trim.getOutterLoopIndex().begin(),trim.getOutterLoopIndex().end(),i) != trim.getOutterLoopIndex().end())
				continue;
    	message << "\t" << "InnerLoop["<<i<<"]"<< endl;
		message << trim.getLoop(i);
	}
    message << "\t" << "---------------------------------End Trimming" << endl;
	return message;
}

Message &operator<<(Message &message, const IGAPatchSurfaceTrimmingLoop &trim) {
    /// output loop
    for(int i=0;i<trim.getIGACurves().size();++i) {
        message << "\t" << "Curve["<<i<<"]"<< endl;
        message << "\t\tpDegree:  " << trim.getIGACurve(i).getIGABasis()->getPolynomialDegree()<< endl;

        message << "\t\tKnots Vector U: \t";
        for (int k = 0; k < trim.getIGACurve(i).getIGABasis()->getNoKnots(); k++)
            message << trim.getIGACurve(i).getIGABasis()->getKnotVector()[k] << "  ";
        message << endl;
        message << "\t\t" << "number of control points: " << trim.getNoControlPoints(i) << endl;

        message << "\t\tControl Points Net: " << endl;
        for (int k = 0; k < trim.getNoControlPoints(i); k++) {
            message << "\t\t";
            message << trim.getIGACurve(i).getControlPoint(k).getX() << ", "
                    << trim.getIGACurve(i).getControlPoint(k).getY() << ", "
                    << trim.getIGACurve(i).getControlPoint(k).getZ() << endl;
        }
    }
    message << "\tLinear Polygon: " << endl;
    for (int k = 0; k < trim.getPolylines().size()/2; k++) {
        message << "\t\t";
        message << trim.getPolylines()[2*k] << ", "
                << trim.getPolylines()[2*k+1]<< endl;
    }
    return message;
}

}/* namespace EMPIRE */
