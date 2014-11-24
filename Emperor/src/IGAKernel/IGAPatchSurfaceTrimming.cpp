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

IGAPatchSurfaceTrimming::IGAPatchSurfaceTrimming():outter(-1) {
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
		if(outter>=0) {
			ERROR_OUT() << "Outter loop boundary for trimming has already been defined" << endl;
			exit(-1);
		}
		outter=loops.size();
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
    // Check that size is not going over allocated during instantiation
    assert(loop.IGACurve.size()<loop.IGACurve.capacity());
    assert(loop.direction.size()<loop.direction.capacity());
    // Add direction of the curve
    loop.direction.push_back(_direction);
    // Create the NURBS or the B-Spline underlying basis
    loop.IGACurve.push_back(new IGAPatchCurve(_IDBasis, _pDegree, _uNoKnots, _uKnotVector,_uNoControlPoints,_controlPointNet));
}

void IGAPatchSurfaceTrimming::linearizeLoops() {
     for(int i=0;i<loops.size();i++) {
    	loops[i]->linearize();
    }
}

IGAPatchSurfaceTrimmingLoop::IGAPatchSurfaceTrimmingLoop(int _numCurves) {
	 IGACurve.reserve(_numCurves);
	 direction.reserve(_numCurves);
	 assert(IGACurve.size()==0);
	 assert(IGACurve.capacity()!=0);
}

IGAPatchSurfaceTrimmingLoop::~IGAPatchSurfaceTrimmingLoop() {
	for(int i=0;i<getNoCurves();i++)
		delete IGACurve[i];
}

void IGAPatchSurfaceTrimmingLoop::linearize() {
    /*
     * Linearize approximation of the nurbs curves
     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
     *
     * Function layout :
     *
     * 1. For every NURBS curve
     * 1.1. For every control points of the curve
     * 1.1.1. Compute Greville abscissae defined Knot space
     * 1.1.2. Compute basis function at the Greville abscissae
     * 1.1.3. Compute position in parametric space
     * 1.1.4. Store point in data structure of polylines
     */
	for(int j=0;j<IGACurve.size();j++) {
		/// Check direction to put points in the right sequence (counter clockwise for outter loop, clockwise for inner)
		if(direction[j]) {
			for(int cpIndex=0;cpIndex<getIGACurve(j).getNoControlPoints();cpIndex++) {
				double knotGreville=getIGACurve(j).getIGABasis()->computeGrevilleAbscissae(cpIndex);
				double parametricCoordinates[2] = {0};
				getIGACurve(j).computeCartesianCoordinates(parametricCoordinates,knotGreville);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
			}
		} else {
			for(int cpIndex=getIGACurve(j).getNoControlPoints()-1;cpIndex>=0;cpIndex--) {
				double knotGreville=getIGACurve(j).getIGABasis()->computeGrevilleAbscissae(cpIndex);
				double parametricCoordinates[2] = {0};
				getIGACurve(j).computeCartesianCoordinates(parametricCoordinates,knotGreville);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
			}
		}
	}
	ClipperAdapter::cleanPolygon(polylines);
}

Message &operator<<(Message &message, const IGAPatchSurfaceTrimming &trim) {
    message << "\t" << "---------------------------------Start Trimming" << endl;
    message << "\t" << "Trimming Info"<< endl;
    message << "\t\t" << "Number of loops: "<<trim.getNumOfLoops()<< endl;
    message << "\t\t" << "Outer loop index: "<<trim.getOutterLoopIndex()<< endl;
    message << "\t" << "Outter Loop["<<trim.getOutterLoopIndex()<<"]"<< endl;
    message << trim.getOutterLoop();
	for(int i=0; i<trim.getNumOfLoops();i++) {
		if(i==trim.getOutterLoopIndex()) continue;
    	message << "\t" << "InnerLoop["<<i<<"]"<< endl;
		message << trim.getLoop(i);
	}
    message << "\t" << "---------------------------------End Trimming" << endl;
	return message;
}

Message &operator<<(Message &message, const IGAPatchSurfaceTrimmingLoop &trim) {
    /// output loop
    for(int i=0;i<trim.getIGACurve().size();++i) {
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
