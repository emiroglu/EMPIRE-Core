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

double IGAPatchSurfaceTrimmingLoop::TOL_LINEARIZATION = 1E-6;

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

void IGAPatchSurfaceTrimmingLoop::linearize() {
//    linearizeUsingGreville();
//    linearizeUsingNCPxP();
    linearizeCombined();
}

void IGAPatchSurfaceTrimmingLoop::linearizeUsingGreville() {
    /*
     * Linear approximation of the nurbs curves
     *
     * Function layout :
     *
     * 1. For every NURBS curve
     * 1.1. For every control points of the curve
     * 1.1.1. Compute Greville abscissae
     * 1.1.2. Compute position in parametric space at Greville abscissae
     * 1.1.3. Store point in data structure of polylines
     * 2. Clean output polygon
     */
    for(int j=0;j<IGACurves.size();j++) {
		/// Check direction to put points in the right sequence (counter clockwise for outter loop, clockwise for inner)
		if(direction[j]) {
            for(int cpIndex=0;cpIndex<getIGACurve(j).getNoControlPoints();cpIndex++) {
				double knotGreville=getIGACurve(j).getIGABasis()->computeGrevilleAbscissae(cpIndex);
				double parametricCoordinates[2] = {0};
				getIGACurve(j).computeCartesianCoordinates(parametricCoordinates,knotGreville);
                polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
                IGACurves.at(j)->addPolylineVertex(parametricCoordinates[0],parametricCoordinates[1],knotGreville);
			}
		} else {
            for(int cpIndex=getIGACurve(j).getNoControlPoints()-1;cpIndex>=0;cpIndex--) {
				double knotGreville=getIGACurve(j).getIGABasis()->computeGrevilleAbscissae(cpIndex);
				double parametricCoordinates[2] = {0};
				getIGACurve(j).computeCartesianCoordinates(parametricCoordinates,knotGreville);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
                IGACurves.at(j)->addPolylineVertex(parametricCoordinates[0],parametricCoordinates[1],knotGreville);
			}
		}
	}
	ClipperAdapter::cleanPolygon(polylines);
}

void IGAPatchSurfaceTrimmingLoop::linearizeUsingNCPxP() {
    /*
     * Linear approximation of the nurbs curves
     *
     * Function layout :
     *
     * 1. For every NURBS curve
     * 1.1. Prepare data
     * 1.2. Compute knot delta
     * 1.1. For (nCP * p) nodes
     * 1.1.1. Compute knot
     * 1.1.2. Compute position in parametric space at knot
     * 1.1.3. Store point in data structure of polylines
     * 2. Clean the output polygon
     */
    for(int j=0;j<IGACurves.size();j++) {
		int nCP = getIGACurve(j).getNoControlPoints();
		int p = getIGACurve(j).getIGABasis()->getPolynomialDegree();
		double u0 = getIGACurve(j).getIGABasis()->getFirstKnot();
		double u1 = getIGACurve(j).getIGABasis()->getLastKnot();
		double du = (u1-u0)/(nCP*p-1);
		/// Check direction to put points in the right sequence (counter clockwise for outter loop, clockwise for inner
		if(direction[j]) {
			for(int i=0;i<nCP*p;i++) {
				double knot = u0 + i*du;
				double parametricCoordinates[2] = {0};
				getIGACurve(j).computeCartesianCoordinates(parametricCoordinates,knot);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
                IGACurves.at(j)->addPolylineVertex(parametricCoordinates[0],parametricCoordinates[1],knot);
            }
		} else {
			for(int i=nCP*p-1;i>=0;i--) {
				double knot = u0 + i*du;
				double parametricCoordinates[2] = {0};
				getIGACurve(j).computeCartesianCoordinates(parametricCoordinates,knot);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
                IGACurves.at(j)->addPolylineVertex(parametricCoordinates[0],parametricCoordinates[1],knot);
			}
		}
	}
	ClipperAdapter::cleanPolygon(polylines);
}

void IGAPatchSurfaceTrimmingLoop::linearizeCombined() {
    /*
     * Linear approximation of the nurbs curves
     *
     * Function layout :
     *
     * 1. For every NURBS curve
     * 1.1. Prepare data for NCPxP method
     * 1.2. Compute knot delta
     * 1.3. For every control points of the curve
     * 1.3.1. Compute and store Greville abscissae
     * 1.4. For (nCP * p) nodes
     * 1.4.1 Compute knot
     * 1.4.2 Find the position in the stored knots where the value is in correct order
     * 1.4.3 Insert the knot if it is not lying in the vicinity of an existing knot with a tolerance=1e-6
     * 1.1.1. Compute knot
     * 1.5. For every computed knot
     * 1.5.1 Store knot in uTildeCurve
     * 1.5.2 Compute position in parametric space
     * 1.5.3 Store point in data structure of polylines and in data structure of polyline inside IGAPatchCurve
     * 2. Clean the output polygon of trimming loop
     */

    for(int j=0;j<IGACurves.size();j++) {
        // NCPxP variables
        int nCP = getIGACurve(j).getNoControlPoints();
        int p = getIGACurve(j).getIGABasis()->getPolynomialDegree();
        double u0 = getIGACurve(j).getIGABasis()->getFirstKnot();
        double u1 = getIGACurve(j).getIGABasis()->getLastKnot();
        double du = (u1-u0)/(nCP*p-1);
        std::vector<double> uTildeCurve;
        /// Check direction to put points in the right sequence (counter clockwise for outter loop, clockwise for inner
        if(direction[j]) {

            // Create the linearization with Greville Abscissae and store the corresponding curve parameters
            for(int cpIndex=0;cpIndex<getIGACurve(j).getNoControlPoints();cpIndex++) {
                double knotGreville=getIGACurve(j).getIGABasis()->computeGrevilleAbscissae(cpIndex);
                uTildeCurve.push_back(knotGreville);
            }

            // Create the linearization with NCPxP and store the corresponding curve parameters considering the ordering
            // NCPxP
            std::vector<double>::iterator iUTildeCurve;
            for(int i=1;i<nCP*p-1;i++) {    // loop excludes the start and the end points
                iUTildeCurve = uTildeCurve.begin();
                double knot = u0 + i*du;

                while (knot>*iUTildeCurve)  iUTildeCurve++;

                if (fabs(*iUTildeCurve-knot)>TOL_LINEARIZATION){
                    uTildeCurve.insert(iUTildeCurve,knot);
                    iUTildeCurve++;
                }
             }
        } else {

            // Create the linearization with Greville Abscissae and store the corresponding curve parameters
            for(int cpIndex=getIGACurve(j).getNoControlPoints()-1;cpIndex>=0;cpIndex--) {
                double knotGreville=getIGACurve(j).getIGABasis()->computeGrevilleAbscissae(cpIndex);
                uTildeCurve.push_back(knotGreville);
            }

            // Create the linearization with NCPxP and store the corresponding curve parameters considering the ordering
            // NCPxP
            std::vector<double>::iterator iUTildeCurve;
            for(int i=nCP*p-2;i>0;i--) {    // loop excludes the start and the end points
                iUTildeCurve = uTildeCurve.begin();
                double knot = u0 + i*du;

                while (knot<*iUTildeCurve)  iUTildeCurve++;

                if (fabs(*iUTildeCurve-knot)>TOL_LINEARIZATION){
                    uTildeCurve.insert(iUTildeCurve,knot);
                    iUTildeCurve++;
                }
             }
        }

        double parametricCoordinates[2] = {0.0, 0.0};
        for (int knotCtr=0; knotCtr<uTildeCurve.size(); knotCtr++){
            getIGACurve(j).computeCartesianCoordinates(parametricCoordinates, uTildeCurve[knotCtr]);
            polylines.push_back(parametricCoordinates[0]);
            polylines.push_back(parametricCoordinates[1]);
            IGACurves.at(j)->addPolylineVertex(parametricCoordinates[0], parametricCoordinates[1], uTildeCurve[knotCtr]);
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
