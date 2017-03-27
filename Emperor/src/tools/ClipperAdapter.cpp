/*  Copyright &copy; 2014, TU Muenchen, Chair of Structural Analysis,
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

#include "ClipperAdapter.h"
#include "Message.h"
#include "assert.h"

using namespace ClipperLib;
using namespace std;

namespace EMPIRE {

ClipperAdapter::ClipperAdapter():
	accuracy(1e-9),
	factor(1e9),
	operation(ctIntersection),
	fillingClipWindow(pftNonZero),
	fillingSubject(pftNonZero) {
}

ClipperAdapter::ClipperAdapter(double _accuracy):
		accuracy(_accuracy),
		factor(1/_accuracy),
		operation(ctIntersection),
		fillingClipWindow(pftNonZero),
		fillingSubject(pftNonZero) {
}

ClipperAdapter::ClipperAdapter(double _accuracy, Operation _operation):
	accuracy(_accuracy),
	factor(1.0/_accuracy),
	fillingClipWindow(pftNonZero),
	fillingSubject(pftNonZero) {
	setOperation(_operation);
}
ClipperAdapter::ClipperAdapter(double _accuracy, Operation _operation, Filling _filling):
	accuracy(_accuracy),
	factor(1.0/_accuracy) {
	setOperation(_operation);
	setFilling(_filling);
}

ClipperAdapter::ClipperAdapter(double _accuracy, Operation _operation, Filling _fillingClipWindow, Filling _fillingSubject):
	accuracy(_accuracy),
	factor(1.0/_accuracy) {
	setOperation(_operation);
	setFilling(_fillingClipWindow, 0);
	setFilling(_fillingSubject, 1);
}

	ClipperAdapter::~ClipperAdapter() {
}

void ClipperAdapter::setAccuracy(double _accuracy) {
		accuracy=_accuracy;
		factor=1.0/accuracy;
}

void ClipperAdapter::setFilling(Filling _filling, int _subject) {
	ClipperLib::PolyFillType filling;
	switch(_filling) {
	case EVENODD : 	filling=pftEvenOdd; break;
	case NONZERO : 	filling=pftNonZero; break;
	case POSITIVE : filling=pftPositive; break;
	case NEGATIVE : filling=pftNegative; break;
	}
	switch(_subject) {
	case 0 : fillingClipWindow=filling;break;
	case 1 : fillingSubject=filling;break;
	default :fillingSubject=filling; fillingClipWindow=filling; break;
	}
}

void ClipperAdapter::setOperation(Operation _operation) {
	switch(_operation) {
	case INTERSECTION : operation=ctIntersection; break;
	case UNION : operation=ctUnion; break;
	case DIFFERENCE : operation=ctDifference; break;
	case XOR : operation=ctXor; break;
	}
}

void ClipperAdapter::getSolution(std::vector<std::vector<double> >& _container){
	_container.resize(solution.size());
	for(int i = 0; i < solution.size(); i++) {
		_container[i].resize(2*solution[i].size());
		for(int p = 0; p < solution[i].size(); p++) {
			_container[i][2*p]	= solution[i][p].X / factor;
			_container[i][2*p+1]= solution[i][p].Y / factor;
		}
	}
}

void ClipperAdapter::getSolution(std::vector<std::vector<std::pair<double,double> > >& _container) {
	_container.resize(solution.size());
	for(int i = 0; i < solution.size(); i++) {
		_container[i].resize(solution[i].size());
		for(int p = 0; p < solution[i].size(); p++) {
			_container[i][p].first	= solution[i][p].X / factor;
			_container[i][p].second = solution[i][p].Y / factor;
		}
	}
}

void ClipperAdapter::addPathClipper(const std::vector<double>& _path) {
	Path clip;
	for(int p=0; p < _path.size()/2; p++) {
		clip<<IntPoint((cInt)(_path[2*p]*factor),(cInt)(_path[2*p+1]*factor));
	}
	clipWindow.push_back(clip);
}

void ClipperAdapter::addPathClipper(const std::vector<std::pair<double,double> >& _path) {
	Path clip;
	for(int p=0; p < _path.size(); p++) {
		clip<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	clipWindow.push_back(clip);
}

void ClipperAdapter::addPathSubject(const std::vector<double>& _path) {
	Path subj;
	for(int p=0; p < _path.size()/2; p++) {
		subj<<IntPoint((cInt)(_path[2*p]*factor),(cInt)(_path[2*p+1]*factor));
	}
	subject.push_back(subj);
}

void ClipperAdapter::addPathSubject(const std::vector<std::pair<double,double> >& _path) {
	Path subj;
	for(int p=0; p < _path.size(); p++) {
		subj<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	subject.push_back(subj);
}

void ClipperAdapter::clip() {
	assert(clipper.AddPaths(subject, ptSubject, true)==true);
	assert(clipper.AddPaths(clipWindow, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, fillingSubject, fillingClipWindow)==true);
}

void ClipperAdapter::clip(int _numNodesPolygonToClip, double* _nodesPolygonToClip, int _numNodesClipper,double* _nodesClipper, int& _numNodesOutputPolygon, double*& _nodesOutputPolygon) {
	Path subj,clip;
	Paths solution;
	for(int p=0; p < _numNodesPolygonToClip; p++) {
		subj<<IntPoint((cInt)(_nodesPolygonToClip[2*p]*factor),(cInt)(_nodesPolygonToClip[2*p+1]*factor));
	}
	for(int p=0;p < _numNodesClipper; p++) {
		clip<<IntPoint((cInt)(_nodesClipper[2*p]*factor),(cInt)(_nodesClipper[2*p+1]*factor));
	}
	assert(clipper.AddPath(subj, ptSubject, true)==true);
	assert(clipper.AddPath(clip, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, pftNonZero, pftNonZero)==true);

	if(solution.empty())
		_numNodesOutputPolygon = 0;
	else
		_numNodesOutputPolygon = solution[0].size();
	_nodesOutputPolygon=new double[2*_numNodesOutputPolygon];
	for(int p=0;p<_numNodesOutputPolygon;p++) {
		_nodesOutputPolygon[2*p] = solution[0][p].X / factor;
		_nodesOutputPolygon[2*p+1] = solution[0][p].Y / factor;
	}
}

void ClipperAdapter::clip(const std::vector<double>& _nodesPolygonToClip, const std::vector<double>& _nodesClipper,std::vector<double>& _nodesOutputPolygon) {
	Path subj,clip;
	Paths solution;
	for(int p=0; p < _nodesPolygonToClip.size()/2; p++) {
		subj<<IntPoint((cInt)(_nodesPolygonToClip[2*p]*factor),(cInt)(_nodesPolygonToClip[2*p+1]*factor));
	}
	for(int p=0;p < _nodesClipper.size()/2; p++) {
		clip<<IntPoint((cInt)(_nodesClipper[2*p]*factor),(cInt)(_nodesClipper[2*p+1]*factor));
	}

	assert(clipper.AddPath(subj, ptSubject, true)==true);
	assert(clipper.AddPath(clip, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, pftNonZero, pftNonZero)==true);

	int numNodesOutputPolygon;
	if(solution.empty())
		numNodesOutputPolygon = 0;
	else
		numNodesOutputPolygon = solution[0].size();	_nodesOutputPolygon.resize(2*numNodesOutputPolygon);
	for(int p=0;p<numNodesOutputPolygon;p++) {
		_nodesOutputPolygon[2*p] = solution[0][p].X / factor;
		_nodesOutputPolygon[2*p+1] = solution[0][p].Y / factor;
	}
}

std::vector<double> ClipperAdapter::clip(const std::vector<double>& _nodesPolygonToClip, const std::vector<double>& _nodesClipper) {
	Path subj,clip;
	Paths solution;
	for(int p=0; p < _nodesPolygonToClip.size()/2; p++) {
		subj<<IntPoint((cInt)(_nodesPolygonToClip[2*p]*factor),(cInt)(_nodesPolygonToClip[2*p+1]*factor));
	}
	for(int p=0;p < _nodesClipper.size()/2; p++) {
		clip<<IntPoint((cInt)(_nodesClipper[2*p]*factor),(cInt)(_nodesClipper[2*p+1]*factor));
	}

	assert(clipper.AddPath(subj, ptSubject, true)==true);
	assert(clipper.AddPath(clip, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, pftNonZero, pftNonZero)==true);

	int numNodesOutputPolygon;
	if(solution.empty())
		numNodesOutputPolygon = 0;
	else
		numNodesOutputPolygon = solution[0].size();	vector<double> nodesOutputPolygon(2*numNodesOutputPolygon);
	for(int p=0; p < numNodesOutputPolygon; p++) {
		nodesOutputPolygon[2*p] = solution[0][p].X / factor;
		nodesOutputPolygon[2*p+1] = solution[0][p].Y / factor;
	}
	return nodesOutputPolygon;
}

void ClipperAdapter::clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip, const std::vector<std::pair<double,double> >& _nodesClipper, std::vector<std::pair<double,double> >& _nodesOutputPolygon) {
	Path subj,clip;
	Paths solution;
	for(int p=0; p < _nodesPolygonToClip.size(); p++) {
		subj<<IntPoint((cInt)(_nodesPolygonToClip[p].first*factor),(cInt)(_nodesPolygonToClip[p].second*factor));
	}
	for(int p=0;p < _nodesClipper.size(); p++) {
		clip<<IntPoint((cInt)(_nodesClipper[p].first*factor),(cInt)(_nodesClipper[p].second*factor));
	}

	assert(clipper.AddPath(subj, ptSubject, true)==true);
	assert(clipper.AddPath(clip, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, pftNonZero, pftNonZero)==true);

	int numNodesOutputPolygon;
	if(solution.empty())
		numNodesOutputPolygon = 0;
	else
		numNodesOutputPolygon = solution[0].size();
	_nodesOutputPolygon.resize(numNodesOutputPolygon);
	for(int p=0; p < numNodesOutputPolygon; p++) {
		_nodesOutputPolygon[p].first = solution[0][p].X / factor;
		_nodesOutputPolygon[p].second = solution[0][p].Y / factor;
	}
}

std::vector<std::pair<double,double> > ClipperAdapter::clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip, const std::vector<std::pair<double,double> >& _nodesClipper) {
	Path subj,clip;
	Paths solution;
	for(int p=0; p < _nodesPolygonToClip.size(); p++) {
		subj<<IntPoint((cInt)(_nodesPolygonToClip[p].first*factor),(cInt)(_nodesPolygonToClip[p].second*factor));
	}
	for(int p=0;p < _nodesClipper.size(); p++) {
		clip<<IntPoint((cInt)(_nodesClipper[p].first*factor),(cInt)(_nodesClipper[p].second*factor));
	}

	assert(clipper.AddPath(subj, ptSubject, true)==true);
	assert(clipper.AddPath(clip, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, pftNonZero, pftNonZero)==true);

	int numNodesOutputPolygon;
	if(solution.empty())
		numNodesOutputPolygon = 0;
	else
		numNodesOutputPolygon = solution[0].size();
	vector<pair<double,double> > nodesOutputPolygon(numNodesOutputPolygon);
	for(int p=0; p < numNodesOutputPolygon; p++) {
		nodesOutputPolygon[p].first = solution[0][p].X / factor;
		nodesOutputPolygon[p].second = solution[0][p].Y / factor;
	}
	return nodesOutputPolygon;
}

int ClipperAdapter::isPointIn(const std::pair<double,double>& _point, const std::vector<std::pair<double,double> >& _path, double _accuracy) {
	double factor = 1 / _accuracy;
	Path subj;
	for(int p=0; p < _path.size(); p++) {
		subj<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	IntPoint point;
	point.X=_point.first*factor;
	point.Y=_point.second*factor;

	return ClipperLib::PointInPolygon(point,subj);
}

void ClipperAdapter::cleanPolygon(std::vector<std::pair<double,double> >& _path, double _accuracy) {
	double factor = 1 / _accuracy;
	Path subj;
	for(int p=0; p < _path.size(); p++) {
		subj<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	CleanPolygon(subj);
	int numNodesOutputPolygon = subj.size();
	_path.resize(numNodesOutputPolygon);
	for(int p=0; p < numNodesOutputPolygon; p++) {
		_path[p].first  = subj[p].X / factor;
		_path[p].second = subj[p].Y / factor;
	}
}
void ClipperAdapter::cleanPolygon(std::vector<double>& _path, double _accuracy) {
	double factor = 1 / _accuracy;
	Path subj;
	for(int p=0; p < _path.size()/2; p++) {
		subj<<IntPoint((cInt)(_path[2*p]*factor),(cInt)(_path[2*p+1]*factor));
	}
	CleanPolygon(subj);
	int numNodesOutputPolygon = subj.size();
	_path.resize(2*numNodesOutputPolygon);
	for(int p=0; p < numNodesOutputPolygon; p++) {
		_path[2*p]   = subj[p].X / factor;
		_path[2*p+1] = subj[p].Y / factor;
	}
}

bool ClipperAdapter::isCounterclockwise(const std::vector<std::pair<double,double> >& _path, double _accuracy) {
	double factor = 1 / _accuracy;
	Path subj;
	for(int p=0; p < _path.size(); p++) {
		subj<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	bool isCounterclockwise = Orientation(subj);
	return isCounterclockwise;
}
void ClipperAdapter::reversePolygon(std::vector<std::pair<double,double> >& _path, double _accuracy){
	double factor = 1 / _accuracy;
	Path subj;
	for(int p=0; p < _path.size(); p++) {
		subj<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	ReversePath(subj);
	for(int p=0; p < _path.size(); p++) {
		_path[p].first  = subj[p].X / factor;
		_path[p].second = subj[p].Y / factor;
	}
}

void ClipperAdapter::simplifyPolygon(std::vector<std::vector<std::pair<double,double> > >& _paths, double _accuracy) {
	double factor = 1 / _accuracy;
	Paths subjs;
	for(int s=0; s< _paths.size(); s++) {
		Path subj;
		for(int p=0; p < _paths[s].size(); p++) {
			subj<<IntPoint((cInt)(_paths[s][p].first*factor),(cInt)(_paths[s][p].second*factor));
		}
		subjs.push_back(subj);
	}
	SimplifyPolygons(subjs,pftNonZero);
	_paths.clear();
	_paths.resize(subjs.size());
	for(int s=0; s < subjs.size(); s++) {
		_paths[s].resize(subjs[s].size());
		for(int p=0; p < subjs[s].size(); p++) {
			_paths[s][p].first  = subjs[s][p].X / factor;
			_paths[s][p].second = subjs[s][p].Y / factor;
		}
	}
}

} /* namespace EMPIRE */
