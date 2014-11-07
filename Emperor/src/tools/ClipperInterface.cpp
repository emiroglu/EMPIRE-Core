/*
 * ClipperInterface.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: fabien
 */

#include "ClipperInterface.h"
#include "Message.h"
#include "assert.h"

using namespace ClipperLib;
using namespace std;


namespace EMPIRE {

ClipperInterface::ClipperInterface():
	accuracy(1e-9),
	factor(1e9),
	operation(ctIntersection) {
}

ClipperInterface::ClipperInterface(double _accuracy):
		accuracy(_accuracy),
		factor(1/_accuracy),
		operation(ctIntersection) {
}

ClipperInterface::ClipperInterface(double _accuracy, Operation _operation):
	accuracy(_accuracy),
	factor(1.0/_accuracy) {
	setOperation(_operation);
}

ClipperInterface::~ClipperInterface() {
}

void ClipperInterface::setAccuracy(double _accuracy) {
		accuracy=_accuracy;
		factor=1.0/accuracy;
}

void ClipperInterface::setOperation(Operation _operation) {
	switch(_operation) {
	case INTERSECTION : operation=ctIntersection; break;
	case UNION : operation=ctUnion; break;
	case DIFFERENCE : operation=ctDifference; break;
	case XOR : operation=ctXor; break;
	}
}
void ClipperInterface::getSolution(std::vector<std::vector<double> >& _container){
	_container.resize(solution.size());
	for(int i = 0; i < solution.size(); i++) {
		_container[i].resize(2*solution[i].size());
		for(int p = 0; p < solution[i].size(); p++) {
			_container[i][2*p]	= solution[i][p].X / factor;
			_container[i][2*p+1]= solution[i][p].Y / factor;
		}
	}
}
void ClipperInterface::getSolution(std::vector<std::vector<std::pair<double,double> > >& _container) {
	_container.resize(solution.size());
	for(int i = 0; i < solution.size(); i++) {
		_container[i].resize(solution[i].size());
		for(int p = 0; p < solution[i].size(); p++) {
			_container[i][p].first	= solution[i][p].X / factor;
			_container[i][p].second = solution[i][p].Y / factor;
		}
	}
}


void ClipperInterface::addPathClipper(const std::vector<double>& _path) {
	Path clip;
	for(int p=0; p < _path.size()/2; p++) {
		clip<<IntPoint((cInt)(_path[2*p]*factor),(cInt)(_path[2*p+1]*factor));
	}
	clipWindow.push_back(clip);
}
void ClipperInterface::addPathClipper(const std::vector<std::pair<double,double> >& _path) {
	Path clip;
	for(int p=0; p < _path.size(); p++) {
		clip<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	clipWindow.push_back(clip);
}
void ClipperInterface::addPathSubject(const std::vector<double>& _path) {
	Path subj;
	for(int p=0; p < _path.size()/2; p++) {
		subj<<IntPoint((cInt)(_path[2*p]*factor),(cInt)(_path[2*p+1]*factor));
	}
	subject.push_back(subj);
}
void ClipperInterface::addPathSubject(const std::vector<std::pair<double,double> >& _path) {
	Path subj;
	for(int p=0; p < _path.size(); p++) {
		subj<<IntPoint((cInt)(_path[p].first*factor),(cInt)(_path[p].second*factor));
	}
	subject.push_back(subj);
}


void ClipperInterface::clip() {
	assert(clipper.AddPaths(subject, ptSubject, true)==true);
	assert(clipper.AddPaths(clipWindow, ptClip, true)==true);
	assert(clipper.Execute(operation, solution, pftNonZero, pftNonZero)==true);
}

void ClipperInterface::clip(int _numNodesPolygonToClip, double* _nodesPolygonToClip, int _numNodesClipper,double* _nodesClipper, int& _numNodesOutputPolygon, double*& _nodesOutputPolygon) {
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

void ClipperInterface::clip(const std::vector<double>& _nodesPolygonToClip, const std::vector<double>& _nodesClipper,std::vector<double>& _nodesOutputPolygon) {
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

std::vector<double> ClipperInterface::clip(const std::vector<double>& _nodesPolygonToClip, const std::vector<double>& _nodesClipper) {
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

void ClipperInterface::clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip, const std::vector<std::pair<double,double> >& _nodesClipper, std::vector<std::pair<double,double> >& _nodesOutputPolygon) {
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

std::vector<std::pair<double,double> > ClipperInterface::clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip, const std::vector<std::pair<double,double> >& _nodesClipper) {
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

void ClipperInterface::cleanPolygon(std::vector<std::pair<double,double> >& _path, double _accuracy) {
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



} /* namespace EMPIRE */
