/*
 * ClipperInterface.h
 *
 *  Created on: Nov 5, 2014
 *      Author: fabien
 */

#ifndef CLIPPERINTERFACE_H_
#define CLIPPERINTERFACE_H_

#include "clipper.hpp"

namespace EMPIRE {

class ClipperInterface {
public:
	typedef enum Operation{INTERSECTION=0, UNION, DIFFERENCE, XOR} Operation;
public:
	ClipperInterface();
	ClipperInterface(double _accuracy);
	ClipperInterface(double _accuracy, Operation _operation);
	virtual ~ClipperInterface();

	void setOperation(Operation _operation);
	void setAccuracy(double _accuracy);
	void getSolution(std::vector<std::vector<double> >& _container);
	void getSolution(std::vector<std::vector<std::pair<double,double> > >& _container);

	void addPathClipper(const std::vector<double>& _path);
	void addPathClipper(const std::vector<std::pair<double,double> >& _path);
	void addPathSubject(const std::vector<double>& _path);
	void addPathSubject(const std::vector<std::pair<double,double> >& _path);


	void clip();
	void clip(int _numNodesPolygonToClip, double* _nodesPolygonToClip, int _numNodesClipper,double* _nodesClipper, int& _numNodesOutputPolygon, double*& _nodesOutputPolygon);
	void clip(const std::vector<double>& _nodesPolygonToClip, const std::vector<double>& _nodesClipper,std::vector<double>& _nodesOutputPolygon);
	std::vector<double> clip(const std::vector<double>& _nodesPolygonToClip, const std::vector<double>& _nodesClipper);
	void clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip, const std::vector<std::pair<double,double> >& _nodesClipper, std::vector<std::pair<double,double> >& _nodesOutputPolygon);
	std::vector<std::pair<double,double> > clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip, const std::vector<std::pair<double,double> >& _nodesClipper);

	static void cleanPolygon(std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);

	inline double getAccuracy() { return accuracy; }
private:
	ClipperLib::Clipper clipper;
	ClipperLib::Paths clipWindow;
	ClipperLib::Paths subject;
	ClipperLib::Paths solution;

	double accuracy;
	double factor;
	ClipperLib::ClipType operation;
};

} /* namespace EMPIRE */
#endif /* CLIPPERINTERFACE_H_ */
