/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
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

#ifndef CLIPPERINTERFACE_H_
#define CLIPPERINTERFACE_H_

#include "clipper.hpp"

namespace EMPIRE {

class ClipperInterface {
public:
	typedef enum Operation{INTERSECTION=0, UNION, DIFFERENCE, XOR} Operation;
public:
    /***********************************************************************************************
     * \brief Constructors
     * \author Fabien Pean
     ***********/
	ClipperInterface();
	ClipperInterface(double _accuracy);
	ClipperInterface(double _accuracy, Operation _operation);
    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean
     ***********/
	virtual ~ClipperInterface();
    /***********************************************************************************************
     * \brief Set the desired clipping operation
     * \param[in] _operation	The operation as defined by the public enum Operation
     * \author Fabien Pean
     ***********/
	void setOperation(Operation _operation);
    /***********************************************************************************************
     * \brief Set the desired accuracy and thus factor
     * \param[in] _accuracy	The accuracy for the desired clipping
     * \author Fabien Pean
     ***********/
	void setAccuracy(double _accuracy);
    /***********************************************************************************************
     * \brief Retrieve the inner solution of the clipper interface
     * \param[in/out] _container	The container where solution is output
     * \author Fabien Pean
     ***********/
	void getSolution(std::vector<std::vector<double> >& _container);
	void getSolution(std::vector<std::vector<std::pair<double,double> > >& _container);
    /***********************************************************************************************
     * \brief Add a clipping window polygon (called path in clipper)
     * \param[in/out] _path	The polygon/curve to add to the clipping window
     * \author Fabien Pean
     ***********/
	void addPathClipper(const std::vector<double>& _path);
	void addPathClipper(const std::vector<std::pair<double,double> >& _path);
    /***********************************************************************************************
     * \brief Add a polygon to be clipped (called path in clipper)
     * \param[in/out] _path	The polygon/curve to add  to the "to be clipped"
     * \author Fabien Pean
     ***********/
	void addPathSubject(const std::vector<double>& _path);
	void addPathSubject(const std::vector<std::pair<double,double> >& _path);
    /***********************************************************************************************
     * \brief Execute the clipping with the clipping window member on the subject member. Results stored in solution
     * \author Fabien Pean
     ***********/
	void clip();
    /***********************************************************************************************
     * \brief Execute a standard intersection clipping by providing everything at once : subject/clipping window/solution
     * \author Fabien Pean
     ***********/
	void clip(int _numNodesPolygonToClip,
			  double* _nodesPolygonToClip,
			  int _numNodesClipper,
			  double* _nodesClipper,
			  int& _numNodesOutputPolygon,
			  double*& _nodesOutputPolygon);

	void clip(const std::vector<double>& _nodesPolygonToClip,
			  const std::vector<double>& _nodesClipper,
			  std::vector<double>& _nodesOutputPolygon);

	std::vector<double> clip(const std::vector<double>& _nodesPolygonToClip,
			                 const std::vector<double>& _nodesClipper);

	void clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip,
			  const std::vector<std::pair<double,double> >& _nodesClipper,
			  std::vector<std::pair<double,double> >& _nodesOutputPolygon);

	std::vector<std::pair<double,double> > clip(const std::vector<std::pair<double,double> >& _nodesPolygonToClip,
			                                    const std::vector<std::pair<double,double> >& _nodesClipper);
    /***********************************************************************************************
     * \brief Static call to the inner clean of clipper library
     * \author Fabien Pean
     ***********/
	static void cleanPolygon(std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);
	static bool isCounterclockwise(const std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);
	static void reversePolygon(std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);

	inline double getAccuracy() { return accuracy; }
private:
	/// Adaptee
	ClipperLib::Clipper clipper;
	/// Inner clipper polygons
	ClipperLib::Paths clipWindow;
	ClipperLib::Paths subject;
	ClipperLib::Paths solution;
	/// Parameters
	double accuracy;
	double factor;
	ClipperLib::ClipType operation;
};

} /* namespace EMPIRE */
#endif /* CLIPPERINTERFACE_H_ */
