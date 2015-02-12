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

#ifndef CLIPPERADAPTER_H_
#define CLIPPERADAPTER_H_

#include "clipper.hpp"

namespace EMPIRE {

/********//**
 * \brief class ClipperInterface is a wrapper to the library clipper made by Angus Johnson @ http://www.angusj.com/delphi/clipper.php
 ***********/
class ClipperAdapter {
public:
    /***********************************************************************************************
     * \brief Enum Operation defines the type of boolean operation required
     ************/
	typedef enum Operation{INTERSECTION=0, UNION, DIFFERENCE, XOR} Operation;
    /***********************************************************************************************
     * \brief Enum Filling defines the type of polygon description
     * 			Details @ http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Types/PolyFillType.htm
     ************/
	typedef enum Filling{EVENODD=0, NONZERO, POSITIVE, NEGATIVE}Filling;
public:
    /***********************************************************************************************
     * \brief Constructors
     * \author Fabien Pean
     ***********/
	ClipperAdapter();
	ClipperAdapter(double _accuracy);
	ClipperAdapter(double _accuracy, Operation _operation);
	ClipperAdapter(double _accuracy, Operation _operation, Filling _filling);
	ClipperAdapter(double _accuracy, Operation _operation, Filling _fillingClipWindow, Filling _fillingSubject);
    /***********************************************************************************************
     * \brief Destructor
     * \author Fabien Pean
     ***********/
	virtual ~ClipperAdapter();
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
     * \brief Set the filling type of the polygon
     * 			See http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Types/PolyFillType.htm
     * \param[in] _filling The desired filling type of the polygon
     * \param[in] _subject If the value is to be set for subject 1 or clipping window 0, -1 for both.
     * \author Fabien Pean
     ***********/
	void setFilling(Filling _filling, int _subject=-1);
    /***********************************************************************************************
     * \brief Set the property strictly simple in clipper
     * 			See http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Classes/Clipper/Properties/StrictlySimple.htm
     * \param[in] _simple Force the parameter if set to 1
     * \author Fabien Pean
     ***********/
	void setStrictlySimple(bool _simple){clipper.StrictlySimple(_simple);};
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
	static void cleanPolygon(std::vector<double>& _path, double _accuracy=1e-9);
    /***********************************************************************************************
     * \brief Static call to the function PointInpolygon in clipper library
     * 			Returns 0 if false, -1 if pt is on poly and +1 if pt is in poly.
     * \author Fabien Pean
     ***********/
	static int isPointIn(const std::pair<double,double>& _point, const std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);
    /***********************************************************************************************
     * \brief Static call to the function Orientation in clipper library
     * 			Returns 0 if clockwise (antitrigonometric), 1 if counter clockwise (trigonometric)
     * \author Fabien Pean
     ***********/
	static bool isCounterclockwise(const std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);
    /***********************************************************************************************
     * \brief Static call to the function ReversePath in clipper library
     * \author Fabien Pean
     ***********/
	static void reversePolygon(std::vector<std::pair<double,double> >& _path, double _accuracy=1e-9);
    /***********************************************************************************************
     * \brief Static call to the function SimplifyPath in clipper library
     * \author Fabien Pean
     ***********/
	static void simplifyPolygon(std::vector<std::vector<std::pair<double,double> > >& _paths, double _accuracy=1e-9);

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
	ClipperLib::PolyFillType fillingSubject;
	ClipperLib::PolyFillType fillingClipWindow;

};

} /* namespace EMPIRE */
#endif /* CLIPPERINTERFACE_H_ */
