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

#ifndef GEOMETRYMATH_H_
#define GEOMETRYMATH_H_

#include <fstream>
#include <vector>
#include <cstdlib>
#include <map>
#include <list>
#include <vector>
#include <assert.h>
#include <typeinfo>
#include "AuxiliaryParameters.h"
#include "MatrixVectorMath.h"

namespace EMPIRE {
namespace MathLibrary {

// Methods
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double distanceLinePlane(double* Pline,double* Uline, double* Pplane,double* Nplane);
double distancePointSegment(double* _P, double* _P1, double* _P2);
// See http://paulbourke.net/geometry/pointlineplane/
double distanceLineLine(double& _ratioA, double& _ratioB, double* _P1, double* _P2,double* P3, double* P4);

/***********************************************************************************************
 * \brief project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \brief The result is the plane which has smallest angle with the unitNormal
 * \param[in] planeNormal normal vector of the plane
 * \return plane to project (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \author Tianyang Wang
 ***********/
int computePlaneToProject(const double *planeNormal);

/***********************************************************************************************
 * \brief Compute the center of a polygon
 * \param[in] polygon the polygon
 * \param[in] num num of points of the polygon
 * \param[out] center the center
 * \author Tianyang Wang
 ***********/
void computePolygonCenter(const double *polygon, int num, double *center);

/***********************************************************************************************
 * \brief Compute polygon area (points should be on the same plane)
 * \param[in] polygon the polygon
 * \param[in] num num of points of the polygon
 * \return the area of the polygon
 * \author Tianyang Wang
 ***********/
double computePolygonArea(const double *polygon, int num);

/***********************************************************************************************
 * \brief Project a number of points to a plane
 * \param[in] pointOnPlane a point on the plane
 * \param[in] unitNormal unit normal of the plane
 * \param[in] points the points to be projected
 * \param[in] num number of points to be projected
 * \param[in] projections the projections
 * \author Tianyang Wang
 ***********/
void projectToPlane(const double *pointOnPlane, const double *unitNormal, const double *points,
        int num, double *projections);

/***********************************************************************************************
 * \brief Compute the normal vector of a triangle
 * \param[in] triangle the triangle
 * \param[in] unitLength whether set the length of the normal vector to 1 or not
 * \param[out] normal normal vector
 * \author Tianyang Wang
 ***********/
void computeNormalOfTriangle(const double *triangle, bool unitLength, double *normal);

/***********************************************************************************************
 * \brief Compute the normal vector of a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] unitLength whether set the length of the normal vector to 1 or not
 * \param[out] normal normal vector
 * \author Tianyang Wang
 ***********/
void computeNormalOfQuad(const double *quad, bool unitLength, double *normal);

/***********************************************************************************************
 * \brief Compute the area of a triangle
 * \param[in] triangle the triangle
 * \return the area
 * \author Tianyang Wang
 ***********/
double computeAreaOfTriangle(const double *triangle);

/***********************************************************************************************
 * \brief Compute the length of a vector
 * \param[in] vec the vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength(const double *vec);

/***********************************************************************************************
 * \brief Compute the length of a 2D vector
 * \param[in] vec the 2D vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength2D(const double *vec);

/***********************************************************************************************
 * \brief Compute the cross product of two vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \param[out] crossProduct the crossProduct
 * \author Tianyang Wang
 ***********/
void computeVectorCrossProduct(const double *vec1, const double *vec2, double *crossProduct);

/***********************************************************************************************
 * \brief Compute the dot product of two vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \return dot product
 * \author Tianyang Wang
 ***********/
double computeVectorDotProduct(const double *vec1, const double *vec2);

/***********************************************************************************************
 * \brief Compute the square of the distance between two points
 * \param[in] p1 the 1st point
 * \param[in] p2 the 2nd point
 * \return the square of the distance between p1 and p2
 * \author Tianyang Wang
 ***********/
double distanceSquare(const double *p1, const double *p2);

/***********************************************************************************************
 * \brief Copy a point
 * \param[in] origin to be copied
 * \param[out] copy the copy
 * \author Tianyang Wang
 ***********/
void copyPoint(const double *origin, double *copy);

/***********************************************************************************************
 * \brief Copy a point
 * \param[in] origin to be copied
 * \param[out] copy the copy
 * \author Tianyang Wang
 ***********/
void copyElem(const double *origin, int size, double *copy);

/***********************************************************************************************
 * \brief Build a triangle
 * \param[in] p0 the 1st point
 * \param[in] p1 the 2nd point
 * \param[in] p2 the 3rd point
 * \param[out] triangle the triangle built by the three pointss
 * \author Tianyang Wang
 ***********/
void buildTrianagle(const double *p0, const double *p1, const double *p2, double *triangle);

/***********************************************************************************************
 * \brief Compute the longest edge length of the element
 * \param[in] elem the element
 * \param[in] size 3--->triangle, 4--->quadrilateral
 * \return the square of the longest length of the element
 * \author Tianyang Wang
 ***********/
double longestEdgeLengthSquare(const double *elem, int size);

/***********************************************************************************************
 * \brief Returns the intersection between a line and a triangle
 * \param[out] Flag whether the intersection has been found
 * \param[in] _coordsTriangle Coordinates of the triangle
 * \param[in] _coordsNode Coordinates of the point
 * \param[in] _direction Direction of the line
 * \param[in/out] _localCoords Local coordinates of the point
 * \author Chenshen Wu
 ***********/
bool computeIntersectionBetweenLineAndTriangle(const double *_coordsTriangle, const double* _coordsNode, const double* _direction, double* _localCoords);

/***********************************************************************************************
 * \brief Returns the intersection between a line and a quadrilateral
 * \param[out] Flag whether the intersection has been found
 * \param[in] _coordsQuad Coordinates of the quadrilateral
 * \param[in] _coordsNode Coordinates of the point
 * \param[in] _direction Direction of the line
 * \param[in/out] _localCoords Local coordinates of the point
 * \author Chenshen Wu
 ***********/
bool computeIntersectionBetweenLineAndQuad(const double *_coordsQuad, const double* _coordsNode, const double* _direction, double* _localCoords);

/***********************************************************************************************
 * \brief Computes the area of a given triangle defined by two vectors
 * \param[in] _x1, _y1, _z1 Vector 1
 * \param[in] _x2, _y2, _z2 Vector 2
 * \param[out] The area of the triangle included by the two vectors
 * \author Chenshen Wu
 ***********/
double computeAreaTriangle(double x1, double y1, double z1, double x2, double y2, double z2);

/***********************************************************************************************
 * \brief Computes the distance between two points
 * \param[in] _x1 Point 1
 * \param[in] _x2 Point 2
 * \param[out] The distance between the two given points
 * \author Chenshen Wu
 ***********/
double computePointDistance(double* _x1, double* _x2);


/***********************************************************************************************
 * \brief Clean the input polygon, remove colinear or identical point
 * \param[in/out] _polygon, vector of points coordinates
 * \author Fabien Pean
 ***********/
void cleanPolygon(std::vector<double>& _polygon);
void cleanPolygon(std::vector<std::pair<double,double> >& polygon);
void cleanPolygon(std::vector<std::pair<double,double> >& polygon,std::vector<std::pair<double,double> >& polygonSlave);



// Classes
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/***********************************************************************************************
 * \brief This class is used to cooperate with the clipping algorithm
 *        The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice,
 *        2nd edition, Addison-Wesley" P237 - P240.
 *        In practice this algorithm can have any strictly convex polygon as clipping window.
 *        This means consecutive points are not authorized to be colinear
 * ***********/
class PolygonClipper {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _polygonWindow the polygon that clips other polygons
     * \param[in] _sizePolygonWindow number of nodes/edges in this polygon
     * \param[in] _planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
     * \author Tianyang Wang
     ***********/
    PolygonClipper(const double *_polygonWindow, int _sizePolygonWindow, int _planeToProject);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~PolygonClipper();
    /***********************************************************************************************
     * \brief Clip.
     * \param[in] polygonToBeClipped polygon to be clipped by the window polygon
     * \param[in] sizePolygonToBeClipped size of the polygon
     * \param[out] the polygonResult results from the clipping
     * \return if there is a clip, return true, otherwise, return false
     * \author Tianyang Wang
     ***********/
    bool clip(const double *polygonToBeClipped, int sizePolygonToBeClipped,
            std::vector<double*> *polygonResult);
    /***********************************************************************************************
     * \brief Compute intersection between two lines. The result is the intersection with tolerance,
     *        which means it may locate outside of both lines. This error should be taken into account
     *        in the future. Make it static so that it is easily tested by some unit test.
     * \param[in] la0 the 1st point on la
     * \param[in] la1 the 2nd point on la
     * \param[in] lb0 the 1st point on lb
     * \param[in] lb1 the 2nd point on lb
     * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
     * \param[out] intersection the intersection of la and lb
     * \return boolean saying whether the intersection is on lb or not (true of false)
     * \author Tianyang Wang
     ***********/
    static bool intersect(const double *la0, const double *la1, const double *lb0, const double *lb1,
            int planeToProject, double *intersection);

private:
    /// the polygon that clips other polygons
    double *polygonWindow;
    /// number of edges of the polygon
    int sizePolygonWindow;
    /// project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
    const int planeToProject;
    /// flags used to determine inside/outside
    bool *insideFlag;
    /// slopes of all edges
    double *edgeSlopes;
    /// sign of whether to use y=f(x) or x=f(y)
    bool *reverseXY;
    // put the magic numbers here, if the number is too small, numerical error is high; if the number is
    // too large, computation error is high. Therefore, set the suitable value by experiment
    /// epsilon 1
    static const double EPS1;
    /// epsilon 2
    static const double EPS2;
    /// unit test class
    friend class TestMortarMath;
    /***********************************************************************************************
     * \brief Decides whether a point is on the "inside" side of the i-th edge of the polygon
     * \param[in] edgeID id of the edge
     * \param[in] point the point
     * \return if inside, return true, otherwise, return false
     * \author Tianyang Wang
     ***********/
    bool inside(int edgeID, double *point);
};


/***********************************************************************************************
 * \brief This class is used to clip fluid element on the IGA patch. The fluid nodes are already
 *  				projected to the IGA in the previous step.
 *        The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice,
 *        2nd edition, Addison-Wesley" P237 - P240.
 * ***********/
class IGAPolygonClipper {

private:
    /// Use the polygon clipper in MortarMath temporarily
    EMPIRE::MathLibrary::PolygonClipper *clipper;

    /// Instersection area between the polygon and window
    double polygonWindow[12];

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _u1, _u2, _v1, _v2 define the IGA element that clips other the fluid element.
     * \author Chenshen Wu
     ***********/
    IGAPolygonClipper(const double _u1, const double _u2, const double _v1, const double _v2);
    /***********************************************************************************************
     * \brief Destructor
     * \author Chenshen Wu
     ***********/
    virtual ~IGAPolygonClipper();
    /***********************************************************************************************
     * \brief Clip.
     * \param[in] polygonToBeClipped polygon to be clipped by the window polygon
     * \param[in] sizePolygonToBeClipped size of the polygon
     * \param[out] the polygonResult results from the clipping
     * \return if there is a clip, return true, otherwise, return false
     * \author Tianyang Wang
     ***********/
    bool clip(const double *_polygonToBeClipped, int _sizePolygonToBeClipped,
            std::vector<double*> *_polygonResult);

    /// unit test class
    friend class TestIGAMortarMath;

};




}
}




#endif /* GEOMETRYMATH_H_ */
