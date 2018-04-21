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

#include "GeometryMath.h"
#include "ConstantsAndVariables.h"
#include "DebugMath.h"

using namespace std;
namespace EMPIRE {
namespace MathLibrary {


// Methods
// %%%%%%%%%%%%%%%%%%%%%%%

double distancePointSegment(double* _P, double* _P1, double* _P2) {
	double distance[3];
	double P1P[3], PP2[3];
	double P1P2[3];
	double normP1P2;
	for(int i=0;i<3;i++){
		P1P[i]=_P[i]-_P1[i];
		PP2[i]=_P2[i]-_P[i];
		P1P2[i]=_P2[i]-_P1[i];
	}
	// Edit Aditya
	//normP1P2=sqrt(EMPIRE::MathLibrary::square2normVector(3,P1P2));
	normP1P2=EMPIRE::MathLibrary::vector2norm(P1P2,3);

    double projP1P = sqrt(EMPIRE::MathLibrary::computeDenseDotProduct(3,P1P,P1P2));
     if ( projP1P <= 0 ){
          //return sqrt(EMPIRE::MathLibrary::square2normVector(3,P1P));
     	  return EMPIRE::MathLibrary::vector2norm(P1P,3);
     }

     if ( normP1P2 <= projP1P ){
         //return sqrt(EMPIRE::MathLibrary::square2normVector(3,PP2));
         return EMPIRE::MathLibrary::vector2norm(PP2,3);
     }

	double t = projP1P / normP1P2;
	double tmp[3];
	for(int i=0;i<3;i++){
	     tmp[i]= _P1[i] + t * P1P2[i];
		 tmp[i]= _P[i] - tmp[i];
	}
    //return sqrt(EMPIRE::MathLibrary::square2normVector(3,tmp));
    return EMPIRE::MathLibrary::vector2norm(tmp,3);
}

double distanceLinePlane(double* Pline,double* Uline, double* Pplane,double* Nplane) {

	double denom=EMPIRE::MathLibrary::computeDenseDotProduct(3,Nplane,Uline);
	if(denom<1e-9) return -1;
	double PpPl[3];
	for(int i=0;i<3;i++){
		PpPl[i]=Pline[i]-Pplane[i];
	}
	return EMPIRE::MathLibrary::computeDenseDotProduct(3,PpPl,Nplane)/denom;
}

double distanceLineLine(double& _ratioA, double& _ratioB, double* _P1, double* _P2,double* _P3, double* _P4) {
	double distance[3];
	double P1P2[3],P1P3[3], P3P4[3];
	double normP1P2, normP3P4;
	for(int i=0;i<3;i++){
		P1P2[i]=_P2[i]-_P1[i];
		P1P3[i]=_P3[i]-_P1[i];
		P3P4[i]=_P4[i]-_P3[i];

	}
	//normP1P2=sqrt(EMPIRE::MathLibrary::square2normVector(3,P1P2));
	//normP3P4=sqrt(EMPIRE::MathLibrary::square2normVector(3,P3P4));
	normP1P2=EMPIRE::MathLibrary::vector2norm(P1P2,3);
	normP3P4=EMPIRE::MathLibrary::vector2norm(P3P4,3);

	if(normP3P4 < EPS)
	    return -1;
	if(normP1P2 < EPS)
	    return -1;

	double d13_34=EMPIRE::MathLibrary::computeDenseDotProduct(3,P1P3,P3P4);
	double d34_12=EMPIRE::MathLibrary::computeDenseDotProduct(3,P3P4,P1P2);
	double d13_12=EMPIRE::MathLibrary::computeDenseDotProduct(3,P1P3,P1P2);
	double d34_34=EMPIRE::MathLibrary::computeDenseDotProduct(3,P3P4,P3P4);
	double d12_12=EMPIRE::MathLibrary::computeDenseDotProduct(3,P1P2,P1P2);

	double numer,denom;
	denom = d12_12 * d34_34 - d34_12 * d34_12;
	if (fabs(denom) < EPS)
	  return -1;
	numer = d13_34 * d34_12 - d13_12 * d34_34;

	_ratioA = numer / denom;
	_ratioB = (d13_34 + d34_12 * _ratioA) / d34_34;

	double Pa[3], Pb[3];
	Pa[0] = _P1[0] + _ratioA * P1P2[0];
	Pa[1] = _P1[1] + _ratioA * P1P2[1];
	Pa[2] = _P1[2] + _ratioA * P1P2[2];
	Pb[0] = _P3[0] + _ratioB * P3P4[0];
	Pb[1] = _P3[1] + _ratioB * P3P4[1];
	Pb[2] = _P3[2] + _ratioB * P3P4[2];
	double PaPb[3];
	for(int i=0;i<3;i++){
		PaPb[i]=Pb[i]-Pa[i];
	}
	//return sqrt(EMPIRE::MathLibrary::square2normVector(3,PaPb));
	return EMPIRE::MathLibrary::vector2norm(PaPb,3);
}

/***********************************************************************************************
 * \brief project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \brief The result is the plane which has smallest angle with the unitNormal
 * \param[in] planeNormal normal vector of the plane
 * \return plane to project (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \author Tianyang Wang
 ***********/
int computePlaneToProject(const double *planeNormal) {
    // find out the biggest entry in planeNormal, which indicates the plane to project
    double max = fabs(planeNormal[0]);
    int imax = 0;
    for (int i = 1; i < 3; i++) {
        if (fabs(planeNormal[i]) > max) {
            max = fabs(planeNormal[i]);
            imax = i;
        }
    }
    return imax;
}

/***********************************************************************************************
 * \brief Compute the center of a polygon
 * \param[in] polygon the polygon
 * \param[in] num num of points of the polygon
 * \param[out] center the center
 * \author Tianyang Wang
 ***********/
void computePolygonCenter(const double *polygon, int num, double *center) {
    for (int i = 0; i < 3; i++)
        center[i] = 0;
    for (int i = 0; i < num; i++)
        for (int j = 0; j < 3; j++)
            center[j] += polygon[i * 3 + j];
    for (int i = 0; i < 3; i++)
        center[i] /= (double) num;
}

/***********************************************************************************************
 * \brief Compute polygon area (points should be on the same plane)
 * \param[in] polygon the polygon
 * \param[in] num num of points of the polygon
 * \return the area of the polygon
 * \author Tianyang Wang
 ***********/
double computePolygonArea(const double *polygon, int num) {
    double center[3];
    computePolygonCenter(polygon, num, center);

    double area = 0.0;
    for (int i = 0; i < num; i++) {
        double clipTriangle[9];
        buildTrianagle(center, &polygon[i * 3], &polygon[(i + 1) % num * 3], clipTriangle);
        area += computeAreaOfTriangle(clipTriangle);
    }
    return area;
}

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
        int num, double *projections) {
    double ridge[3];
    double distance;
    double path[3];
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < 3; j++)
            ridge[j] = pointOnPlane[j] - points[i * 3 + j];
        distance = computeVectorDotProduct(unitNormal, ridge);
        for (int j = 0; j < 3; j++)
            path[j] = distance * unitNormal[j];

        for (int j = 0; j < 3; j++)
            projections[i * 3 + j] = points[i * 3 + j] + path[j];
    }
}

/***********************************************************************************************
 * \brief Compute the normal vector of a triangle
 * \param[in] triangle the triangle
 * \param[in] unitLength whether set the length of the normal vector to 1 or not
 * \param[out] normal normal vector
 * \author Tianyang Wang
 ***********/
void computeNormalOfTriangle(const double *triangle, bool unitLength, double *normal) {
    const int NUM = 3; // three points in a triangle
    double length = 0; // length of the normal vector
    normal[0] = (triangle[2 * NUM + 1] - triangle[0 * NUM + 1])
            * (triangle[2 * NUM + 2] - triangle[1 * NUM + 2])
            - (triangle[2 * NUM + 1] - triangle[1 * NUM + 1])
                    * (triangle[2 * NUM + 2] - triangle[0 * NUM + 2]);
    normal[1] = -((triangle[2 * NUM + 0] - triangle[0 * NUM + 0])
            * (triangle[2 * NUM + 2] - triangle[1 * NUM + 2])
            - (triangle[2 * NUM + 0] - triangle[1 * NUM + 0])
                    * (triangle[2 * NUM + 2] - triangle[0 * NUM + 2]));
    normal[2] = (triangle[2 * NUM + 0] - triangle[0 * NUM + 0])
            * (triangle[2 * NUM + 1] - triangle[1 * NUM + 1])
            - (triangle[2 * NUM + 0] - triangle[1 * NUM + 0])
                    * (triangle[2 * NUM + 1] - triangle[0 * NUM + 1]);

    if (!unitLength)
        return;

    length = computeVectorLength(normal);

    for (int i = 0; i < 3; i++)
        normal[i] /= length;
}

/***********************************************************************************************
 * \brief Compute the normal vector of a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] unitLength whether set the length of the normal vector to 1 or not
 * \param[out] normal normal vector
 * \author Tianyang Wang
 ***********/
void computeNormalOfQuad(const double *quad, bool unitLength, double *normal) {
    // d_N_d_xi[4] contains the partial derivative w.r.t. xi of the shape functions at (0, 0)
    double d_N_d_xi[4];
    // d_N_d_eta[4] contains the partial derivative w.r.t. eta of the shape functions at (0, 0)
    double d_N_d_eta[4];

    d_N_d_xi[0] = -0.25;
    d_N_d_xi[1] = 0.25;
    d_N_d_xi[2] = 0.25;
    d_N_d_xi[3] = -0.25;

    d_N_d_eta[0] = -0.25;
    d_N_d_eta[1] = -0.25;
    d_N_d_eta[2] = 0.25;
    d_N_d_eta[3] = 0.25;
    // g1 and g2 are the local basis vectors, and det(J)=||g1 x g2||
    double g1[3] = { 0, 0, 0 };
    double g2[3] = { 0, 0, 0 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            g1[i] += quad[j * 3 + i] * d_N_d_xi[j];
            g2[i] += quad[j * 3 + i] * d_N_d_eta[j];
        }
    }
    double length;
    normal[0] = g1[1] * g2[2] - g2[1] * g1[2];
    normal[1] = -(g1[0] * g2[2] - g2[0] * g1[2]);
    normal[2] = g1[0] * g2[1] - g2[0] * g1[1];

    if (!unitLength)
        return;

    length = computeVectorLength(normal);

    for (int i = 0; i < 3; i++)
        normal[i] /= length;
}

/***********************************************************************************************
 * \brief Compute the area of a triangle
 * \param[in] triangle the triangle
 * \return the area
 * \author Tianyang Wang
 ***********/
double computeAreaOfTriangle(const double *triangle) {
    double area = 0.0;
    double normal[3];
    computeNormalOfTriangle(triangle, false, normal); // 2 times the length of the normal vector equals the area
    area = computeVectorLength(normal);
    area = 0.5 * area;
    return area;
}

/***********************************************************************************************
 * \brief Compute the length of a vector
 * \param[in] vec the vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength(const double *vec);

/***********************************************************************************************
 * \brief Compute the length of a vector
 * \param[in] vec the vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength(const double *vec) {
    double length = 0;
    for (int i = 0; i < 3; i++)
        length += vec[i] * vec[i];
    length = sqrt(length);
    return length;
}

/***********************************************************************************************
 * \brief Compute the length of a 2D vector
 * \param[in] vec the 2D vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength2D(const double *vec) {
    double length = 0;
    for (int i = 0; i < 2; i++)
        length += vec[i] * vec[i];
    length = sqrt(length);
    return length;
}

/***********************************************************************************************
 * \brief Compute the cross product of two vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \param[out] crossProduct the crossProduct
 * \author Tianyang Wang
 ***********/
void computeVectorCrossProduct(const double *vec1, const double *vec2, double *crossProduct) {
    crossProduct[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
    crossProduct[1] = -(vec1[0] * vec2[2] - vec2[0] * vec1[2]);
    crossProduct[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];
}

/***********************************************************************************************
 * \brief Compute the dot product of two vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \return dot product
 * \author Tianyang Wang
 ***********/
double computeVectorDotProduct(const double *vec1, const double *vec2) {
    double product = 0;
    for (int i = 0; i < 3; i++)
        product += vec1[i] * vec2[i];
    return product;
}

/***********************************************************************************************
 * \brief Compute the square of the distance between two points
 * \param[in] p1 the 1st point
 * \param[in] p2 the 2nd point
 * \return the square of the distance between p1 and p2
 * \author Tianyang Wang
 ***********/
double distanceSquare(const double *p1, const double *p2) {
    double d;
    d = (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
            + (p1[2] - p2[2]) * (p1[2] - p2[2]);
    return d;
}

/***********************************************************************************************
 * \brief Copy a point
 * \param[in] origin to be copied
 * \param[out] copy the copy
 * \author Tianyang Wang
 ***********/
void copyPoint(const double *origin, double *copy) {
    for (int i = 0; i < 3; i++)
        copy[i] = origin[i];
}

/***********************************************************************************************
 * \brief Copy a point
 * \param[in] origin to be copied
 * \param[out] copy the copy
 * \author Tianyang Wang
 ***********/
void copyElem(const double *origin, int size, double *copy) {
    for (int i = 0; i < 3 * size; i++)
        copy[i] = origin[i];
}

/***********************************************************************************************
 * \brief Build a triangle
 * \param[in] p0 the 1st point
 * \param[in] p1 the 2nd point
 * \param[in] p2 the 3rd point
 * \param[out] triangle the triangle built by the three pointss
 * \author Tianyang Wang
 ***********/
void buildTrianagle(const double *p0, const double *p1, const double *p2, double *triangle) {
    for (int i = 0; i < 3; i++) {
        triangle[i] = p0[i];
        triangle[3 + i] = p1[i];
        triangle[6 + i] = p2[i];
    }
}

/***********************************************************************************************
 * \brief Compute the longest edge length of the element
 * \param[in] elem the element
 * \param[in] size 3--->triangle, 4--->quadrilateral
 * \return the square of the longest length of the element
 * \author Tianyang Wang
 ***********/
double longestEdgeLengthSquare(const double *elem, int size) {
    double l[size];
    for (int i = 0; i < size; i++) {
        l[i] = distanceSquare(&elem[i * 3], &elem[(i + 1) % size * 3]);
    }
    double max = l[0];
    for (int i = 1; i < size; i++) {
        if (l[i] > max)
            max = l[i];
    }
    return max;
}

/***********************************************************************************************
 * \brief Returns the intersection between a line and a triangle
 * \param[out] Flag whether the intersection has been found
 * \param[in] _coordsTriangle Coordinates of the triangle
 * \param[in] _coordsNode Coordinates of the point
 * \param[in] _direction Direction of the line
 * \param[in/out] _localCoords Local coordinates of the point
 * \author Chenshen Wu
 ***********/
bool computeIntersectionBetweenLineAndTriangle(const double *_X, const double* _X0,
        const double* _n, double* _localCoords) {
    double A[9];
    //  A(1:3,1) = X1-X3;
    for (int i = 0; i < 3; i++)
        A[i * 3] = _X[i] - _X[i + 6];

    //  A(1:3,2) = X2-X3;
    for (int i = 0; i < 3; i++)
        A[i * 3 + 1] = _X[i + 3] - _X[i + 6];

    //  A(1:3,3) = -n;
    for (int i = 0; i < 3; i++)
        A[i * 3 + 2] = -_n[i];

    double b[3];
    // b = X0 - X3
    b[0] = _X0[0] - _X[6];
    b[1] = _X0[1] - _X[7];
    b[2] = _X0[2] - _X[8];

    solve3x3LinearSystem(A, b, EPS_IVERTIBILITYOFSQUAREMATRICES);
    _localCoords[0] = b[0];
    _localCoords[1] = b[1];
    if (fabs(_localCoords[0] - 0.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[0] = 0.0;
    if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[0] = 1.0;
    if (fabs(_localCoords[0] - 0.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[1] = 0.0;
    if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[1] = 1.0;

    if (_localCoords[0] >= 0.0 && _localCoords[0] <= 1.0 && _localCoords[1] >= 0.0
            && _localCoords[1] <= 1.0)
        return true;
    else
        return false;

}

/***********************************************************************************************
 * \brief Returns the intersection between a line and a quadrilateral
 * \param[out] Flag whether the intersection has been found
 * \param[in] _coordsQuad Coordinates of the quadrilateral
 * \param[in] _coordsNode Coordinates of the point
 * \param[in] _direction Direction of the line
 * \param[in/out] _localCoords Local coordinates of the point
 * \author Chenshen Wu
 ***********/
bool computeIntersectionBetweenLineAndQuad(const double *_X, const double* _X0, const double* _n,
        double* _localCoords) {
    double x[3] = { 0.0, 0.0, 0.0 };
    double f[3];
    double df[9];
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 3; j++) {
            // f = 0.25*((-X1+X2+X3-X4)*x(1)+(-X1-X2+X3+X4)*x(2)+(X1-X2+X3-X4)*x(1)*x(2)+X1+X2+X3+X4)-x(3)*n-X0;
            f[j] = 0.25
                    * ((-_X[j] + _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[0]
                            + (-_X[j] - _X[j + 3] + _X[j + 6] + _X[j + 9]) * x[1]
                            + (_X[j] - _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[0] * x[1] + _X[j]
                            + _X[j + 3] + _X[j + 6] + _X[j + 9]) - x[2] * _n[j] - _X0[j];
            //df = [.25*((-X1+X2+X3-X4)+(X1-X2+X3-X4)*x(2)) 0.25*((-X1-X2+X3+X4)+(X1-X2+X3-X4)*x(1)) -n];
            df[j * 3] = 0.25
                    * ((-_X[j] + _X[j + 3] + _X[j + 6] - _X[j + 9])
                            + (_X[j] - _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[1]);
            df[j * 3 + 1] = 0.25
                    * ((-_X[j] - _X[j + 3] + _X[j + 6] + _X[j + 9])
                            + (_X[j] - _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[0]);
            df[j * 3 + 2] = -_n[j];
        }
        if (sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]) > 1e-12) {
            solve3x3LinearSystem(df, f, EPS_IVERTIBILITYOFSQUAREMATRICES);

            for (int j = 0; j < 3; j++)
                x[j] -= f[j];
        } else {
            _localCoords[0] = x[0];
            _localCoords[1] = x[1];
            // Clamp w
            if (fabs(_localCoords[0] + 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[0] = -1.0;
            if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[0] = 1.0;
            // Clamp z
            if (fabs(_localCoords[1] + 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[1] = -1.0;
            if (fabs(_localCoords[1] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[1] = 1.0;
            // Return validity criteria
            if (_localCoords[0] >= -1.0 && _localCoords[0] <= 1.0 && _localCoords[1] >= -1.0
                    && _localCoords[1] <= 1.0)
                return true;
            else
                return false;
        }
    }
    return false;

}

/***********************************************************************************************
 * \brief Computes the area of a given triangle defined by two vectors
 * \param[in] _x1, _y1, _z1 Vector 1
 * \param[in] _x2, _y2, _z2 Vector 2
 * \param[out] The area of the triangle included by the two vectors
 * \author Chenshen Wu
 ***********/
double computeAreaTriangle(double _x1, double _y1, double _z1, double _x2, double _y2, double _z2) {
    double x = _y1 * _z2 - _y2 * _z1;
    double y = _z1 * _x2 - _z2 * _x1;
    double z = _x1 * _y2 - _x2 * _y1;
    return sqrt(x * x + y * y + z * z) / 2;

}

/***********************************************************************************************
 * \brief Computes the distance between two points
 * \param[in] _x1 Point 1
 * \param[in] _x2 Point 2
 * \param[out] The distance between the two given points
 * \author Chenshen Wu
 ***********/
double computePointDistance(double* _x1, double* _x2) {
    return sqrt(
            (_x1[0] - _x2[0]) * (_x1[0] - _x2[0]) + (_x1[1] - _x2[1]) * (_x1[1] - _x2[1])
                    + (_x1[2] - _x2[2]) * (_x1[2] - _x2[2]));
}

/***********************************************************************************************
 * \brief Clean the input polygon, remove colinear or identical point
 * \param[in/out] _polygon, vector of points coordinates
 * \author Fabien Pean
 ***********/
void cleanPolygon(std::vector<double>& polygon) {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polygon.size()/2;
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polygon[p1*2]==polygon[p2*2] && polygon[p1*2+1]==polygon[p2*2+1];

		double v1x=polygon[p2*2]-polygon[p1*2];
		double v1y=polygon[p2*2+1]-polygon[p1*2+1];
		double v2x=polygon[p3*2]-polygon[p1*2];
		double v2y=polygon[p3*2+1]-polygon[p1*2+1];
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		// Edit Aditya Ghantasala
		double cProd[3];
		//double v=computeCrossProduct2D(v1x,v1y,v2x,v2y);
		crossProduct(cProd, v1, v2);
		double v= EMPIRE::MathLibrary::vector2norm(cProd,3);

		// Result of cross product only in Z direction
		bool isColinear=fabs(v)<1e-9?true:false;
		//double n1=v1x*v1x+v1y*v1y;
		//double n2=v2x*v2x+v2y*v2y;
		//bool isWrong=n1>3*fmin(n1,n2);
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polygon.erase(polygon.begin()+p2*2+1);
			polygon.erase(polygon.begin()+p2*2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polygon.size()/2;
			k--;
		}
	}
}

/***********************************************************************************************
 * \brief Clean the input polygon, remove colinear or identical point
 * \param[in/out] _polygon, vector of points coordinates
 * \author Fabien Pean
 ***********/
void cleanPolygon(std::vector<pair<double,double> >& polygon) {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polygon.size();
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polygon[p1]==polygon[p2];

		double v1x=polygon[p2].first-polygon[p1].first;
		double v1y=polygon[p2].second-polygon[p1].second;
		double v2x=polygon[p3].first-polygon[p1].first;
		double v2y=polygon[p3].second-polygon[p1].second;
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		//double v=computeCrossProduct2D(v1x,v1y,v2x,v2y);
		// Edit Aditya
		double cProd[3];
		EMPIRE::MathLibrary::crossProduct(cProd, v1, v2);
		double v = EMPIRE::MathLibrary::vector2norm(cProd,3);
		// Result of cross product only in Z direction
		bool isColinear=fabs(v)<1e-9?true:false;
		//double n1=v1x*v1x+v1y*v1y;
		//double n2=v2x*v2x+v2y*v2y;
		//bool isWrong=n1>3*fmin(n1,n2);
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polygon.erase(polygon.begin()+p2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polygon.size();
			k--;
		}
	}
}

/***********************************************************************************************
 * \brief Clean the input polygon, remove colinear or identical point
 * \param[in/out] _polygon, vector of points coordinates
 * \author Fabien Pean
 ***********/
void cleanPolygon(std::vector<pair<double,double> >& polygon,std::vector<pair<double,double> >& polygonSlave) {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polygon.size();
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polygon[p1]==polygon[p2];

		double v1x=polygon[p2].first-polygon[p1].first;
		double v1y=polygon[p2].second-polygon[p1].second;
		double v2x=polygon[p3].first-polygon[p1].first;
		double v2y=polygon[p3].second-polygon[p1].second;
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		//double v=computeCrossProduct2D(v1x,v1y,v2x,v2y);
		// Edit Aditya
		double cProd[3];
		EMPIRE::MathLibrary::crossProduct(cProd,v1, v2);
		double v = EMPIRE::MathLibrary::vector2norm(cProd,3);
		// Result of cross product only in Z direction
		bool isColinear=fabs(v)<1e-9?true:false;
		//double n1=v1x*v1x+v1y*v1y;
		//double n2=v2x*v2x+v2y*v2y;
		//bool isWrong=n1>3*fmin(n1,n2);
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polygon.erase(polygon.begin()+p2);
			polygonSlave.erase(polygonSlave.begin()+p2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polygon.size();
			k--;
		}
	}
}

bool findIfPointIsInside2DPolygon(int _numVertices, std::vector<double> _polygon, double* _point) {
    /*
     * Returns a flag on whether the given point lies inside the given polygon in 2D.
     *
     *    Input :
     *   _point : double vector with 2 elements
     * _polygon : standard vector of doubles containing the x and y components of each point comprising the polygon sequentially. The points must be ordered coherently.
     */

    int i, j, c = 0;
    for (i = 0, j = _numVertices-1; i < _numVertices; j = i++) {
        if ( ((_polygon[2*i + 1]>_point[1]) != (_polygon[2*j + 1]>_point[1])) &&
            (_point[0] < (_polygon[2*j]-_polygon[2*i]) * (_point[1]-_polygon[2*i + 1]) / (_polygon[2*j + 1]-_polygon[2*i + 1]) + _polygon[2*i]) )
            c = !c;
    }
    return c;
}

// Class Methods
// %%%%%%%%%%%%%%%%%%%%%%
const double PolygonClipper::EPS1 = 1e-10;
const double PolygonClipper::EPS2 = 1e-10;
/***********************************************************************************************
 * \brief Constructor
 * \param[in] _polygonWindow the polygon that clips other polygons
 * \param[in] _sizePolygonWindow number of nodes/edges in this polygon
 * \param[in] _planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \author Tianyang Wang
 ***********/
PolygonClipper::PolygonClipper(const double *_polygonWindow, int _sizePolygonWindow,
        int _planeToProject) :
        planeToProject(_planeToProject) {
    assert(_sizePolygonWindow>2);
    // Michael Andre begin: correcting collapsed nodes so clipping algorithm works properly
    const double eps = 1e-15;
    double lastNoSkip[3];
    int iLastNoSkip;
    bool *skip = new bool[_sizePolygonWindow];
    int numSkip;
    int i, ii;

    sizePolygonWindow = _sizePolygonWindow;
    skip[0] = false;
    iLastNoSkip = 0;
    copyPoint(&_polygonWindow[0], lastNoSkip);

    for (i = 1; i < _sizePolygonWindow; i++) {
        skip[i] = (distanceSquare(lastNoSkip, &_polygonWindow[i * 3]) < eps);
        if (!skip[i]) {
            copyPoint(&_polygonWindow[i * 3], lastNoSkip);
            iLastNoSkip = i;
        }
    }

    if (distanceSquare(&_polygonWindow[0], &_polygonWindow[iLastNoSkip * 3]) < eps)
        skip[iLastNoSkip] = true;

    numSkip = 0;
    for (i = 0; i < _sizePolygonWindow; i++)
        if (skip[i])
            numSkip++;

    sizePolygonWindow = _sizePolygonWindow - numSkip;
    assert(sizePolygonWindow > 2);

    polygonWindow = new double[3 * sizePolygonWindow];

    ii = 0;
    for (i = 0; i < _sizePolygonWindow; i++) {
        if (!skip[i]) {
            copyPoint(&_polygonWindow[i * 3], &polygonWindow[ii * 3]);
            ii++;
        }
    }
    delete[] skip;
    // Michael Andre end

    // numbering the edge and point in the following way:
    // The edge i is the edge between point i and i+1
    insideFlag = new bool[sizePolygonWindow];
    edgeSlopes = new double[sizePolygonWindow];
    reverseXY = new bool[sizePolygonWindow];

    const int XX = (planeToProject + 1) % 3;
    const int YY = (planeToProject + 2) % 3;

    for (int i = 0; i < sizePolygonWindow; i++) {
        // i is the ID of the edge
        double x1, y1, x2, y2, x3, y3;
        int p1_pos, p2_pos, p3_pos; // index of points of the edge in the triangle
        p1_pos = i;
        p2_pos = (i + 1) % sizePolygonWindow;
        p3_pos = (i + 2) % sizePolygonWindow;
        x1 = polygonWindow[p1_pos * 3 + XX];
        y1 = polygonWindow[p1_pos * 3 + YY];
        x2 = polygonWindow[p2_pos * 3 + XX];
        y2 = polygonWindow[p2_pos * 3 + YY];
        x3 = polygonWindow[p3_pos * 3 + XX];
        y3 = polygonWindow[p3_pos * 3 + YY];
        if (fabs(x2 - x1) < fabs(y2 - y1)) {
            reverseXY[i] = true;
            edgeSlopes[i] = (x2 - x1) / (y2 - y1);
            insideFlag[i] = (x3 - x1) > edgeSlopes[i] * (y3 - y1);
        } else {
            reverseXY[i] = false;
            edgeSlopes[i] = (y2 - y1) / (x2 - x1);
            insideFlag[i] = (y3 - y1) > edgeSlopes[i] * (x3 - x1);
        }
    }
}

/***********************************************************************************************
 * \brief Destructor
 * \author Tianyang Wang
 ***********/
PolygonClipper::~PolygonClipper() {
    delete[] insideFlag;
    delete[] edgeSlopes;
    delete[] reverseXY;
    delete[] polygonWindow;
}

/***********************************************************************************************
 * \brief Clip.
 * \param[in] polygonToBeClipped polygon to be clipped by the window polygon
 * \param[in] sizePolygonToBeClipped size of the polygon
 * \param[out] the polygonResult results from the clipping
 * \return if there is a clip, return true, otherwise, return false
 * \author Tianyang Wang
 ***********/
bool PolygonClipper::clip(const double *polygonToBeClipped, int sizePolygonToBeClipped,
        std::vector<double*> *polygonResult) {
    // numbering the edge and point in the following way:
    // The edge i is the edge between point i and i+1.
    // The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice, 2nd edition, Addison-Wesley" P237 - P240.
    list<double*> *listInput = new list<double*>;
    list<double*> *listOutput = new list<double*>;

    // initialize listOutput to be the slave nodes
    for (int i = 0; i < sizePolygonToBeClipped; i++) {
        double *point = new double[3];
        copyPoint(&polygonToBeClipped[i * 3], point);
        listOutput->push_back(point);
    }
    //printElem(polygonWindow, sizePolygonWindow);
    //printElem(polygonToBeClipped, sizePolygonToBeClipped);
    const double TOL_SQR = EPS1 * EPS1 * longestEdgeLengthSquare(polygonWindow, sizePolygonWindow);

    for (int i = 0; i < sizePolygonWindow; i++) {
        // i is the edge ID
        // 1. make the output the new input
        for (list<double*>::iterator it = listInput->begin(); it != listInput->end(); it++)
            delete[] *it;
        listInput->clear();
        for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) {
            double *point = new double[3];
            copyPoint(*it, point);
            listInput->push_back(point);
        }
        for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++)
            delete[] *it;
        listOutput->clear();

        const double *edgeP1 = &polygonWindow[i * 3];
        const double *edgeP2 = &polygonWindow[(i + 1) % sizePolygonWindow * 3];
        int size = listInput->size();
        // 2. clip the input list by edge, create new output list
        for (list<double*>::iterator it = listInput->begin(); it != listInput->end(); it++) {
            double *p1 = *it;
            list<double*>::iterator it2 = it;
            advance(it2, 1);
            if (it2 == listInput->end())
                it2 = listInput->begin();
            double *p2 = *it2;
            bool inside1 = inside(i, p1);
            bool inside2 = inside(i, p2);
            if (!inside1 && inside2) { // there should be intersection, unless two lines are overlapped
                double *point1 = new double[3];
                if (intersect(edgeP1, edgeP2, p1, p2, planeToProject, point1))
                    listOutput->push_back(point1);
                double *point2 = new double[3];
                copyPoint(p2, point2);
                listOutput->push_back(point2);
            } else if (inside1 && inside2) {
                double *point = new double[3];
                copyPoint(p2, point);
                listOutput->push_back(point);
            } else if (inside1 && !inside2) { // there should be intersection, unless two lines are overlapped
                double *point = new double[3];
                if (intersect(edgeP1, edgeP2, p1, p2, planeToProject, point))
                    listOutput->push_back(point);
            } else {
                // do nothing
            }
        }

        /*cout << "===================================" << i << endl;
         printPoint(edgeP1);
         printPoint(edgeP2);
         double tmp[listOutput->size() * 3];
         int count = 0;
         for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) {
         for (int jj = 0; jj < 3; jj++)
         tmp[count * 3 + jj] = (*it)[jj];
         count++;
         }
         printElem(tmp, listOutput->size());
         cout << "===================================" << endl;*/

        for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) { // remove overlapped points in listOutput
            if (listOutput->size() == 1)
                break; // otherwise, segmentation fault
            double *p1 = *it;
            list<double*>::iterator it2 = it;
            advance(it2, 1);
            if (it2 == listOutput->end())
                it2 = listOutput->begin();
            double *p2 = *it2;
            if (distanceSquare(p1, p2) < TOL_SQR) {
                delete[] *it2; // otherwise, memory leak
                listOutput->erase(it2); // listOutput->erase(it) will cause segmentation fault
            }
        }
    }

    for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) {
        double *point = new double[3];
        copyPoint(*it, point);
        polygonResult->push_back(point);
    }

    // release the storage
    for (list<double*>::iterator it = listInput->begin(); it != listInput->end(); it++)
        delete[] *it;
    delete listInput;

    for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++)
        delete[] *it;
    delete listOutput;

    // judge whether there is really overlapped area or not
    int size = polygonResult->size();

    if (size < 3)
        return false;

    double *tmp = new double[size * 3];
    for (int i = 0; i < size; i++)
        for (int j = 0; j < 3; j++)
            tmp[i * 3 + j] = polygonResult->at(i)[j];
    if (computePolygonArea(tmp, size) < TOL_SQR) {
        delete[] tmp;
        return false; // it could happen the points are almost on the same line
    }
    delete[] tmp;

    return true;
    // since there is error when computing the intersection, now we have to check whether the center
    // of the points is inside both elements
    /*double polygon[polygonResult->size() * 3];
     for (int i = 0; i < polygonResult->size(); i++)
     for (int j = 0; j < 3; j++)
     polygon[i * 3 + j] = polygonResult->at(i)[j];
     double polygonCenter[3];
     computePolygonCenter(polygon, polygonResult->size(), polygonCenter);
     if (sizePolygonWindow == 3) {
     double localCoor[3];
     if (!computeLocalCoorInTriangle(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else if (sizePolygonWindow == 4) {
     double localCoor[2];
     if (!computeLocalCoorInQuad(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else {
     assert(false);
     }
     if (sizePolygonWindow == 3) {
     double localCoor[3];
     if (!computeLocalCoorInTriangle(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else if (sizePolygonWindow == 4) {
     double localCoor[2];
     if (!computeLocalCoorInQuad(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else {
     assert(false);
     }*/
}

/***********************************************************************************************
 * \brief Decides whether a point is on the "inside" side of the i-th edge of the polygon
 * \param[in] edgeID id of the edge
 * \param[in] point the point
 * \return if inside, return true, otherwise, return false
 * \author Tianyang Wang
 ***********/
bool PolygonClipper::inside(int edgeID, double *point) {
    const int XX = (planeToProject + 1) % 3;
    const int YY = (planeToProject + 2) % 3;
    int p1_pos = edgeID;
    double x1 = polygonWindow[p1_pos * 3 + XX];
    double y1 = polygonWindow[p1_pos * 3 + YY];
    double x4 = point[XX];
    double y4 = point[YY];
    if (reverseXY[edgeID])
        return ((x4 - x1) > edgeSlopes[edgeID] * (y4 - y1)) == insideFlag[edgeID];
    else
        return ((y4 - y1) > edgeSlopes[edgeID] * (x4 - x1)) == insideFlag[edgeID];
}

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
bool PolygonClipper::intersect(const double *la0, const double *la1, const double *lb0,
        const double *lb1, int _planeToProject, double *intersection) {
    // equation (in vector): alpha a0a1 - beta b0b1 = a0b0
    double unknown[2];
    double a0a1[3];
    double b0b1[3];
    double a0b0[3];

    for (int i = 0; i < 3; i++) {
        a0a1[i] = la1[i] - la0[i];
        b0b1[i] = lb1[i] - lb0[i];
        a0b0[i] = lb0[i] - la0[i];
    }

    double A[4];
    double a0a1_2D[2];
    double b0b1_2D[2];

    int XX = (_planeToProject + 1) % 3;
    int YY = (_planeToProject + 2) % 3;

    A[0] = a0a1_2D[0] = a0a1[XX];
    A[2] = b0b1_2D[0] = -b0b1[XX];
    unknown[0] = a0b0[XX];
    A[1] = a0a1_2D[1] = a0a1[YY];
    A[3] = b0b1_2D[1] = -b0b1[YY];
    unknown[1] = a0b0[YY];

    // if the system is not solvable, it means la and lb are parallel
    if (!solve2x2LinearSystem(A, unknown, EPS2))
        return false;

    /* intersection point should be on lb (0=<unknown[1]<=1), except:
     * 1. if la and lb are 'very' parallel, the numerical error could be large, therefore
     *    the intersection point may be outside of lb.
     *   (if la and lb are 'very' parallel, we can return any point on lb, or do not return
     *    any point, which won't cause problems)
     * 2. if end point of lb locates exactly on la, due to numerical error, the intersection
     *    point may locate 'a little bit' out of lb.
     * For both cases, it is enough to say, if unknown[1] is out of range, pick up one end
     * point of lb as the intersection point.
     */

    if (unknown[1] < 0.0)
        unknown[1] = 0.0;
    else if (unknown[1] > 1.0)
        unknown[1] = 1.0;

    for (int i = 0; i < 3; i++)
        intersection[i] = lb0[i] + unknown[1] * b0b1[i];
    /*for (int i = 0; i < 3; i++) {
     double tmp = intersection[i];
     assert(tmp==intersection[i]);//check whether intersection[i] is NAN
     }*/

    return true;
}

/***********************************************************************************************
 * \brief Constructor
 * \param[in] _u1, _u2, _v1, _v2 define the IGA element that clips other the fluid element.
 * \author Chenshen Wu
 ***********/
IGAPolygonClipper::IGAPolygonClipper(const double _u1, const double _u2, const double _v1,
        const double _v2) {

    polygonWindow[0] = _u1;
    polygonWindow[1] = _v1;
    polygonWindow[2] = 0;
    polygonWindow[3] = _u2;
    polygonWindow[4] = _v1;
    polygonWindow[5] = 0;
    polygonWindow[6] = _u2;
    polygonWindow[7] = _v2;
    polygonWindow[8] = 0;
    polygonWindow[9] = _u1;
    polygonWindow[10] = _v2;
    polygonWindow[11] = 0;
    clipper = new EMPIRE::MathLibrary::PolygonClipper(polygonWindow, 4, 2);
}

/***********************************************************************************************
 * \brief Destructor
 * \author Chenshen Wu
 ***********/
IGAPolygonClipper::~IGAPolygonClipper() {
    delete clipper;
}

/***********************************************************************************************
 * \brief Clip.
 * \param[in] polygonToBeClipped polygon to be clipped by the window polygon
 * \param[in] sizePolygonToBeClipped size of the polygon
 * \param[out] the polygonResult results from the clipping
 * \return if there is a clip, return true, otherwise, return false
 * \author Tianyang Wang
 ***********/
bool IGAPolygonClipper::clip(const double *_polygonToBeClipped, int _sizePolygonToBeClipped,
        vector<double*> *_polygonResult) {
    double polygonToBeClipped[_sizePolygonToBeClipped * 3];
    for (int i = 0; i < _sizePolygonToBeClipped; i++) {
        polygonToBeClipped[i * 3] = _polygonToBeClipped[i * 2];
        polygonToBeClipped[i * 3 + 1] = _polygonToBeClipped[i * 2 + 1];
        polygonToBeClipped[i * 3 + 2] = 0;
    }
    return clipper->clip(polygonToBeClipped, _sizePolygonToBeClipped, _polygonResult);
}



}
}

