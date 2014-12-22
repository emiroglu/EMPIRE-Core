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
/***********************************************************************************************//**
 * \file KinematicMotion.h
 * This file holds the class KinematicMotion
 * \date 12/22/2014
 **************************************************************************************************/
#ifndef KINEMATICMOTION_H_
#define KINEMATICMOTION_H_

namespace EMPIRE {
/********//**
 * \brief Class KinematicMotion performs rotation and translation of points in 3D: x = R*x + t. It can
 * also perform coordinate transformation between different coordinate systems (see function move()).
 ***********/
class KinematicMotion {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Tianyang Wang
     ***********/
    KinematicMotion();
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~KinematicMotion();
    /***********************************************************************************************
     * \brief Reset the translation vector to {0,0,0}
     * \author Tianyang Wang
     ***********/
    void resetTranslation();
    /***********************************************************************************************
     * \brief Reset the rotation matrix to the identity matrix
     * \author Tianyang Wang
     ***********/
    void resetRotation();
    /***********************************************************************************************
     * \brief Return the translation vector
     * \author Tianyang Wang
     ***********/
    const double *getTranslationVector() const;
    /***********************************************************************************************
     * \brief Return the rotational matrix
     * \author Tianyang Wang
     ***********/
    const double *getRotationMatrix() const;
    /***********************************************************************************************
     * \brief Add another motion to the current motion: R'*(R+t) + t' = R'*R + (R'*t+t')
     * \param[in] _numNodesA number of nodes of A
     * \author Tianyang Wang
     ***********/
    void addKinematicMotion(const KinematicMotion *kinematicMotion);
    /***********************************************************************************************
     * \brief Return the inverse kinematic motion of this one: x' = R*x+t <-> x = R^T*(x'-t) = R^T*x' - R^T*t
     * \param[in] _numNodesA number of nodes of A
     * \author Tianyang Wang
     ***********/
    KinematicMotion *newInverse();
    /***********************************************************************************************
     * \brief Add to the current translation vector: t = t + t'
     * \param[in] _translationVector number of nodes of A
     * \author Tianyang Wang
     ***********/
    void addTranslation(const double *_translationVector);
    /***********************************************************************************************
     * \brief Add rotation to the current motion with given axis and angle: R'*(R*x + t) = R'*R*x +  R'*t
     * \param[in] axis axis of rotation
     * \param[in] normalized whether the axis vector is normalized or not
     * \param[in] angle angle of rotation
     * \author Tianyang Wang
     ***********/
    void addRotation(const double *axis, bool normalized, double angle);
    /***********************************************************************************************
     * \brief Add rotation to the current motion with given two vectors. Rotate from vec1 to vec2, the solution is not unique,
     * so choose the axis orthogonal to the plane of vec1 and vec2: R'*(R*x + t) = R'*R*x +  R'*t
     * \param[in] vec1 vec1
     * \param[in] vec2 vec2
     * \param[in] normalized whether vec1 and vec2 are normalized
     * \author Tianyang Wang
     ***********/
    void addRotation(const double *vec1, const double *vec2, bool normalized);
    /***********************************************************************************************
     * \brief Add rotation to the current motion by rotating to the new coordinate system: R'*(R*x + t) = R'*R*x +  R'*t
     * \param[in] xAxisNew new x axis
     * \param[in] yAxisNew new y axis
     * \param[in] zAxisNew new z axis
     * \param[in] normalized whether the axis vector is normalized or not
     * \author Tianyang Wang
     ***********/
    void addRotation(const double *xAxisNew, const double *yAxisNew, const double *zAxisNew,
            bool normalized);
    /***********************************************************************************************
     * \brief Add rotation to the current motion with given rotation matrix: R'*(R*x + t) = R'*R*x +  R'*t
     * \param[in] _rotationMatrix the rotation matrix to be added
     * \author Tianyang Wang
     ***********/
    void addRotation(const double *_rotationMatrix);
    /***********************************************************************************************
     * \brief Assume this defines the kinematic motion from coordinate system A to B, given a point, if the input is
     * its coordinates in A, then the point is moved according to the motion; if the input is its coordinates in B, then
     * its coordinates in B is obtained.
     * \param[in] coordinates coordinates of a point
     * \author Tianyang Wang
     ***********/
    void move(double *coordinates) const;
    /***********************************************************************************************
     * \brief Print the rotation matrix and the translation vector
     * \author Tianyang Wang
     ***********/
    void print() const;
protected:
    /// the rotation matrix
    double *rotationMatrix;
    /// the translation vector
    double *translationVector;
    /// unit test class
    friend class TestKinematicMotion;

    /***********************************************************************************************
     * \brief Normalize a vector
     * \param[in] vector a vector
     * \author Tianyang Wang
     ***********/
    static void normalizeVector(double *vector);
    /***********************************************************************************************
     * \brief Matrix multiplication: B=A*B
     * \param[in] A marix A
     * \param[in,out] B marix B
     * \author Tianyang Wang
     ***********/
    static void matrixProduct(const double *A, double *B);
    /***********************************************************************************************
     * \brief Copy matrix M to MCopy
     * \param[in] M marix M
     * \param[in,out] MCopy marix MCopy
     * \author Tianyang Wang
     ***********/
    static void copyMatrix(const double *M, double *MCopy);
    /***********************************************************************************************
     * \brief Copy vector v to vCopy
     * \param[in] v vector v
     * \param[in,out] vCopy vector vCopy
     * \author Tianyang Wang
     ***********/
    static void copyVector(const double *v, double *vCopy);
    /***********************************************************************************************
     * \brief Get inner product of two vectors
     * \param[in] v1 vector v1
     * \param[in] v2 vector v2
     * \author Tianyang Wang
     ***********/
    static double vectorProduct(const double *v1, const double *v2);
    /***********************************************************************************************
     * \brief Perform matrix vector multiplication v=M*v
     * \param[in] M matrix M
     * \param[in,out] v vector v
     * \author Tianyang Wang
     ***********/
    static void matrixVectorProduct(const double *M, double *v);
    /***********************************************************************************************
     * \brief Add two vectors (v1 = v1 + v2)
     * \param[in,out] v1 vector v1
     * \param[in] v2 vector v2
     * \author Tianyang Wang
     ***********/
    static void vectorAddition(double *v1, const double *v2);
    /***********************************************************************************************
     * \brief Transpose a matrix
     * \param[in,out] M matrix M
     * \author Tianyang Wang
     ***********/
    static void matrixTanspose(double *M);
    /***********************************************************************************************
     * \brief Inverse a vector (v = -v)
     * \param[in,out] v vector v
     * \author Tianyang Wang
     ***********/
    static void vectorInverse(double *v);
    /***********************************************************************************************
     * \brief Print the matrix
     * \param[in] M matrix M
     * \author Tianyang Wang
     ***********/
    static void printMatrix(const double *M);
    /***********************************************************************************************
     * \brief Print the vector
     * \param[in] v vector v
     * \author Tianyang Wang
     ***********/
    static void printVector(const double *v);
};

} /* namespace EMPIRE */

#endif /* KINEMATICMOTION_H_ */
