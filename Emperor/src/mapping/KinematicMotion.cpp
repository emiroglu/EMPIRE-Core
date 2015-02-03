#include "KinematicMotion.h"
#include <math.h>
#include <assert.h>
#include <iostream>

using namespace std;

namespace EMPIRE {

KinematicMotion::KinematicMotion() {
    translationVector = new double[3];
    resetTranslation();
    rotationMatrix = new double[9];
    resetRotation();
}

KinematicMotion::~KinematicMotion() {
    delete[] translationVector;
    delete[] rotationMatrix;
}

void KinematicMotion::resetTranslation() {
    translationVector[0] = 0.0;
    translationVector[1] = 0.0;
    translationVector[2] = 0.0;
}

void KinematicMotion::resetRotation() {
    for (int i = 0; i < 9; i++)
        rotationMatrix[i] = 0.0;
    rotationMatrix[0] = 1.0;
    rotationMatrix[1 * 3 + 1] = 1.0;
    rotationMatrix[2 * 3 + 2] = 1.0;
}

const double *KinematicMotion::getTranslationVector() const {
    return translationVector;
}

const double *KinematicMotion::getRotationMatrix() const {
    return rotationMatrix;
}

void KinematicMotion::getRotationVector(double *rot) const {
    // see Non-linear Modeling and Analysis of Solids and Structures (Steen Krenk 2009) P52
    double angle = rotationMatrix[0] + rotationMatrix[4] + rotationMatrix[8] - 1.0;
    angle /= 2.0;
    angle = acos(angle); // between 0 and pi

    const double EPS = 1E-10;
    if (angle < EPS) {
        rot[0] = 0.0;
        rot[1] = 0.0;
        rot[2] = 0.0;
        return;
    } else if ((M_PI - angle) < EPS) {
        const double product11 = (rotationMatrix[0 * 3 + 0] + 1.0) / 2.0;
        const double product22 = (rotationMatrix[1 * 3 + 1] + 1.0) / 2.0;
        const double product33 = (rotationMatrix[2 * 3 + 2] + 1.0) / 2.0;
        const double product12 = (rotationMatrix[0 * 3 + 1] + 1.0) / 2.0;
        const double product23 = (rotationMatrix[1 * 3 + 2] + 1.0) / 2.0;
        const double product13 = (rotationMatrix[0 * 3 + 2] + 1.0) / 2.0;
        const double tmp1 = sqrt(product11);
        const double tmp2 = sqrt(product22);
        const double tmp3 = sqrt(product33);

        { // case 1 +++:
            rot[0] = tmp1;
            rot[1] = tmp2;
            rot[2] = tmp3;
            const double tmp12 = rot[0] * rot[1];
            const double tmp13 = rot[0] * rot[2];
            const double tmp23 = rot[1] * rot[2];
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rot[0] *= M_PI;
                        rot[1] *= M_PI;
                        rot[2] *= M_PI;
                        return;
                    }
        }
        { // case 2 +--:
            rot[0] = tmp1;
            rot[1] = -tmp2;
            rot[2] = -tmp3;
            const double tmp12 = rot[0] * rot[1];
            const double tmp13 = rot[0] * rot[2];
            const double tmp23 = rot[1] * rot[2];
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rot[0] *= M_PI;
                        rot[1] *= M_PI;
                        rot[2] *= M_PI;
                        return;
                    }
        }
        { // case 3 -+-:
            rot[0] = -tmp1;
            rot[1] = tmp2;
            rot[2] = -tmp3;
            const double tmp12 = rot[0] * rot[1];
            const double tmp13 = rot[0] * rot[2];
            const double tmp23 = rot[1] * rot[2];
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rot[0] *= M_PI;
                        rot[1] *= M_PI;
                        rot[2] *= M_PI;
                        return;
                    }
        }
        { // case 4 --+:
            rot[0] = -tmp1;
            rot[1] = -tmp2;
            rot[2] = tmp3;
            const double tmp12 = rot[0] * rot[1];
            const double tmp13 = rot[0] * rot[2];
            const double tmp23 = rot[1] * rot[2];
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rot[0] *= M_PI;
                        rot[1] *= M_PI;
                        rot[2] *= M_PI;
                        return;
                    }
        }
        assert(false);
    }

    double tmp = angle / 2.0 / sin(angle);
    rot[0] = -(rotationMatrix[1 * 3 + 2] - rotationMatrix[2 * 3 + 1]) * tmp;
    rot[1] = (rotationMatrix[0 * 3 + 2] - rotationMatrix[2 * 3 + 0]) * tmp;
    rot[2] = -(rotationMatrix[0 * 3 + 1] - rotationMatrix[1 * 3 + 0]) * tmp;
}

void KinematicMotion::addKinematicMotion(const KinematicMotion *kinematicMotion) {
    addRotation(kinematicMotion->getRotationMatrix()); // step 1
    addTranslation(kinematicMotion->getTranslationVector()); // step 2, order must not be changed
}

KinematicMotion *KinematicMotion::newInverse() {
    // x2 = R*x1 + T  ==>   x1 = R^{-1} * (x2-T) = R^{-1}*x2 - R^{-1}*T
    double rotationMatrixInverse[9];
    copyMatrix(rotationMatrix, rotationMatrixInverse);
    matrixTanspose(rotationMatrixInverse);
    double translationVectorInverse[3];
    copyVector(translationVector, translationVectorInverse);
    vectorInverse(translationVectorInverse);
    KinematicMotion *t_inverse = new KinematicMotion();
    t_inverse->addTranslation(translationVectorInverse); // step 1
    t_inverse->addRotation(rotationMatrixInverse); // step 2, order cannot be changed
    return t_inverse;
}

void KinematicMotion::addTranslation(const double *_translationVector) {
    // translation is defined in the CURRENT coordinates system
    translationVector[0] += _translationVector[0];
    translationVector[1] += _translationVector[1];
    translationVector[2] += _translationVector[2];
}

void KinematicMotion::addRotation(const double *axis, bool normalized, double angle) {
    // rotation is defined in the CURRENT coordinates system
    double u[3];
    copyVector(axis, u);
    if (!normalized) {
        normalizeVector(u);
    }
    // algorithm from http://en.wikipedia.org/wiki/Rotation_matrix
    double c = cos(angle);
    double s = sin(angle);

    double newRotation[9];

    newRotation[0] = c + u[0] * u[0] * (1.0 - c);
    newRotation[1] = u[0] * u[1] * (1.0 - c) - u[2] * s;
    newRotation[2] = u[0] * u[2] * (1.0 - c) + u[1] * s;

    newRotation[3] = u[1] * u[0] * (1.0 - c) + u[2] * s;
    newRotation[4] = c + u[1] * u[1] * (1.0 - c);
    newRotation[5] = u[1] * u[2] * (1.0 - c) - u[0] * s;

    newRotation[6] = u[2] * u[0] * (1.0 - c) - u[1] * s;
    newRotation[7] = u[2] * u[1] * (1.0 - c) + u[0] * s;
    newRotation[8] = c + u[2] * u[2] * (1.0 - c);

    //R'*(R*x + t) = R'*R*x +  R'*t
    matrixProduct(newRotation, rotationMatrix);
    matrixVectorProduct(newRotation, translationVector);
}

void KinematicMotion::addRotation(const double *vec0, const double *vec1, bool normalized) { // not unique, but shortest path
    // rotation is defined in the CURRENT coordinates system
    double e0[3];
    copyVector(vec0, e0);
    double e1[3];
    copyVector(vec1, e1);
    if (!normalized) {
        normalizeVector(e0);
        normalizeVector(e1);
    }
    // algorithm form Non-linear Modeling and Analysis of Solids and Structures (Steen Krenk) P54
    double e2[3];
    e2[0] = e0[0] + e1[0];
    e2[1] = e0[1] + e1[1];
    e2[2] = e0[2] + e1[2];
    double e2_square = vectorProduct(e2, e2);

    double newRotation[9];

    newRotation[0] = 1.0 + 2.0 * e1[0] * e0[0] - 2.0 / e2_square * e2[0] * e2[0];
    newRotation[1] = 0.0 + 2.0 * e1[0] * e0[1] - 2.0 / e2_square * e2[0] * e2[1];
    newRotation[2] = 0.0 + 2.0 * e1[0] * e0[2] - 2.0 / e2_square * e2[0] * e2[2];

    newRotation[3] = 0.0 + 2.0 * e1[1] * e0[0] - 2.0 / e2_square * e2[1] * e2[0];
    newRotation[4] = 1.0 + 2.0 * e1[1] * e0[1] - 2.0 / e2_square * e2[1] * e2[1];
    newRotation[5] = 0.0 + 2.0 * e1[1] * e0[2] - 2.0 / e2_square * e2[1] * e2[2];

    newRotation[6] = 0.0 + 2.0 * e1[2] * e0[0] - 2.0 / e2_square * e2[2] * e2[0];
    newRotation[7] = 0.0 + 2.0 * e1[2] * e0[1] - 2.0 / e2_square * e2[2] * e2[1];
    newRotation[8] = 1.0 + 2.0 * e1[2] * e0[2] - 2.0 / e2_square * e2[2] * e2[2];

    //R'*(R*x + t) = R'*R*x +  R'*t
    matrixProduct(newRotation, rotationMatrix);
    matrixVectorProduct(newRotation, translationVector);
}

void KinematicMotion::addRotation(const double *xAxisNew, const double *yAxisNew,
        const double *zAxisNew, bool normalized) {
    double e0[3];
    double e1[3];
    double e2[3];
    copyVector(xAxisNew, e0);
    copyVector(yAxisNew, e1);
    copyVector(zAxisNew, e2);
    if (!normalized) {
        normalizeVector(e0);
        normalizeVector(e1);
        normalizeVector(e2);
    }

    double newRotation[9];
    for (int i = 0; i < 3; i++) {
        newRotation[i * 3 + 0] = e0[i];
        newRotation[i * 3 + 1] = e1[i];
        newRotation[i * 3 + 2] = e2[i];
    }

    //R'*(R*x + t) = R'*R*x +  R'*t
    matrixProduct(newRotation, rotationMatrix);
    matrixVectorProduct(newRotation, translationVector);
}

void KinematicMotion::addRotation(const double *_rotationMatrix) {
    //R'*(R*x + t) = R'*R*x +  R'*t
    matrixProduct(_rotationMatrix, rotationMatrix);
    matrixVectorProduct(_rotationMatrix, translationVector);
}

void KinematicMotion::move(double *coordinates) const {
    // rotate first, then translate: x=R*x + T. (Remark: R*(x+T) is wrong)
    matrixVectorProduct(rotationMatrix, coordinates);
    vectorAddition(coordinates, translationVector);
}

void KinematicMotion::checkRotationCorrectness() const {
    double M[9];
    double M_inv[9];
    copyMatrix(rotationMatrix, M);
    copyMatrix(rotationMatrix, M_inv);
    matrixTanspose(M_inv);
    matrixProduct(M_inv, M); // M = M_inv * M
    const double TOL = 1E-10;
    //printMatrix(M);
    assert(fabs(M[0 * 3 + 0] - 1.0) < TOL);
    assert(fabs(M[0 * 3 + 1] - 0.0) < TOL);
    assert(fabs(M[0 * 3 + 2] - 0.0) < TOL);
    assert(fabs(M[1 * 3 + 0] - 0.0) < TOL);
    assert(fabs(M[1 * 3 + 1] - 1.0) < TOL);
    assert(fabs(M[1 * 3 + 2] - 0.0) < TOL);
    assert(fabs(M[2 * 3 + 0] - 0.0) < TOL);
    assert(fabs(M[2 * 3 + 1] - 0.0) < TOL);
    assert(fabs(M[2 * 3 + 2] - 1.0) < TOL);
}

void KinematicMotion::print() const {
    cout << "***Transfromation***: " << endl;
    printMatrix(rotationMatrix);
    printVector(translationVector);
    cout << "********************* " << endl << endl;
}

void KinematicMotion::normalizeVector(double *vector) {
    double length = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    assert(length > 1E-20); // 1E-20 is a physical "small" length
    vector[0] /= length;
    vector[1] /= length;
    vector[2] /= length;
}

void KinematicMotion::matrixProduct(const double *A, double *B) { // row-major
    double C[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i * 3 + j] = 0.0;
            for (int k = 0; k < 3; k++) {
                C[i * 3 + j] += A[i * 3 + k] * B[k * 3 + j];
            }
        }
    }
    for (int i = 0; i < 9; i++)
        B[i] = C[i];
}

void KinematicMotion::copyMatrix(const double *M, double *MCopy) {
    for (int i = 0; i < 9; i++)
        MCopy[i] = M[i];
}

void KinematicMotion::copyVector(const double *v, double *vCopy) {
    for (int i = 0; i < 3; i++)
        vCopy[i] = v[i];
}

double KinematicMotion::vectorProduct(const double *v1, const double *v2) {
    double result = 0.0;
    for (int i = 0; i < 3; i++)
        result += v1[i] * v2[i];
    return result;
}

void KinematicMotion::matrixVectorProduct(const double *M, double *v) {
    double result[] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i] += M[i * 3 + j] * v[j];
        }
    }
    for (int i = 0; i < 3; i++) {
        v[i] = result[i];
    }
}

void KinematicMotion::vectorAddition(double *v1, const double *v2) {
    for (int i = 0; i < 3; i++)
        v1[i] += v2[i];
}

void KinematicMotion::matrixTanspose(double *M) {
    double MCopy[9];
    copyMatrix(M, MCopy);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            MCopy[i * 3 + j] = M[j * 3 + i];
    copyMatrix(MCopy, M);
}

void KinematicMotion::vectorInverse(double *v) {
    for (int i = 0; i < 3; i++)
        v[i] = -v[i];
}

void KinematicMotion::printMatrix(const double *M) {
    cout << "----------Matrix----------" << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << '\t' << M[i * 3 + j];
        }
        cout << endl;
    }
    cout << "--------------------------" << endl;
}

void KinematicMotion::printVector(const double *v) {
    cout << "----------Vector----------" << endl;
    for (int i = 0; i < 3; i++) {
        cout << '\t' << v[i];
    }
    cout << endl;
    cout << "--------------------------" << endl;
}

} /* namespace EMPIRE */
