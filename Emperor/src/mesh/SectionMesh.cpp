#include "SectionMesh.h"
#include "KinematicMotion.h"

namespace EMPIRE {

SectionMesh::SectionMesh(std::string _name, int _numNodes, int _numElems) :
        FEMesh(_name, _numNodes, _numElems) {
    type = EMPIRE_Mesh_SectionMesh;
    translationGlobal2Root = new double[3];
    rotationGlobal2Root = new double[9];
}

SectionMesh::~SectionMesh() {
    delete[] translationGlobal2Root;
    delete[] rotationGlobal2Root;
}

int SectionMesh::getNumSections() {
    return numSections;
}

int SectionMesh::getNumRootSectionNodes() {
    return numRootSectionNodes;
}

int SectionMesh::getNumNormalSectionNodes() {
    return numNormalSectionNodes;
}

int SectionMesh::getNumTipSectionNodes() {
    return numTipSectionNodes;
}

const double *SectionMesh::getTranslationGlobal2Root() {
    return translationGlobal2Root;
}

const double *SectionMesh::getRotationGlobal2Root() {
    return rotationGlobal2Root;
}

void SectionMesh::setNumSections(int _numSections) {
    numSections = _numSections;
}

void SectionMesh::setNumRootSectionNodes(int _numRootSectionNodes) {
    numRootSectionNodes = _numRootSectionNodes;
}

void SectionMesh::setNumNormalSectionNodes(int _numNormalSectionNodes) {
    numNormalSectionNodes = _numNormalSectionNodes;
}

void SectionMesh::setNumTipSectionNodes(int _numTipSectionNodes) {
    numTipSectionNodes = _numTipSectionNodes;
}

void SectionMesh::setTranslationGlobal2Root(double *_translationGlobal2Root) {
    for (int i = 0; i < 3; i++)
        translationGlobal2Root[i] = _translationGlobal2Root[i];
}
void SectionMesh::setRotationGlobal2Root(double *_rotationGlobal2Root) {
    for (int i = 0; i < 9; i++)
        rotationGlobal2Root[i] = _rotationGlobal2Root[i];
}

} /* namespace EMPIRE */
