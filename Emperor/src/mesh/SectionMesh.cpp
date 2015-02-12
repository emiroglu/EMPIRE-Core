#include "SectionMesh.h"
#include "KinematicMotion.h"
#include "Message.h"

namespace EMPIRE {

using namespace std;

SectionMesh::SectionMesh(std::string _name, int _numNodes, int _numElems, bool _triangulateAll) :
        FEMesh(_name, _numNodes, _numElems, _triangulateAll) {
    type = EMPIRE_Mesh_SectionMesh;
    rotationGlobal2Root = new double[9];
    translationGlobal2Root = new double[3];
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

Message &operator<<(Message &message, SectionMesh &mesh) {
    // First validate the mesh
    mesh.validateMesh();
    // Display the mesh information
    message << "\t+" << "FEMesh name: " << mesh.name << endl;
    message << "\t\t+" << "no. of nodes: " << mesh.numNodes << endl;
    message << "\t\t+" << "no. of elements: " << mesh.numElems << endl;
    message << "\t\t+" << "Nodes:" << endl;
    for (int i = 0; i < mesh.numNodes; i++) {
        message << "\t\t+" << "\t" << mesh.nodeIDs[i];
        for (int j = 0; j < 3; j++) {
            message << "\t" << mesh.nodes[i * 3 + j];
        }
        message << endl;
    }
    message << "\t\t+" << "Elements:" << endl;
    int count = 0;
    for (int i = 0; i < mesh.numElems; i++) {
        message << "\t\t+" << '\t' << i + 1;
        for (int j = 0; j < mesh.numNodesPerElem[i]; j++) {
            message << "\t" << mesh.elems[count + j];
        }
        count += mesh.numNodesPerElem[i];
        message << endl;
    }
    message << "\t\t+" << "number of sections: " << mesh.getNumSections() << endl;
    message << "\t\t+" << "number of root section nodes: " << mesh.getNumRootSectionNodes() << endl;
    message << "\t\t+" << "number of normal section nodes: " << mesh.getNumNormalSectionNodes() << endl;
    message << "\t\t+" << "number of tip section nodes: " << mesh.getNumTipSectionNodes() << endl;
    cout << "\t+" << "----------Rotation matrix----------" << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << '\t' << mesh.getRotationGlobal2Root()[i * 3 + j];
        }
        cout << endl;
    }
    cout << "\t+" << "-----------------------------------" << endl;
    cout << "\t+" << "----------Translation vector----------" << endl;
    for (int i = 0; i < 3; i++) {
        cout << '\t' << mesh.getTranslationGlobal2Root()[i];
    }
    cout << endl;
    cout << "\t+" << "--------------------------------------" << endl;
    message() << "\t+" << "---------------------------------" << endl;

    return message;
}

} /* namespace EMPIRE */
