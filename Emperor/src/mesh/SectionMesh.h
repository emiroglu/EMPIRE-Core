#ifndef SECTIONMESH_H_
#define SECTIONMESH_H_

#include "FEMesh.h"
#include <string>

namespace EMPIRE {

class SectionMesh : public FEMesh {
public:
    SectionMesh(std::string _name, int _numNodes, int _numElems);
    virtual ~SectionMesh();

    int getNumRootSectionNodes();
    int getNumNormalSectionNodes();
    int getNumTipSectionNodes();

    void setNumRootSectionNodes(int _numRootSectionNodes);
    void setNumNormalSectionNodes(int _numNormalSectionNodes);
    void setNumTipSectionNodes(int _numTipSectionNodes);
    void setTranslationGlobal2Root(double *_translationGlobal2Root);
    void setRotationGlobal2Root(double *_rotationGlobal2Root);

private:
    int numRootSectionNodes;
    int numNormalSectionNodes;
    int numTipSectionNodes;
    double *translationGlobal2Root;
    double *rotationGlobal2Root;

};

} /* namespace EMPIRE */

#endif /* SECTIONMESH_H_ */
