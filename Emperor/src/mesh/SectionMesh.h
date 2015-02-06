#ifndef SECTIONMESH_H_
#define SECTIONMESH_H_

#include "FEMesh.h"
#include <string>

namespace EMPIRE {
class Message;
class SectionMesh : public FEMesh {
public:
    SectionMesh(std::string _name, int _numNodes, int _numElems, bool _triangulateAll = false);
    virtual ~SectionMesh();

    int getNumSections();
    int getNumRootSectionNodes();
    int getNumNormalSectionNodes();
    int getNumTipSectionNodes();
    const double *getTranslationGlobal2Root();
    const double *getRotationGlobal2Root();


    void setNumSections(int _numSections);
    void setNumRootSectionNodes(int _numRootSectionNodes);
    void setNumNormalSectionNodes(int _numNormalSectionNodes);
    void setNumTipSectionNodes(int _numTipSectionNodes);
    void setTranslationGlobal2Root(double *_translationGlobal2Root);
    void setRotationGlobal2Root(double *_rotationGlobal2Root);

private:
    int numSections;
    int numRootSectionNodes;
    int numNormalSectionNodes;
    int numTipSectionNodes;
    double *rotationGlobal2Root;
    double *translationGlobal2Root;
};
/***********************************************************************************************
 * \brief Allows for nice debug output later
 * \author Stefan Sicklinger
 ***********/
Message &operator<<(Message &message, SectionMesh &mesh);

} /* namespace EMPIRE */

#endif /* SECTIONMESH_H_ */
