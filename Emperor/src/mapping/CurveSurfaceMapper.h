#ifndef CURVESURFACEMAPPER_H_
#define CURVESURFACEMAPPER_H_

#include "AbstractMapper.h"
#include <map>

namespace EMPIRE {

class KinematicMotion;

class CurveSurfaceMapper : public AbstractMapper {
public:
    CurveSurfaceMapper(int _curveNumNodes, int _curveNumElements, const double *_curveNodeCoors,
            const int *_curveNodeIDs, const int *_curveElems, int _surfaceNumNodes,
            const double *_surfaceNodeCoors, int _surfaceNumSections, int _surfaceNumRootSectionNodes,
            int _surfaceNumNormalSectionNodes, int _surfaceNumTipSectionNodes,
            const double *rotation_O_Q, const double *translation_O_Q);

    virtual ~CurveSurfaceMapper();
    /***********************************************************************************************
     * \brief Map deformation from curve to surface
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *curveDispRot, double *surfaceDisp);
    /***********************************************************************************************
     * \brief Map loads form surface to curve
     * \param[in] surfaceForce forces on surface nodes
     * \param[out] curveForceMoment forces and moments on curve nodes
     * \author Tianyang Wang
     * ***********/
    void conservativeMapping(const double *surfaceForce, double *curveForceMoment);
private:
    /// number of nodes of a curve/beam
    int curveNumNodes;
    /// number of elements of a curve/beam
    int curveNumElements;
    //const double *curveNodeCoors;
    //const int *curveNodeIDs;
    /// connectivity table of curve/beam elements
    const int *curveElems;
    /// number of surface nodes
    int surfaceNumNodes;
    /// node coordinates of surface nodes
    const double *surfaceNodeCoors;
    /// number of sections of the surface
    int surfaceNumSections;
    /// number of root section nodes of the surface
    int surfaceNumRootSectionNodes;
    /// number of normal section nodes of the surface
    int surfaceNumNormalSectionNodes;
    /// number of the tip section nodes of the surface
    int surfaceNumTipSectionNodes;
    /// After sorting the surface nodes into sections, map from sorted position to its position in node array
    int *sortedPosToUnsortedPos;
    /// Map a curve/beam node ID to the its position in node array
    std::map<int, int> *curveNodeIDToPos;
    /// Map a section to the corresponding curve/beam element
    int *sectionToCurveElem;
    /// cross point between section and curve/beam
    double *sectionP;
    /// Shape functions and their derivatives of a section in the corresponding curve/beam element
    double *shapeFuncOfSection;
    /// Rotation from the global system the curve/beam element local system
    KinematicMotion **ROT_O_ELEM;
    /// unit test class
    friend class TestCurveSurfaceMapper;
};

} /* namespace EMPIRE */

#endif /* CURVESURFACEMAPPER_H_ */
