#ifndef CURVESURFACECOROTATE3DMAPPER_H_
#define CURVESURFACECOROTATE3DMAPPER_H_

#include "AbstractCurveSurfaceMapper.h"
#include <map>

namespace EMPIRE {

class KinematicMotion;

class CurveSurfaceCorotate3DMapper : public AbstractCurveSurfaceMapper {
public:
    CurveSurfaceCorotate3DMapper(int _curveNumNodes, int _curveNumElements, const double *_curveNodeCoors,
            const int *_curveNodeIDs, const int *_curveElems, int _surfaceNumNodes,
            const double *_surfaceNodeCoors, int _surfaceNumSections, int _surfaceNumRootSectionNodes,
            int _surfaceNumNormalSectionNodes, int _surfaceNumTipSectionNodes,
            const double *rotation_O_Q, const double *translation_O_Q);

    virtual ~CurveSurfaceCorotate3DMapper();
    /***********************************************************************************************
     * \brief Map deformation from curve to surface
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *curveDispRot, double *surfaceDisp);

};


} /* namespace EMPIRE */

#endif /* CURVESURFACECOROTATE3DMAPPER_H_ */
