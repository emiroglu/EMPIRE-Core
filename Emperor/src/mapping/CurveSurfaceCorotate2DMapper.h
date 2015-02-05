#ifndef CURVESURFACECOROTATE2DMAPPER_H_
#define CURVESURFACECOROTATE2DMAPPER_H_

#include "AbstractCurveSurfaceMapper.h"
#include <map>

namespace EMPIRE {

class KinematicMotion;

class CurveSurfaceCorotate2DMapper: public AbstractCurveSurfaceMapper {
public:
    CurveSurfaceCorotate2DMapper(int _curveNumNodes, int _curveNumElements,
            const double *_curveNodeCoors, const int *_curveNodeIDs, const int *_curveElems,
            int _surfaceNumNodes, const double *_surfaceNodeCoors, int _surfaceNumSections,
            int _surfaceNumRootSectionNodes, int _surfaceNumNormalSectionNodes,
            int _surfaceNumTipSectionNodes, const double *rotation_O_Q,
            const double *translation_O_Q);

    virtual ~CurveSurfaceCorotate2DMapper();
    /***********************************************************************************************
     * \brief Map deformation from curve to surface
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *curveDispRot, double *surfaceDisp);

private:
    /// Rotation from the root section system the curve/beam element local system
    KinematicMotion **ROT_Q_ELEM;
    KinematicMotion *KM_O_Q;
    double *angle_Q_ELEM;
    double *curveNodeCoorsInQ;
};

} /* namespace EMPIRE */

#endif /* CURVESURFACECOROTATE2DMAPPER_H_ */
