#ifndef CURVESURFACELINEARMAPPER_H_
#define CURVESURFACELINEARMAPPER_H_

#include "AbstractCurveSurfaceMapper.h"
#include <map>

namespace EMPIRE {

class KinematicMotion;

class CurveSurfaceLinearMapper : public AbstractCurveSurfaceMapper {
public:
    CurveSurfaceLinearMapper(int _curveNumNodes, int _curveNumElements, const double *_curveNodeCoors,
            const int *_curveNodeIDs, const int *_curveElems, int _surfaceNumNodes,
            const double *_surfaceNodeCoors, int _surfaceNumSections, int _surfaceNumRootSectionNodes,
            int _surfaceNumNormalSectionNodes, int _surfaceNumTipSectionNodes,
            const double *rotation_O_Q, const double *translation_O_Q);

    virtual ~CurveSurfaceLinearMapper();
    /***********************************************************************************************
     * \brief Map deformation from curve to surface
     * \param[in] curveDispRot displacements and rotations on curve nodes
     * \param[out] surfaceDisp displacements on surface nodes
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *curveDispRot, double *surfaceDisp);
private:
    /// length of curve elements
    double *curveElemLength;
};

} /* namespace EMPIRE */

#endif /* CURVESURFACELINEARMAPPER_H_ */
