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
#ifndef SCALINGFILTER_H_
#define SCALINGFILTER_H_

#include "AbstractFilter.h"
#include <assert.h>
#include <stdlib.h>

namespace EMPIRE {
/********//**
 * \brief Class ScalingFilter scales signals and fields by a constant factor
 * \author Tianyang Wang
 ***********/
class ScalingFilter : public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _factor for scaling
     * \author Tianyang Wang
     ***********/
    ScalingFilter(double _factor);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~ScalingFilter();
    /***********************************************************************************************
     * \brief Filtering
     * \author Tianyang Wang
     ***********/
    void filtering();
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Tianyang Wang
     ***********/
    void init();
private:
    /// scaling factor
    double factor;
};

} /* namespace EMPIRE */
#endif /* SCALINGFILTER_H_ */
