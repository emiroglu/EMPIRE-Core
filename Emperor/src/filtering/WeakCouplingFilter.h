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
#ifndef WEAKCOUPLINGFILTER_H_
#define WEAKCOUPLINGFILTER_H_

#include "AbstractFilter.h"
#include <assert.h>
#include <stdlib.h>

namespace EMPIRE {
/********//**
 * \brief Class WeakCouplingFilter predicts the traction force at the interface following:
 * master thesis of Ivan Hanzlicek, 2014
 ***********/
class WeakCouplingFilter : public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _beta for weak coupling weighting
     * \author Ivan Hanzlicek
     ***********/
	WeakCouplingFilter(double _beta);
    /***********************************************************************************************
     * \brief Destructor
     * \author Ivan Hanzlicek
     ***********/
    virtual ~WeakCouplingFilter();
    /***********************************************************************************************
     * \brief Filtering
     * \author Ivan Hanzlicek
     ***********/
    void filtering();
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Ivan Hanzlicek
     ***********/
    void init();
private:
    // weighting parameter
    const double beta;
    // data storage of oldest traction forces (t1 -> t2 -> current time step)
    double* t1;
    // data storage of old traction forces
    double* t2;
    // data storage of the traction force predictor
    double* tp;
};

} /* namespace EMPIRE */
#endif /* WEAKCOUPLINGFILTER_H_ */
