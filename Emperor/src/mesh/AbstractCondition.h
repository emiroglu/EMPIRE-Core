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
#ifndef ABSTRACTCONDITION_H_
#define ABSTRACTCONDITION_H_

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class Message;
class IGAPatchSurface;
/********//**
 * \brief Class AbstractCondition is the superclass of all other conditions
 ***********/
class AbstractCondition {

protected:

    /// ID of the condition
    int ID;

public:
    /***********************************************************************************************
     * \brief Constructor, allocating the storage of the mesh
     * \param[in] _ID ID of the condition
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    AbstractCondition(int _ID = 0);
    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    virtual ~AbstractCondition(){};

    /***********************************************************************************************
     * \brief Create GP Data
     * \author Andreas Apostolatos, Altug Emiroglu
     ***********/
    virtual void createGPData(const std::vector<IGAPatchSurface*>& _surfacePatches){};

    /// type of the condition
    EMPIRE_Condition_type type;
};

} /* namespace EMPIRE */

#endif /* ABSTRACTCONDITION_H_ */
