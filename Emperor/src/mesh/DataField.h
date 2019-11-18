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
/***********************************************************************************************//**
 * \file DataField.h
 * This file holds the class DataField
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef DATAFIELD_H_
#define DATAFIELD_H_

#include <string>
#include "EMPEROR_Enum.h"
#include "Message.h"

namespace EMPIRE {

/********//**
 * \brief Class DataField stores a data field (vector/scalar) on a mesh
 ***********/
class DataField {
public:
    /***********************************************************************************************
     * \brief Constructor, allocating the storage instead of initialize everything
     * \param[in] _name name of the data field
     * \param[in] _numLocations number of locations (e.g. nodes or Gauss points)
     * \param[in] _dimension could be vector or scalar
     * \author Tianyang Wang
     ***********/
    DataField(std::string _name, EMPIRE_DataField_location _location, int _numLocations,
            EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~DataField();
    /***********************************************************************************************
     * \brief Method to write the data field out to a ASCII file
     * \param[in] _name name of the file to which it should be written
     * \param[in] _header Header which should be put on top of the file
     * \param[in] _footer Header which should be put on top of the file
     * \author Aditya Ghantasala
     ***********/
     void writeToFile(std::string _name,std::string _header, std::string _footer);

    /// name of the data field
    const std::string name;
    /// where the data are located, at node or at element centroid
    const EMPIRE_DataField_location location;
    /// number of locations (e.g. nodes or Gauss points)
    const int numLocations;
    /// dimension could be vector or scalar
    const EMPIRE_DataField_dimension dimension;
    /// typeOfQuantity could be field or fieldIntegral
    const EMPIRE_DataField_typeOfQuantity typeOfQuantity;
    /// pointer to the data storage
    double * const data;

};

/***********************************************************************************************
 * \brief Allows for nice debug output later
 * \author Stefan Sicklinger
 ***********/
Message &operator<<(Message &message, DataField &dataField);

} /* namespace EMPIRE */
#endif /* DATAFIELD_H_ */
