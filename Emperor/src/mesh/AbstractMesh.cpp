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
#include "AbstractMesh.h"
#include "Message.h"
#include "DataField.h"
#include <iostream>
#include <assert.h>

using namespace std;

namespace EMPIRE {

AbstractMesh::AbstractMesh(std::string _name, EMPIRE_Mesh_type _type) :
        name(_name), type(_type) {
}

AbstractMesh::~AbstractMesh() {
    for (map<string, DataField*>::iterator it = nameToDataFieldMap.begin();
            it != nameToDataFieldMap.end(); it++) {
        delete it->second;
    }
}

DataField *AbstractMesh::getDataFieldByName(std::string dataFieldName) {
    assert(nameToDataFieldMap.size()!=0);
    assert(nameToDataFieldMap.find(dataFieldName)!=nameToDataFieldMap.end());
    return nameToDataFieldMap.at(dataFieldName);
}

} /* namespace EMPIRE */
