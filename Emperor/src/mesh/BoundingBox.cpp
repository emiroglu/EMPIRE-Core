/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Fabien Pean, Munich
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

#include "BoundingBox.h"
#include <assert.h>

namespace EMPIRE {

BoundingBox::BoundingBox() {
	// TODO Auto-generated constructor stub
}

BoundingBox::~BoundingBox() {
	// TODO Auto-generated destructor stub
}

AABB::AABB():
		box(),
		flagComputed(false) {
	// TODO Auto-generated constructor stub
}

AABB::~AABB() {
	// TODO Auto-generated destructor stub
}

bool AABB::isPointInside(const double* P, double offset) const {
	assert(P != NULL);
	if(P[0] >= box[0] - offset && P[0] <= box[1] + offset
			&& P[1] >= box[2] - offset && P[1] <= box[3] + offset
			&& P[2] >= box[4] - offset && P[2] <= box[5] + offset)
		return true;
	return false;
}

Message &operator<<(Message &message, AABB &boundingBox) {
	using std::endl;
    message << "\t+" << "Axis Aligned Bounding box: " << endl;
    message << "\t\t+" << '\t' << "xmin: " << boundingBox[0] << endl;
    message << "\t\t+" << '\t' << "xmax: " << boundingBox[1] << endl;
    message << "\t\t+" << '\t' << "ymin: " << boundingBox[2] << endl;
    message << "\t\t+" << '\t' << "ymax: " << boundingBox[3] << endl;
    message << "\t\t+" << '\t' << "zmin: " << boundingBox[4] << endl;
    message << "\t\t+" << '\t' << "zmax: " << boundingBox[5] << endl;
    message() << "\t+" << "---------------------------------" << endl;
    return message;
}

} /* namespace EMPIRE */
