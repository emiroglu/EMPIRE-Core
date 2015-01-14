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

#ifndef BOUNDINGBOX_H_
#define BOUNDINGBOX_H_

#include "Message.h"

namespace EMPIRE {

class BoundingBox {
public:
	BoundingBox();
	virtual ~BoundingBox();

	virtual bool isPointInside(const double* P, double offset=0) const = 0;
};

class AABB : public BoundingBox {
public:
	AABB();
	virtual ~AABB();

	bool isPointInside(const double* P, double offset=0) const;

	inline double getXmin() const {return box[0];}
	inline double getXmax() const {return box[1];}
	inline double getYmin() const {return box[2];}
	inline double getYmax() const {return box[3];}
	inline double getZmin() const {return box[4];}
	inline double getZmax() const {return box[5];}
	inline void setXmin(double xmin) {box[0] = xmin;}
	inline void setXmax(double xmax) {box[1] = xmax;}
	inline void setYmin(double ymin) {box[2] = ymin;}
	inline void setYmax(double ymax) {box[3] = ymax;}
	inline void setZmin(double zmin) {box[4] = zmin;}
	inline void setZmax(double zmax) {box[5] = zmax;}
	inline bool isComputed() {return flagComputed;}
	inline void isComputed(bool flag) {flagComputed=flag;}

	double& operator[](int id) {return box[id];}
	const double& operator[](int id) const {return box[id];}
private:
	double box[6];
	bool flagComputed;
};

/***********************************************************************************************
 * \brief Output bounding box
 * \author Tianyang Wang
 ***********/
Message &operator<<(Message &message, AABB &aabb);


} /* namespace EMPIRE */
#endif /* BOUNDINGBOX_H_ */
