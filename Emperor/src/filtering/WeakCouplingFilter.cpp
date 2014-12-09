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
#include "WeakCouplingFilter.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"

namespace EMPIRE {

WeakCouplingFilter::WeakCouplingFilter(double _beta) :
		AbstractFilter(), beta(_beta) {
}

WeakCouplingFilter::~WeakCouplingFilter() {
}

void WeakCouplingFilter::filtering() {
	INDENT_OUT(1, "WeakCouplingFilter: Calculating the predictor ...", infoOut);

	EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
	if (IOType == EMPIRE_ConnectionIO_DataField) {
		DataField *inDataField = inputVec[0]->dataField;
		for (int i = 0; i < inDataField->numLocations * inDataField->dimension; i++) {
			t2[i] = beta * inDataField->data[i] + (1.0 - beta) * tp[i];
			tp[i] = 2.0 * t2[i] - t1[i];
			t1[i] = t2[i];
			inDataField->data[i] = tp[i];
		}

	} else if (IOType == EMPIRE_ConnectionIO_Signal) {
		Signal *inSignal = inputVec[0]->signal;
		for (int i = 0; i < inSignal->size; i++) {
			t2[i] = beta * inSignal->array[i] + (1.0 - beta) * tp[i];
			tp[i] = 2.0 * t2[i] - t1[i];
			t1[i] = t2[i];
			inSignal->array[i] = tp[i];
		}
	} else {
		assert(false);
	}
}

void WeakCouplingFilter::init() {
	assert(inputVec.size() == 1);
	assert(outputVec.size() == 1);
	assert(inputVec[0]->type == outputVec[0]->type);
	EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
	if (IOType == EMPIRE_ConnectionIO_DataField) {
		assert(inputVec[0]->dataField == outputVec[0]->dataField);
		DataField *dataField = inputVec[0]->dataField;
		t1 = new double[dataField->dimension * dataField->numLocations]();
		t2 = new double[dataField->dimension * dataField->numLocations]();
		tp = new double[dataField->dimension * dataField->numLocations]();
	} else if (IOType == EMPIRE_ConnectionIO_Signal) {
		assert(inputVec[0]->signal == outputVec[0]->signal);
		Signal *signal = inputVec[0]->signal;
		t1 = new double[signal->size]();
		t2 = new double[signal->size]();
		tp = new double[signal->size]();
	} else {
		assert(false);
	}
}

} /* namespace EMPIRE */
