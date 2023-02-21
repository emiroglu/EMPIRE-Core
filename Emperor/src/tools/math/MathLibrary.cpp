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

#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>

namespace EMPIRE {
namespace MathLibrary {

// only very very important stuff goes here
int ceil(double _num) {
    int inum = (int)_num;
    if (_num == (double)inum) {
        return inum;
    }
    return inum + 1;
}

void sortRemoveDuplicates(std::vector<double>& _vec) {
    if (!_vec.empty()) {
        std::sort(_vec.begin(), _vec.end());
        _vec.erase(std::unique(_vec.begin(), _vec.end()), _vec.end());
    }
}

void mergeSortRemoveDuplicates(std::vector<double>& _vec0, std::vector<double> _vec1, std::vector<double> _vec2) {

    if (!_vec1.empty() && !_vec2.empty()) {
        sortRemoveDuplicates(_vec1);
        sortRemoveDuplicates(_vec2);
        _vec0.resize(_vec1.size()+_vec2.size());
        std::merge(_vec1.begin(), _vec1.end(), _vec2.begin(), _vec2.end(), _vec0.begin());
        sortRemoveDuplicates(_vec0);
    } else if (!_vec1.empty() && _vec2.empty()) {
        sortRemoveDuplicates(_vec1);
        _vec0 = _vec1;
    } else if (!_vec2.empty() && _vec1.empty()) {
        sortRemoveDuplicates(_vec2);
        _vec0 = _vec2;
    }
}

void sortRemoveSimilar(std::vector<double>& _vec, double _tol) {

    if (!_vec.empty()) {
        std::sort(_vec.begin(), _vec.end());
        std::vector<double>::iterator iVec = _vec.begin();
        while (iVec != _vec.end()-1 ) {
            if (fabs (*iVec - *(iVec+1)) < _tol ){
                _vec.erase(iVec);
                iVec = _vec.begin();
            } else iVec++;
        }
    }
}

} /* namespace Math */
} /* namespace EMPIRE */
