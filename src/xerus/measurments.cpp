// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2015 Benjamin Huber and Sebastian Wolf. 
// 
// Xerus is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// Xerus is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with Xerus. If not, see <http://www.gnu.org/licenses/>.
//
// For further information on Xerus visit https://libXerus.org 
// or contact us at contact@libXerus.org.

/**
 * @file
 * @brief Implementation of the measurment classes class.
 */

#include <xerus/measurments.h>


namespace xerus {
	// --------------------- SinglePointMeasurment -----------------
	
	SinglePointMeasurment::SinglePointMeasurment(const std::vector<size_t>& _positions, const value_t _value) : positions(_positions), value(_value) {}
	
	bool SinglePointMeasurment::Comparator::operator()(const SinglePointMeasurment &_lhs, const SinglePointMeasurment &_rhs) const {
		REQUIRE(_lhs.positions.size() == _rhs.positions.size(), "");
		for (size_t i=0; i<split_position-1 && i<_lhs.positions.size(); ++i) {
			if (_lhs.positions[i] < _rhs.positions[i]) return true;
			if (_lhs.positions[i] > _rhs.positions[i]) return false;
		}
		for (size_t i=_lhs.positions.size(); i>split_position; ++i) {
			if (_lhs.positions[i-1] < _rhs.positions[i-1]) return true;
			if (_lhs.positions[i-1] > _rhs.positions[i-1]) return false;
		}
		return false; // equality
	}
	
}
