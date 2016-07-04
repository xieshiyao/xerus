// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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
 * @brief Implementation of the example measure sets for tensor completion.
 */

#include <xerus/examples/tensorCompletion.h>
#include <xerus/measurments.h>

namespace xerus { namespace examples { namespace completion {
	
	void inverse_index_norm(SinglePointMeasurementSet& _measurements, const value_t _additiveConst) {
		for (size_t i = 0; i < _measurements.size(); ++i) {
			value_t normSqr = 0;
			for (const auto idx : _measurements.positions[i]) {
				normSqr += misc::sqr(static_cast<value_t>(idx) + _additiveConst);
			}
			_measurements.measuredValues[i] = 1/std::sqrt(normSqr);
		}
	}
	
	
	void inverse_index_ratios(SinglePointMeasurementSet& _measurements, const value_t _additiveConst) {
		for (size_t i = 0; i < _measurements.size(); ++i) {
			value_t sum = 0;
			for (size_t j = 0; j+1 < _measurements.positions[i].size(); ++j) {
				sum += (static_cast<value_t>(_measurements.positions[i][j]) + 1.0) / (static_cast<value_t>(_measurements.positions[i][j+1]) + _additiveConst);
			}
			_measurements.measuredValues[i] = 1/(_additiveConst + sum);
		}
	}
} // namespace completion
} // namespace examples
} // namespace xerus

