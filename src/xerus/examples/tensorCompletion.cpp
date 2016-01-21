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
	
	void inverse_index_norm(std::vector<SinglePointMeasurment> &_measurements, value_t _additiveConst) {
		for (SinglePointMeasurment &meas : _measurements) {
			value_t normSqr=0;
			for (size_t i=0; i<meas.positions.size(); ++i) {
				normSqr += misc::sqr(double(meas.positions[i])+_additiveConst);
			}
			meas.value = 1/std::sqrt(normSqr);
		}
	}
	
	
	
	void inverse_index_ratios(std::vector<SinglePointMeasurment> &_measurements, value_t _additiveConst) {
		for (SinglePointMeasurment &meas : _measurements) {
			value_t sum=0;
			for (size_t i=0; i<meas.positions.size()-1; ++i) {
				sum += (double(meas.positions[i]) + 1.0) / (double(meas.positions[i+1]) + _additiveConst);
			}
			meas.value = 1/(_additiveConst + sum);
		}
	}
	
}}}

