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
 * @brief Implementation of the IHT variants.
 */

#include <xerus.h>

namespace xerus {
	double iht(TTTensor& _x, const SinglePointMeasurmentSet& _measurments, PerformanceData& _perfData) {
		const size_t numMeasurments = _measurments.size();
		const size_t degree = _x.degree();
		const size_t USER_MEASUREMENTS_PER_ITR = numMeasurments;
		

		TTTensor largeX(_x);
		// 		SinglePointMeasurmentSet currentValues(_measurments);
		std::vector<value_t> currentValues(numMeasurments);
		std::vector<size_t> measurementOrder(numMeasurments);
		std::iota(measurementOrder.begin(), measurementOrder.end(), 0);
		
		_perfData.start();
		
		Index i1, i2, i3, i4, i5;
		
		value_t normMeasuredValues = 0.0;
		for(size_t i = 0; i < numMeasurments; ++i) {
			normMeasuredValues += misc::sqr(_measurments.measuredValues[i]);
		}
		normMeasuredValues = std::sqrt(normMeasuredValues);
		
		std::mt19937_64 rnd(0x5EED0);
		
		for(size_t iteration = 0; iteration < 1000000; ++iteration) {
// 			_x.measure(currentValues);
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				currentValues[i] = _x[_measurments.positions[i]];
			}
			
// 			std::shuffle(measurementOrder.begin(), measurementOrder.end(), rnd);
			std::sort(measurementOrder.begin(), measurementOrder.end(), [&](size_t a, size_t b){
				return std::abs(_measurments.measuredValues[a] - currentValues[a]) > std::abs(_measurments.measuredValues[b] - currentValues[b]);
			});
			
			
			// Calculate residual for perfdata
			if (_perfData) {
				_perfData.stop_timer();
				value_t residual = 0;
				for(size_t i = 0; i < numMeasurments; ++i) {
					residual += misc::sqr(_measurments.measuredValues[i] - currentValues[i]);
				}
				_perfData.continue_timer();
				_perfData.add(iteration, std::sqrt(residual)/normMeasuredValues, _x, 0);
			}
			
			
			// Build the largeX
			for(size_t d = 0; d < degree; ++d) {
				Tensor& currComp = _x.component(d);
				SparseTensor newComp({d == 0 ? 1 : (currComp.dimensions[0]+USER_MEASUREMENTS_PER_ITR), currComp.dimensions[1], d == degree-1 ? 1 : (currComp.dimensions[2]+USER_MEASUREMENTS_PER_ITR)});
				
				// Copy dense part
				for(size_t r1 = 0; r1 < currComp.dimensions[0]; ++r1) {
					for(size_t n = 0; n < currComp.dimensions[1]; ++n) {
						for(size_t r2 = 0; r2 < currComp.dimensions[2]; ++r2) {
							newComp[{r1, n, r2}] = currComp[{r1, n, r2}];
						}
					}
				}
				
				// Copy sparse part
				if (d==0) {
					for(size_t i = 0; i < USER_MEASUREMENTS_PER_ITR; ++i) {
						newComp[{0, _measurments.positions[measurementOrder[i]][d], i+currComp.dimensions[2]}] = (_measurments.measuredValues[measurementOrder[i]] - currentValues[measurementOrder[i]]);
					}
				} else if (d!=degree-1) {
					for(size_t i = 0; i < USER_MEASUREMENTS_PER_ITR; ++i) {
						newComp[{i + currComp.dimensions[0], _measurments.positions[measurementOrder[i]][d], i+currComp.dimensions[2]}] = 1.0;
					}
				} else {
					// d == degree-1
					for(size_t i = 0; i < USER_MEASUREMENTS_PER_ITR; ++i) {
						newComp[{i + currComp.dimensions[0], _measurments.positions[measurementOrder[i]][d], 0}] = 1.0;
					}
				}
				
				
				largeX.set_component(d, newComp);
			}
			FullTensor fullX(largeX);
			
			_x = TTTensor(fullX, 0.0, _x.ranks());
			
// 			decomposition_als(_x, fullX);
			
// 			largeX = TTTensor(fullX);
			/*
			FullTensor doubleComp, U, Vt;
			SparseTensor S;
			for(size_t d = 0; d+1 < degree; ++d) {
				doubleComp(i1^2 , i3^2) = largeX.get_component(d)(i1^2, i2) * largeX.get_component(d+1)(i2, i3^2); 
				(U(i1^2, i2), S(i2, i3), Vt(i3, i4^2)) = SVD(doubleComp(i1^2, i4^2), 2*_x.rank(d));
				largeX.set_component(d, U);
				doubleComp(i1, i2^2) = S(i1, i3) * Vt(i3, i2^2);
				largeX.set_component(d+1, doubleComp);
			}
			
			largeX.assume_core_position(degree-1);
			largeX.cannonicalize_left();
			largeX.round(_x.ranks());*/
// 			_x = largeX;
		}
		
		value_t residual = 0;
		for(size_t i = 0; i < numMeasurments; ++i) {
			residual += misc::sqr(_measurments.measuredValues[i] - currentValues[i]);
		}
		return std::sqrt(residual);
	}
}
