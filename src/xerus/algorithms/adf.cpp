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
 * @brief Implementation of the ADF variants.
 */

#include <xerus/algorithms/adf.h>
#include <xerus/misc/check.h>
#include <xerus/indexedTensor_TN_operators.h>

namespace xerus {
	
	double ADFVariant::solve(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments) const {
		REQUIRE(_x.is_valid_tt(), "");
		REQUIRE(_measurments.size() > 0, "Need at very least one measurment.");
		REQUIRE(_measurments[0].positions.size() == _x.degree(), "Measurment order must coincide with _x order.");
		
		const size_t degree = _x.degree();
		Index i1, i2, i3, i4, r1, r2;
		
		std::vector<SinglePointMeasurment> ourMeasurments(_measurments);
		
		_x.move_core(0, true);
		
		std::vector<TensorNetwork> stack;
		
		for(size_t iteration = 0; iteration < maxInterations; ++iteration) {
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				// Build initial stack 
				stack.resize(1);
				std::pair<TensorNetwork, TensorNetwork> choppedTN = _x.chop(corePosition);
				stack.back()(i1^(corePosition+1), i2&(corePosition+1)) = choppedTN.first(i1&0)*choppedTN.second(i2&0);
				
				// Build translation array
				std::vector<size_t> orderToPos;
				orderToPos.reserve(degree);
				for(size_t i = 0; i < corePosition; ++i) {
					orderToPos.emplace_back(i);
				}
				for(size_t i = degree; i > corePosition; --i) {
					orderToPos.emplace_back(i-1);
				}
				
				// Sort measurements according to core position.
				std::sort(ourMeasurments.begin(), ourMeasurments.end(), SinglePointMeasurment::Comparator(corePosition));
				
				// Create update tensor
				FullTensor updateDirection({_x.get_component(corePosition).dimensions[1], _x.get_component(corePosition).dimensions[0], _x.get_component(corePosition).dimensions[2]});
				const size_t blockSize = _x.get_component(corePosition).dimensions[0] * _x.get_component(corePosition).dimensions[2];
				
				std::vector<size_t> previousPosition;
				bool firstStep = true;
				for(const SinglePointMeasurment& measurment : ourMeasurments) {
					// Find the maximal recyclable stack position
					size_t rebuildIndex = ~0ul;
					if(firstStep) {
						firstStep = false;
						rebuildIndex = 0;
					} else {
						for(size_t i = 0; i < degree; ++i) {
							if(previousPosition[orderToPos[i]] != measurment.positions[orderToPos[i]]) {
								rebuildIndex = i;
								break;
							}
						}
						REQUIRE(rebuildIndex != ~0ul, "There were two identical measurements? pos: " << previousPosition);
					}
					previousPosition = measurment.positions;
					
					// Trash stack that is not needed anymore
					while(stack.size() > rebuildIndex+1) { stack.pop_back(); }
					
					// Rebuild stack
					for(size_t i = rebuildIndex; i < degree-1; ++i) {
						REQUIRE(orderToPos[i] != corePosition, "Internal Error.");
						
						if (orderToPos[i] < corePosition) {
							stack.emplace_back(stack.back());
							stack.back().fix_slate(0, measurment.positions[orderToPos[i]]);
							stack.back().reduce_representation();
						} else {
							stack.emplace_back(stack.back());
							stack.back().fix_slate(stack.back().degree()-1, measurment.positions[orderToPos[i]]);
							stack.back().reduce_representation();
						}
					}
					REQUIRE(stack.size() == degree && stack.back().degree() == 2, "Internal Error");
					
					size_t bla = 0;
					for(const TensorNetwork& tn : stack) {
						tn.draw(std::string("stack_lvl") + misc::to_string(bla++));
					}
					
					// Calculate current value at measured position
					FullTensor missingComp;
					missingComp(i1, i2) = _x.get_component(corePosition)(i1, measurment.positions[corePosition], i2);
					const value_t factor = measurment.value - value_t(missingComp(i1^2)*stack.back()(i1^2));
					
					FullTensor X;
					X = FullTensor(stack.back());
					REQUIRE(X.size == blockSize, "Internal Error.");
					
					misc::array_add(updateDirection.data.get()+measurment.positions[corePosition]*blockSize, factor, X.data.get(), blockSize);
					
				}
			}
			
			// Check for termination
		}
		return 0.0;
	}
	
}
