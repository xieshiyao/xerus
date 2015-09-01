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
		const size_t numMeasurments = _measurments.size();
		
		// Calculate the norm of all measured entries.
		value_t normMeasuredValues = 0;
		for(const SinglePointMeasurment& measurement : _measurments) {
			normMeasuredValues += misc::sqr(measurement.value);
		}
		normMeasuredValues = std::sqrt(normMeasuredValues);
		
		Index i1, i2, i3, i4, r1, r2;
		
		// Create the stack
		std::vector<std::vector<FullTensor>> stack(degree+2, std::vector<FullTensor>(numMeasurments, Tensor::ones({1})));
		
		value_t residual=1;
		
		for(size_t iteration = 0; iteration < maxInterations; ++iteration) {
			// Move core back to position zero
			_x.move_core(0, true);
			
			// Rebuild the lower part of the stack
			for(size_t corePosition = degree - 1; corePosition > 0; --corePosition) {
				
				// Prepare the single slates of the current component
				std::vector<FullTensor> fixedComponents(_x.dimensions[corePosition], FullTensor(_x.get_component(corePosition)));
				for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
					fixedComponents[i].fix_slate(1, i);
				}
				
				// Reset the FullTensors
				for(size_t i = 0; i < numMeasurments; ++i) {
					stack[corePosition+1][i](r2) = fixedComponents[_measurments[i].positions[corePosition]](r2, r1) * stack[corePosition+2][i](r1);
				}
			}
			
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				// Prepare the single slates of the current component
				std::vector<FullTensor> fixedComponents(_x.dimensions[corePosition], FullTensor(_x.get_component(corePosition)));
				for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
					fixedComponents[i].fix_slate(1, i);
				}
				
				FullTensor deltaPlus({_x.get_component(corePosition).dimensions[1], _x.get_component(corePosition).dimensions[0], _x.get_component(corePosition).dimensions[2]});
				FullTensor entryAddition({_x.get_component(corePosition).dimensions[0], _x.get_component(corePosition).dimensions[2]});
				FullTensor currentValue;
				residual = 0;
				
				for(size_t i = 0; i < numMeasurments; ++i) {
					entryAddition(r1, r2) = stack[corePosition][i](r1) * stack[corePosition+2][i](r2);
					currentValue() = entryAddition(r1, r2) * fixedComponents[_measurments[i].positions[corePosition]](r1, r2);
					const value_t currentDifference = _measurments[i].value-currentValue[0];
					misc::array_add(deltaPlus.data.get()+_measurments[i].positions[corePosition]*entryAddition.size, currentDifference, entryAddition.data.get(), entryAddition.size);
					residual += misc::sqr(currentDifference);
				}
				
				residual = sqrt(residual)/normMeasuredValues;
				LOG(ADF, iteration << " " << residual);
				
				// Update current component
				_x.component(corePosition)(r1, i1, r2) = _x.component(corePosition)(r1, i1, r2) + deltaPlus(i1, r1, r2);
				
				if(corePosition+1 < degree) {
					_x.move_core(corePosition+1, true);
					
					// We need updated fixedComponents
					for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
						fixedComponents[i](r1, r2) = _x.component(corePosition)(r1, i, r2);
					}
					
					// Update the stack
					for(size_t i = 0; i < numMeasurments; ++i) {
						stack[corePosition+1][i](r2) = stack[corePosition][i](r1) * fixedComponents[_measurments[i].positions[corePosition]](r1, r2);
					}
				}
			}
		}
		return residual;
	}
	
	/*
	
	double ADFVariant::solve(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments) const {
		REQUIRE(_x.is_valid_tt(), "");
		REQUIRE(_measurments.size() > 0, "Need at very least one measurment.");
		REQUIRE(_measurments[0].positions.size() == _x.degree(), "Measurment order must coincide with _x order.");
		
		const size_t degree = _x.degree();
		Index i1, i2, i3, i4, r1, r2;
		
		std::vector<SinglePointMeasurment> ourMeasurments(_measurments);
		
		std::vector<std::vector<FullTensor>> upperStack, lowerStack;
		
		value_t normMeasuredValues = 0;
		for(const SinglePointMeasurment& measurement : ourMeasurments) {
			normMeasuredValues += misc::sqr(measurement.value);
		}
		normMeasuredValues = std::sqrt(normMeasuredValues);
		
		std::vector<TensorNetwork> stack;
		value_t residual=1;
		
		for(size_t iteration = 0; iteration < maxInterations; ++iteration) {
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				_x.move_core(corePosition, true);
				// Build initial stack 
				stack.resize(1);
				std::pair<TensorNetwork, TensorNetwork> choppedTN = _x.chop(corePosition);
				stack.back()(i1^(corePosition+1), i2&(corePosition+1)) = choppedTN.first(i1&0)*choppedTN.second(i2&0);
				stack.back().reduce_representation();
				
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
				FullTensor componentUpdate({_x.get_component(corePosition).dimensions[1], _x.get_component(corePosition).dimensions[0], _x.get_component(corePosition).dimensions[2]});
				const size_t blockSize = _x.get_component(corePosition).dimensions[0] * _x.get_component(corePosition).dimensions[2];
				
				std::vector<size_t> previousPosition;
				bool firstStep = true;
				residual=0;
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
					
// 					size_t bla = 0;
// 					for(const TensorNetwork& tn : stack) {
// 						tn.draw(std::string("stack_lvl") + misc::to_string(bla++));
// 					}
					
					// Calculate current value at measured position
					FullTensor missingComp;
					missingComp(i1, i2) = _x.get_component(corePosition)(i1, measurment.positions[corePosition], i2);
					const value_t factor = measurment.value - value_t(missingComp(i1^2)*stack.back()(i1^2));
					residual += misc::sqr(factor);
					
					FullTensor X;
					X = FullTensor(stack.back());
					REQUIRE(X.size == blockSize, "Internal Error.");
					
					misc::array_add(componentUpdate.data.get()+measurment.positions[corePosition]*blockSize, factor, X.data.get(), blockSize);
				}
				residual = std::sqrt(residual)/ normMeasuredValues;
				
				if (printProgress) {
					std::cout << "in iteration: " << iteration << " core pos: " << corePosition << " residual: " << residual << std::endl;
					std::cout << "\b                                                                                          \b";
				}
				
				// calculate ||Ay||^2 for step size
				std::vector<SinglePointMeasurment> measCopy(ourMeasurments);
				TTTensor changeDirection(_x);
				// Reshuffle the component update
				componentUpdate(i1,i2,i3) = componentUpdate(i2,i1,i3);
				changeDirection.set_component(corePosition, componentUpdate);
				changeDirection.measure(measCopy);
				value_t norm_ay = 0;
				for (const SinglePointMeasurment &m : measCopy) {
					norm_ay += misc::sqr(m.value);
				}
				value_t alpha = misc::sqr(frob_norm(componentUpdate))/norm_ay;
				
				_x.component(corePosition)(i1,i2,i3) = _x.component(corePosition)(i1,i2,i3) + alpha*componentUpdate(i1,i2,i3);
			}
			
			// Check for termination
			// TODO relative residual
			if (residual < convergenceEpsilon) {
				return residual;
			}
		}
		return residual;
	}
	*/
	
}
