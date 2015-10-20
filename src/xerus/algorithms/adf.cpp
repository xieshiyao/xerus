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
#include <xerus/sparseTensor.h>
#include <xerus/misc/check.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/blasLapackWrapper.h>

namespace xerus {
	
	double calculate_norm_of_measured_values(const std::vector<value_t>& _measurments) {
		value_t normMeasuredValues = 0;
		for(const value_t measurement : _measurments) {
			normMeasuredValues += misc::sqr(measurement);
		}
		return std::sqrt(normMeasuredValues);
	}
	
	void construct_forward_stacks(std::unique_ptr<FullTensor[]>& stackSaveSlots, std::vector<bool>& forwardUpdates, std::unique_ptr<FullTensor*[]>& forwardStackMem, TTTensor& _x, const SinglePointMeasurmentSet& _measurments) {
		const size_t degree = _x.degree();
		const size_t numMeasurments = _measurments.size();
		
		// Temporary map
		std::vector<size_t> forwardCalculationMap(degree*numMeasurments);
		
		// Count how many FullTensors we need for the stacks
		size_t numUniqueStackEntries = 0;
	
		// TODO ensure that the measurments are ordered according to their position (atm this is a prerequisite)
		
		// Create the forward map
		for(size_t corePosition = 0; corePosition+1 < degree; ++corePosition) {
			forwardCalculationMap[0 + corePosition*numMeasurments] = 0;
			forwardUpdates[0 + corePosition*numMeasurments] = true;
			++numUniqueStackEntries;
		}
		
		for(size_t i = 1; i < numMeasurments; ++i) {
			size_t corePosition = 0;
			for( ; corePosition+1 < degree && _measurments.positions[i][corePosition] == _measurments.positions[i-1][corePosition]; corePosition++) {
				forwardCalculationMap[i + corePosition*numMeasurments] = forwardCalculationMap[i-1 + corePosition*numMeasurments];
				forwardUpdates[i + corePosition*numMeasurments] = false;
			}
			for( ; corePosition+1 < degree; ++corePosition) {
				forwardCalculationMap[i + corePosition*numMeasurments] = i;
				forwardUpdates[i + corePosition*numMeasurments] = true;
				++numUniqueStackEntries;
			}
		}
		
		// Create the stack
		stackSaveSlots.reset(new FullTensor[numUniqueStackEntries+1]); // +1 for the special positions -1 and degree.
		size_t usedSlots = 0; // Zero is reserved for the the position -1 and degree stacks
		stackSaveSlots[usedSlots++] = Tensor::ones({1});
		
		for(size_t i = 0; i < numMeasurments; ++i) {
			forwardStackMem[i + 0] = &stackSaveSlots[0];
		}
			
		for(size_t corePosition = 0; corePosition+1 < degree; ++corePosition) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(forwardCalculationMap[i + corePosition*numMeasurments] == i) {
					stackSaveSlots[usedSlots].reset({_x.rank(corePosition)}, DONT_SET_ZERO());
					forwardStackMem[i + (corePosition+1)*numMeasurments] = &stackSaveSlots[usedSlots++];
				} else {
					forwardStackMem[i + (corePosition+1)*numMeasurments] = forwardStackMem[forwardCalculationMap[i + corePosition*numMeasurments] + (corePosition+1)*numMeasurments];
				}
			}
		}
		
		REQUIRE(usedSlots == numUniqueStackEntries+1, "Internal Error.");
		
		LOG(ADF, "We have " << numUniqueStackEntries << " unique forward stack entries. There are " << numMeasurments*(degree-1)+1 << " virtual forward stack entries.");
	}
	
	void construct_backward_stacks(std::unique_ptr<FullTensor[]>& stackSaveSlots, std::vector<bool>& backwardUpdates, std::unique_ptr<FullTensor*[]>& backwardStackMem, TTTensor& _x, const SinglePointMeasurmentSet& _measurments) {
		const size_t degree = _x.degree();
		const size_t numMeasurments = _measurments.size();
		
		// Temporary map
		std::vector<size_t> backwardCalculationMap(degree*numMeasurments);
		
		// Count how many FullTensors we need for the stacks
		size_t numUniqueStackEntries = 0;
	
		// TODO ensure that the measurments are ordered according to their position (atm this is a prerequisite)
		
		// Create the backward map
		std::vector<size_t> reorderedMeasurments(numMeasurments);
		std::iota(reorderedMeasurments.begin(), reorderedMeasurments.end(), 0);
		
		std::sort(reorderedMeasurments.begin(), reorderedMeasurments.end(), 
			[&](const size_t _a, const size_t _b) {
				for (size_t j = degree; j > 0; --j) {
					if (_measurments.positions[_a][j-1] < _measurments.positions[_b][j-1]) return true;
					if (_measurments.positions[_a][j-1] > _measurments.positions[_b][j-1]) return false;
				}
				return false;
			}
		);
		
		for(size_t corePosition = 1; corePosition < degree; ++corePosition) {
			const size_t realId = reorderedMeasurments[0];
			backwardCalculationMap[realId + corePosition*numMeasurments] = realId;
			backwardUpdates[realId + corePosition*numMeasurments] = true;
			++numUniqueStackEntries;
		}
		
		for(size_t i = 1; i < numMeasurments; ++i) {
			const size_t realId = reorderedMeasurments[i];
			const size_t realPreviousId = reorderedMeasurments[i-1];
			
			size_t corePosition = degree-1;
			for( ; corePosition > 0 && _measurments.positions[realId][corePosition] == _measurments.positions[realPreviousId][corePosition]; --corePosition) {
				size_t otherId = realPreviousId;
				while(true) {
					if( otherId < realId ) {
						backwardCalculationMap[realId + corePosition*numMeasurments] = backwardCalculationMap[otherId + corePosition*numMeasurments];
						backwardUpdates[realId + corePosition*numMeasurments] = false;
						break;
					} else if(otherId == backwardCalculationMap[otherId + corePosition*numMeasurments]) {
						backwardCalculationMap[otherId + corePosition*numMeasurments] = realId;
						backwardUpdates[otherId + corePosition*numMeasurments] = false;
						backwardCalculationMap[realId + corePosition*numMeasurments] = realId;
						backwardUpdates[realId + corePosition*numMeasurments] = true;
						break;
					} else if( realId < backwardCalculationMap[otherId + corePosition*numMeasurments]) {
						size_t nextOther = backwardCalculationMap[otherId + corePosition*numMeasurments];
						backwardCalculationMap[otherId + corePosition*numMeasurments] = realId;
						REQUIRE(backwardUpdates[otherId + corePosition*numMeasurments] == false, "IE");
						otherId = nextOther;
					} else {
						otherId = backwardCalculationMap[otherId + corePosition*numMeasurments];
					}
				}
			}
			
			for( ; corePosition > 0; --corePosition) {
				backwardCalculationMap[realId + corePosition*numMeasurments] = realId;
				backwardUpdates[realId + corePosition*numMeasurments] = true;
				++numUniqueStackEntries;
			}
		}
		
		// Create the stack
		stackSaveSlots.reset(new FullTensor[numUniqueStackEntries+1]); // +1 for the special positions -1 and degree.
		size_t usedSlots = 0; // Zero is reserved for the the position -1 and degree stacks
		stackSaveSlots[usedSlots++] = Tensor::ones({1});
		
		for(size_t i = 0; i < numMeasurments; ++i) {
			backwardStackMem[i + (degree+1)*numMeasurments] = &stackSaveSlots[0];
		}
		
		for(size_t corePosition = 1; corePosition < degree; ++corePosition) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(backwardCalculationMap[i + corePosition*numMeasurments] == i) {
					stackSaveSlots[usedSlots].reset({_x.rank(corePosition-1)}, DONT_SET_ZERO());
					backwardStackMem[i + (corePosition+1)*numMeasurments] = &stackSaveSlots[usedSlots++];
				} else {
					backwardStackMem[i + (corePosition+1)*numMeasurments] = backwardStackMem[backwardCalculationMap[i + corePosition*numMeasurments] + (corePosition+1)*numMeasurments];
				}
			}
		}
		
		REQUIRE(usedSlots == numUniqueStackEntries+1, "Internal Error.");
		
		LOG(ADF, "We have " << numUniqueStackEntries << " unique backward stack entries. There are " << numMeasurments*(degree-1)+1 << " virtual backward stack entries.");
	}
	
	double ADFVariant::solve(TTTensor& _x, const SinglePointMeasurmentSet& _measurments) const {
		REQUIRE(_x.is_valid_tt(), "_x must be a valid TT-Tensor.");
		
		const size_t degree = _x.degree();
		const size_t numMeasurments = _measurments.size();
		
		REQUIRE(numMeasurments > 0, "Need at very least one measurment.");
		REQUIRE(_measurments.degree() == degree, "Measurment order must coincide with _x order.");
		
		const value_t normMeasuredValues = calculate_norm_of_measured_values(_measurments.measuredValues);
		
		// We want to construct two times numMeasurments stacks of size degree+2 containing the 
		// left (forward) and right (backward) pre contracted stack for each measurment. (technically numMeasurments stacks
		// would be sufficent, but there is nearly no overhead since its only pointers)
		
		// Bool arrays signaaling whether a specific stack entry is unique and has to be calculated
		std::vector<bool> forwardUpdates(degree*numMeasurments);
		std::vector<bool> backwardUpdates(degree*numMeasurments);
		
		// Array that will contain all the unqiue stack entries.
		std::unique_ptr<FullTensor[]> forwardStackSaveSlots;
		std::unique_ptr<FullTensor[]> backwardStackSaveSlots;
		
		// Arrays mapping the stack entry to the corresponding unique entry in stackSaveSlots (NOTE that both allow position -1 and degree.
		std::unique_ptr<FullTensor*[]> forwardStackMem( new FullTensor*[numMeasurments*(degree+2)]);
		FullTensor* const * const forwardStack = forwardStackMem.get()+numMeasurments;
		
		std::unique_ptr<FullTensor*[]> backwardStackMem( new FullTensor*[numMeasurments*(degree+2)]);
		FullTensor* const * const backwardStack = backwardStackMem.get()+numMeasurments;
		
		// We need _x to be canonicalized in the sense that there is no edge with more than maximal rank. (NOTE move_core(0) is called either way at the begining of the iterations)
		_x.cannonicalize_left();
		
		construct_forward_stacks(forwardStackSaveSlots, forwardUpdates, forwardStackMem, _x, _measurments);
		construct_backward_stacks(backwardStackSaveSlots, backwardUpdates, backwardStackMem, _x, _measurments);
		
		

		Index i1, i2, r1, r2;
		
		value_t residual = 1.0, lastResidual = 1.0;
		size_t smallResidualCount = 0;
		
		std::vector<FullTensor> fixedComponents(misc::max(_x.dimensions));
		std::vector<value_t> currentDifferences(numMeasurments);
		
		for(size_t iteration = 0; iteration < maxInterations; ++iteration) {
			// Move core back to position zero
			_x.move_core(0, true);
			
			// Rebuild the lower part of the stack
			for(size_t corePosition = degree-1; corePosition > 0; --corePosition) {
				
				// Prepare the single slates of the current component 
				for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
					fixedComponents[i](r1, r2) = _x.get_component(corePosition)(r1, i, r2);
				}
				
				// Reset the FullTensors
				for(size_t i = 0; i < numMeasurments; ++i) {
					if(backwardUpdates[i + corePosition*numMeasurments]) {
						contract(*backwardStack[i + corePosition*numMeasurments], fixedComponents[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
					}
				}
			}
			
			// Sweep from the first to the last component
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				const size_t localN = _x.get_component(corePosition).dimensions[1];
				const size_t localLeftRank = _x.get_component(corePosition).dimensions[0];
				const size_t localRightRank = _x.get_component(corePosition).dimensions[2];
				
				// Prepare the single slates of the current component 
				for(size_t i = 0; i < localN; ++i) {
					fixedComponents[i](r1, r2) = _x.get_component(corePosition)(r1, i, r2);
				}
				
				std::vector<FullTensor> deltas(localN, FullTensor({localLeftRank, localRightRank}));
				FullTensor entryAddition({localLeftRank, localRightRank});
				FullTensor currentValue({});
				
				for(size_t i = 0; i < numMeasurments; ++i) {
					contract(entryAddition, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 0);
					contract(currentValue, entryAddition, false, fixedComponents[_measurments.positions[i][corePosition]], false, 2);
					currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
					deltas[_measurments.positions[i][corePosition]] += currentDifferences[i]*entryAddition;
				}
				
				// Calculate ||P(y)||^2 for each slice seperately, where y is the update direction.
				std::vector<value_t> PyPys(localN, 0.0);
				FullTensor halfPy({localLeftRank});
				for(size_t i = 0; i < numMeasurments; ++i) {
					contract(halfPy, deltas[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
					contract(currentValue, *forwardStack[i + (corePosition-1)*numMeasurments], false, halfPy, false, 1);
					PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
				}
				
				// Update each slice seperately
				for(size_t j = 0; j < localN; ++j) {
					// Calculate <P(y), P(X-B)> = ||deltaPlus||^2.
					const value_t PyR = misc::sqr(frob_norm(deltas[j]));
					
					// Update
					_x.component(corePosition)(r1, i1, r2) = _x.component(corePosition)(r1, i1, r2) + (PyR/PyPys[j])*Tensor::dirac({localN}, {j})(i1)*deltas[j](r1, r2);
				}
				
				
				// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
				if(corePosition+1 < degree) {
					_x.move_core(corePosition+1, true);
					
					// Prepare the single slates of the current component 
					for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
						fixedComponents[i](r1, r2) = _x.get_component(corePosition)(r1, i, r2);
					}
					
					// Update the stack
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(forwardUpdates[i + corePosition*numMeasurments]) {
							contract(*forwardStack[i + corePosition*numMeasurments] , *forwardStack[i + (corePosition-1)*numMeasurments], false, fixedComponents[_measurments.positions[i][corePosition]], false, 1);
						}
					}
				}
			}
			
			// Calculate the residual (actually the residual of the last iteration because this comes more or less for free).
			lastResidual = residual;
			residual = 0;
			for(size_t i = 0; i < numMeasurments; ++i) {
				residual += misc::sqr(currentDifferences[i]);
			}
			residual = std::sqrt(residual)/normMeasuredValues;
			
			if(residual/lastResidual > 1.0 - 1e-3) {
				smallResidualCount++;
			} else {
				smallResidualCount = 0;
			}
			
			LOG(ADF, "Itr: " << iteration << " Residual: " << std::scientific << residual << " Rel. Residual change: " << residual/lastResidual);
			
			if(residual <= convergenceEpsilon || smallResidualCount > 3) {
				return residual;
			}
		}
		return residual;
	}
}
