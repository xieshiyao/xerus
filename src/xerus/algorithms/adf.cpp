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
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_TN_operators.h>

namespace xerus {
	
	inline double calculate_norm_of_measured_values(const std::vector<value_t>& _measurments) {
		value_t normMeasuredValues = 0;
		for(const value_t measurement : _measurments) {
			normMeasuredValues += misc::sqr(measurement);
		}
		return std::sqrt(normMeasuredValues);
	}
	
	void construct_stacks(std::unique_ptr<FullTensor[]>& _stackSaveSlot, std::vector<std::vector<size_t>>& _updates, std::unique_ptr<FullTensor*[]>& _stackMem, const TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const bool _forward) {
		const size_t degree = _x.degree();
		const size_t numMeasurments = _measurments.size();
		
		// Temporary map. For each stack entrie (i.e. measurement number + corePosition) gives a measurement number of a stack entrie (at same corePosition) that shall have an equal value (or its own number otherwise). 
		std::vector<size_t> calculationMap(degree*numMeasurments);
		
		// Count how many FullTensors we need for the stacks
		size_t numUniqueStackEntries = 0;
		
		// Create a reordering map
		std::vector<size_t> reorderedMeasurments(numMeasurments);
		std::iota(reorderedMeasurments.begin(), reorderedMeasurments.end(), 0);
		
		std::sort(reorderedMeasurments.begin(), reorderedMeasurments.end(), 
			[&](const size_t _a, const size_t _b) {
				if(_forward) {
					for (size_t j = 0; j < degree; ++j) {
						if (_measurments.positions[_a][j] < _measurments.positions[_b][j]) return true;
						if (_measurments.positions[_a][j] > _measurments.positions[_b][j]) return false;
					}
				} else {
					for (size_t j = degree; j > 0; --j) {
						if (_measurments.positions[_a][j-1] < _measurments.positions[_b][j-1]) return true;
						if (_measurments.positions[_a][j-1] > _measurments.positions[_b][j-1]) return false;
					}
				}
				LOG(fatal, "Measurments must not appear twice.");
				return false;
			}
		);
		
		// Create the entries for the first measurement (this one allways unqiue).
		for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
			const size_t realId = reorderedMeasurments[0];
			calculationMap[realId + corePosition*numMeasurments] = realId;
			++numUniqueStackEntries;
		} 
		
		for(size_t i = 1; i < numMeasurments; ++i) {
			const size_t realId = reorderedMeasurments[i];
			const size_t realPreviousId = reorderedMeasurments[i-1];
			
			size_t position = 0;
			size_t corePosition = _forward ? position : degree-1-position;
			for( ; position < degree && _measurments.positions[realId][corePosition] == _measurments.positions[realPreviousId][corePosition]; ++position, corePosition = _forward ? position : degree-1-position) {
				size_t otherId = realPreviousId;
				while(true) {
					if( otherId < realId ) {
						calculationMap[realId + corePosition*numMeasurments] = calculationMap[otherId + corePosition*numMeasurments];
						break;
					} else if(otherId == calculationMap[otherId + corePosition*numMeasurments]) {
						calculationMap[otherId + corePosition*numMeasurments] = realId;
						calculationMap[realId + corePosition*numMeasurments] = realId;
						break;
					} else if( realId < calculationMap[otherId + corePosition*numMeasurments]) {
						const size_t nextOther = calculationMap[otherId + corePosition*numMeasurments];
						calculationMap[otherId + corePosition*numMeasurments] = realId;
						otherId = nextOther;
					} else {
						otherId = calculationMap[otherId + corePosition*numMeasurments];
					}
				}
			}
			
			for( ; position < degree; ++position, corePosition = _forward ? position : degree-1-position) {
				calculationMap[realId + corePosition*numMeasurments] = realId;
				++numUniqueStackEntries;
			}
		}
		
		// Create the stack
		numUniqueStackEntries++; // +1 for the special positions -1 and degree.
		_stackSaveSlot.reset(new FullTensor[numUniqueStackEntries]); 
		size_t usedSlots = 0; 
		_stackSaveSlot[usedSlots++] = Tensor::ones({1}); // Special slot reserved for the the position -1 and degree stacks
		
		// NOTE that _stackMem contains (degree+2)*numMeasurments entries and has an offset of numMeasurments (to have space for corePosition -1).
		
		// Set links for the special entries -1 and degree
		for(size_t i = 0; i < numMeasurments; ++i) {
			_stackMem[i] = &_stackSaveSlot[0];
			_stackMem[i + (degree+1)*numMeasurments] = &_stackSaveSlot[0];
		}
		
		for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(calculationMap[i + corePosition*numMeasurments] == i) {
					_updates[corePosition].emplace_back(i);
					_stackSaveSlot[usedSlots].reset({_x.rank(corePosition - (_forward ? 0 : 1))}, DONT_SET_ZERO());
					_stackMem[i + (corePosition+1)*numMeasurments] = &_stackSaveSlot[usedSlots++];
				} else {
					_stackMem[i + (corePosition+1)*numMeasurments] = _stackMem[calculationMap[i + corePosition*numMeasurments] + (corePosition+1)*numMeasurments];
				}
			}
		}
		
		REQUIRE(usedSlots == numUniqueStackEntries, "Internal Error.");
		LOG(ADF, "We have " << numUniqueStackEntries << " unique stack entries. There are " << numMeasurments*degree+1 << " virtual stack entries.");
	}
	
	
	inline void rebuild_backward_stack(FullTensor* const * const _backwardStack, const TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const std::vector<std::vector<size_t>>& _backwardUpdates) {
		const size_t numMeasurments = _measurments.size();
		
		const Index r1, r2;
		std::vector<FullTensor> fixedComponents(misc::max(_x.dimensions));
		
		for(size_t corePosition = _x.degree()-1; corePosition > 0; --corePosition) {
			// Prepare the single slates of the current component 
			for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
				fixedComponents[i](r1, r2) = _x.get_component(corePosition)(r1, i, r2);
			}
			
			// Reset the FullTensors
			for(const size_t i : _backwardUpdates[corePosition]) {
				contract(*_backwardStack[i + corePosition*numMeasurments], fixedComponents[_measurments.positions[i][corePosition]], false, *_backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
			}
		}
	}
	
	/**
	 * @brief: Calculates the difference between the current and measured values at the measured positions.
	 * @note: Abuses the stack infrastructe for its caluclation. Changes only the stacks at corePosition.
	 */
	inline std::vector<value_t> calculate_current_differences(
									const TTTensor& _x,
									const SinglePointMeasurmentSet& _measurments,
									const std::vector<std::vector<size_t>>& _forwardUpdates,
									const std::vector<std::vector<size_t>>& _backwardUpdates,
									FullTensor* const * const _forwardStack, 
									FullTensor* const * const _backwardStack, 
									const size_t _corePosition) {
		const size_t numMeasurments = _measurments.size();
		const Index r1, r2;
		
		std::vector<FullTensor> fixedComponents(_x.dimensions[_corePosition]);
		FullTensor currentValue({});
		
		// Prepare the single slates of the current component 
		for(size_t i = 0; i < _x.dimensions[_corePosition]; ++i) {
			fixedComponents[i](r1, r2) = _x.get_component(_corePosition)(r1, i, r2);
		}
		
		std::vector<value_t> currentDifferences(numMeasurments);
		
		// Look which side of the stack needs less calculations
		if(_forwardUpdates[_corePosition].size() < _backwardUpdates[_corePosition].size()) {
			for(const size_t i : _forwardUpdates[_corePosition]) {
				contract(*_forwardStack[i + _corePosition*numMeasurments] , *_forwardStack[i + (_corePosition-1)*numMeasurments], false, fixedComponents[_measurments.positions[i][_corePosition]], false, 1);
			}
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *_forwardStack[i + _corePosition*numMeasurments], false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
				currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
			}
		} else {
			for(const size_t i : _backwardUpdates[_corePosition]) {
				contract(*_backwardStack[i + _corePosition*numMeasurments], fixedComponents[_measurments.positions[i][_corePosition]], false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
			}
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *_forwardStack[i + (_corePosition-1)*numMeasurments], false, *_backwardStack[i + _corePosition*numMeasurments], false, 1);
				currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
			}
		}
		
		return currentDifferences;
	}
	
	inline std::vector<FullTensor> calculate_deltas(
									const TTTensor& _x, 
									const SinglePointMeasurmentSet& _measurments,
									const FullTensor* const * const _forwardStack, 
									const FullTensor* const * const _backwardStack,
									const std::vector<value_t>& _currentDifferences,
									const size_t _corePosition) {
		const size_t localLeftRank = _x.get_component(_corePosition).dimensions[0];
		const size_t localRightRank = _x.get_component(_corePosition).dimensions[2];
		const size_t numMeasurments = _measurments.size();
		
		std::vector<FullTensor> deltas(_x.dimensions[_corePosition], FullTensor({localLeftRank, localRightRank}));

		for(size_t i = 0; i < _x.dimensions[_corePosition]; ++i) {
			deltas[i].ensure_own_data();
		}
		
		for(size_t i = 0; i < numMeasurments; ++i) {
			// Interestingly writing a dyadic product on our own turns out to be faster than blas...
			const value_t* const leftPtr = _forwardStack[i + (_corePosition-1)*numMeasurments]->data.get();
			const value_t* const rightPtr = _backwardStack[i + (_corePosition+1)*numMeasurments]->data.get();
			value_t* const deltaPtr = deltas[_measurments.positions[i][_corePosition]].data.get();
			const value_t factor = _currentDifferences[i]*_forwardStack[i + (_corePosition-1)*numMeasurments]->factor*_backwardStack[i + (_corePosition+1)*numMeasurments]->factor;
			for(size_t k = 0; k < localLeftRank; ++k) {
				for(size_t j = 0; j < localRightRank; ++j) {
					deltaPtr[k*localRightRank+j] += factor * leftPtr[k] * rightPtr[j];
				}
			}
		}
		
		return deltas;
	}
	
	inline std::vector<value_t> calculate_slicewise_norm_Py(
												const std::vector<FullTensor>& _deltas,
												const SinglePointMeasurmentSet& _measurments,
												const std::vector<std::vector<size_t>>& _forwardUpdates,
												const std::vector<std::vector<size_t>>& _backwardUpdates,
												FullTensor* const * const _forwardStack, 
												FullTensor* const * const _backwardStack,
												const size_t _corePosition) {
		const size_t numMeasurments = _measurments.size();
		
		std::vector<value_t> PyPys(_deltas.size());
		
		FullTensor currentValue({});
		
		// Look which side of the stack needs less calculations
		if(_forwardUpdates[_corePosition].size() < _backwardUpdates[_corePosition].size()) {
			for(const size_t i : _forwardUpdates[_corePosition]) {
				contract(*_forwardStack[i + _corePosition*numMeasurments] , *_forwardStack[i + (_corePosition-1)*numMeasurments], false, _deltas[_measurments.positions[i][_corePosition]], false, 1);
			}
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *_forwardStack[i + _corePosition*numMeasurments], false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
				PyPys[_measurments.positions[i][_corePosition]] += misc::sqr(currentValue[0]);
			}
		} else {
			for(const size_t i : _backwardUpdates[_corePosition]) {
				contract(*_backwardStack[i + _corePosition*numMeasurments], _deltas[_measurments.positions[i][_corePosition]], false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
			}
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *_forwardStack[i + (_corePosition-1)*numMeasurments], false, *_backwardStack[i + _corePosition*numMeasurments], false, 1);
				PyPys[_measurments.positions[i][_corePosition]] += misc::sqr(currentValue[0]);
			}
		}
		
		return PyPys;
	}
	
	
	inline void update_forward_stack( FullTensor* const * const _forwardStack, const TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const std::vector<std::vector<size_t>>& _forwardUpdates, const size_t _corePosition) {
		const size_t numMeasurments = _measurments.size();
		const Index r1, r2;
		
		std::vector<FullTensor> fixedComponents(_x.dimensions[_corePosition]);
		
		// Prepare the single slates of the current component 
		for(size_t i = 0; i < _x.dimensions[_corePosition]; ++i) {
			fixedComponents[i](r1, r2) = _x.get_component(_corePosition)(r1, i, r2);
		}
		
		// Update the stack
		for(const size_t i : _forwardUpdates[_corePosition]) {
			contract(*_forwardStack[i + _corePosition*numMeasurments] , *_forwardStack[i + (_corePosition-1)*numMeasurments], false, fixedComponents[_measurments.positions[i][_corePosition]], false, 1);
		}
	}
	
	void ADFVariant::solve_with_current_ranks(TTTensor& _x, 
												size_t& _iteration,
												double& _residual,
												double& _lastResidual,
												const SinglePointMeasurmentSet& _measurments,
												const value_t _normMeasuredValues,
												FullTensor* const * const _forwardStack, 
												const std::vector<std::vector<size_t>>& _forwardUpdates, 
												FullTensor* const * const _backwardStack,
												const std::vector<std::vector<size_t>>& _backwardUpdates
												) const {
		const size_t numMeasurments = _measurments.size();
		const size_t degree = _x.degree();
		
		double resDec1 = 1.0, resDec2 = 1.0, resDec3 = 1.0;
		
		std::vector<value_t> currentDifferences(numMeasurments);
			
		for(; _iteration < maxInterations; ++_iteration) {
			// Move core back to position zero
			_x.move_core(0, true);
			
			rebuild_backward_stack(_backwardStack, _x, _measurments, _backwardUpdates);
			
			// Calculate the current residual
			currentDifferences = calculate_current_differences(_x, _measurments, _forwardUpdates, _backwardUpdates, _forwardStack, _backwardStack, 0);
			
			_lastResidual = _residual;
			_residual = 0;
			for(size_t i = 0; i < numMeasurments; ++i) {
				_residual += misc::sqr(currentDifferences[i]);
			}
			_residual = std::sqrt(_residual)/_normMeasuredValues;
			
			LOG(ADF, "Rank " << _x.ranks() << " Itr: " << _iteration << " Residual: " << std::scientific << _residual << " Rel. Residual change: " << (_lastResidual-_residual)/_lastResidual);
			
			// Check for termination criteria
			resDec3 = resDec2; resDec2 = resDec1;
			resDec1 = (_lastResidual-_residual)/_lastResidual;
			if(_residual < targetResidual || resDec1+resDec2+resDec3 < 3*minimalResidualDecrese) { break; }
			
			// Sweep from the first to the last component
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				if(corePosition > 0) { // For corePosition 0 this calculation is allready done in the calculation of the residual.
					currentDifferences = calculate_current_differences(_x, _measurments, _forwardUpdates, _backwardUpdates, _forwardStack, _backwardStack, corePosition);
				}
				
				const std::vector<FullTensor> deltas = calculate_deltas(_x, _measurments, _forwardStack, _backwardStack, currentDifferences, corePosition);
				
				const std::vector<value_t> PyPys = calculate_slicewise_norm_Py(deltas, _measurments, _forwardUpdates, _backwardUpdates, _forwardStack, _backwardStack, corePosition);
				
				// Update each slice seperately
				for(size_t j = 0; j < _x.dimensions[corePosition]; ++j) {
					// Calculate <P(y), P(X-B)> = ||deltaPlus||^2.
					const value_t PyR = misc::sqr(frob_norm(deltas[j]));
					
					// Update
					_x.component(corePosition)(r1, i1, r2) = _x.component(corePosition)(r1, i1, r2) + (PyR/PyPys[j])*Tensor::dirac({_x.dimensions[corePosition]}, {j})(i1)*deltas[j](r1, r2);
				}
				
				// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
				if(corePosition+1 < degree) {
					_x.move_core(corePosition+1, true);
					
					update_forward_stack( _forwardStack, _x, _measurments, _forwardUpdates, corePosition);
				}
			}
		}
	}
	
	
	double ADFVariant::solve(TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const std::vector<size_t>& _maxRanks) const {
		REQUIRE(_x.is_valid_tt(), "_x must be a valid TT-Tensor.");
		
		const size_t degree = _x.degree();
		const size_t numMeasurments = _measurments.size();
		
		REQUIRE(numMeasurments > 0, "Need at very least one measurment.");
		REQUIRE(_measurments.degree() == degree, "Measurment degree must coincide with _x degree.");
		
		const value_t normMeasuredValues = calculate_norm_of_measured_values(_measurments.measuredValues);
		
		// Vectors containing for each core position, the IDs of all unqiue stack entries that need updating.
		std::vector<std::vector<size_t>> forwardUpdates(degree);
		std::vector<std::vector<size_t>> backwardUpdates(degree);
		
		// Array that will contain all the unqiue stack entries.
		std::unique_ptr<FullTensor[]> forwardStackSaveSlots;
		std::unique_ptr<FullTensor[]> backwardStackSaveSlots;
		
		// Arrays mapping the stack entry to the corresponding unique entry in stackSaveSlots (NOTE that both allow position -1 and degree).
		std::unique_ptr<FullTensor*[]> forwardStackMem( new FullTensor*[numMeasurments*(degree+2)]);
		FullTensor* const * const forwardStack = forwardStackMem.get()+numMeasurments;
		
		std::unique_ptr<FullTensor*[]> backwardStackMem( new FullTensor*[numMeasurments*(degree+2)]);
		FullTensor* const * const backwardStack = backwardStackMem.get()+numMeasurments;
		
		// We need _x to be canonicalized in the sense that there is no edge with more than maximal rank.
		_x.cannonicalize_left();
		
		construct_stacks(forwardStackSaveSlots, forwardUpdates, forwardStackMem, _x, _measurments, true);
		construct_stacks(backwardStackSaveSlots, backwardUpdates, backwardStackMem, _x, _measurments, false);
		
		size_t iteration = 0;
		double residual = std::numeric_limits<double>::max(), lastResidual;
		
		// One inital run
		solve_with_current_ranks(_x, iteration, residual, lastResidual, _measurments, normMeasuredValues, forwardStack, forwardUpdates, backwardStack, backwardUpdates);
		
		// If we follow a rank increasing strategie, increase the ransk until we reach the targetResidual or the maxRanks.
		while(residual > targetResidual && _x.ranks() != _maxRanks) {
			
			// Increase the ranks
			_x = _x+((1e-6*frob_norm(_x)/misc::fp_product(_x.dimensions))*TTTensor::ones(_x.dimensions));
			
			// Resize stacks
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				for(const size_t i : forwardUpdates[corePosition]) {
					forwardStack[i + corePosition*numMeasurments]->reset({_x.rank(corePosition)}, DONT_SET_ZERO());
				}
				for(const size_t i : backwardUpdates[corePosition]) {
					backwardStack[i + corePosition*numMeasurments]->reset({_x.rank(corePosition - 1)}, DONT_SET_ZERO());
				}
			}
			
			solve_with_current_ranks(_x, iteration, residual, lastResidual, _measurments, normMeasuredValues, forwardStack, forwardUpdates, backwardStack, backwardUpdates);
			
			// Ensure maximal ranks are not exceeded (may happen if non uniform).
			_x.round(_maxRanks); // TODO do not round edges that do not need it (i.e. do not even do the SVD!)
		}
		return residual;
	}
}
