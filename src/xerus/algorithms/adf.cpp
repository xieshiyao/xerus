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
	
	template<class MeasurmentSet>
	double ADFVariant::InternalSolver<MeasurmentSet>::calculate_norm_of_measured_values(const MeasurmentSet& _measurments) {
		value_t normMeasuredValues = 0;
		for(const value_t measurement : _measurments.measuredValues) {
			normMeasuredValues += misc::sqr(measurement);
		}
		return std::sqrt(normMeasuredValues);
	}
	
	template<class MeasurmentSet>
	class MeasurmentComparator {
		const bool forward;
		const MeasurmentSet& measurments;
		const size_t degree;
	public:
		MeasurmentComparator(const MeasurmentSet& _measurments, const bool _forward) : forward(_forward), measurments(_measurments), degree(_measurments.degree()) {}
		bool operator()(const size_t _a, const size_t _b) const;
	};
	
	template<>
	bool MeasurmentComparator<SinglePointMeasurmentSet>::operator()(const size_t _a, const size_t _b) const {
		if(forward) {
			for (size_t j = 0; j < degree; ++j) {
				if (measurments.positions[_a][j] < measurments.positions[_b][j]) return true;
				if (measurments.positions[_a][j] > measurments.positions[_b][j]) return false;
			}
		} else {
			for (size_t j = degree; j > 0; --j) {
				if (measurments.positions[_a][j-1] < measurments.positions[_b][j-1]) return true;
				if (measurments.positions[_a][j-1] > measurments.positions[_b][j-1]) return false;
			}
		}
		LOG(fatal, "Measurments must not appear twice.");
		return false;
	}
	
	template<>
	bool MeasurmentComparator<RankOneMeasurmentSet>::operator()(const size_t _a, const size_t _b) const {
		if(forward) {
			for (size_t j = 0; j < degree; ++j) {
				for(size_t k = 0; k < measurments.positions[_a][j].size; ++k) {
					if (measurments.positions[_a][j][k] < measurments.positions[_b][j][k]) return true;
					if (measurments.positions[_a][j][k] > measurments.positions[_b][j][k]) return false;
				}
			}
		} else {
			for (size_t j = degree; j > 0; --j) {
				for(size_t k = 0; k < measurments.positions[_a][j].size; ++k) {
					if (measurments.positions[_a][j-1][k] < measurments.positions[_b][j-1][k]) return true;
					if (measurments.positions[_a][j-1][k] > measurments.positions[_b][j-1][k]) return false;
				}
			}
		}
		LOG(fatal, "Measurments must not appear twice.");
		return false;
	}
	
	
	template<class MeasurmentSet>
	void ADFVariant::InternalSolver<MeasurmentSet>::construct_stacks(std::unique_ptr<FullTensor[]>& _stackSaveSlot, std::vector<std::vector<size_t>>& _updates, std::unique_ptr<FullTensor*[]>& _stackMem, const bool _forward) {
		// Temporary map. For each stack entrie (i.e. measurement number + corePosition) gives a measurement number of a stack entrie (at same corePosition) that shall have an equal value (or its own number otherwise). 
		std::vector<size_t> calculationMap(degree*numMeasurments);
		
		// Count how many FullTensors we need for the stacks
		size_t numUniqueStackEntries = 0;
		
		// Create a reordering map
		std::vector<size_t> reorderedMeasurments(numMeasurments);
		std::iota(reorderedMeasurments.begin(), reorderedMeasurments.end(), 0);
		
		std::sort(reorderedMeasurments.begin(), reorderedMeasurments.end(), MeasurmentComparator<MeasurmentSet>(measurments, _forward) );
		
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
			for( ; position < degree && approx_equal(measurments.positions[realId][corePosition], measurments.positions[realPreviousId][corePosition]); ++position, corePosition = _forward ? position : degree-1-position) {
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
					_stackSaveSlot[usedSlots].reset({x.rank(corePosition - (_forward ? 0 : 1))}, DONT_SET_ZERO());
					_stackMem[i + (corePosition+1)*numMeasurments] = &_stackSaveSlot[usedSlots++];
				} else {
					_stackMem[i + (corePosition+1)*numMeasurments] = _stackMem[calculationMap[i + corePosition*numMeasurments] + (corePosition+1)*numMeasurments];
				}
			}
		}
		
		REQUIRE(usedSlots == numUniqueStackEntries, "Internal Error.");
		LOG(ADF, "We have " << numUniqueStackEntries << " unique stack entries. There are " << numMeasurments*degree+1 << " virtual stack entries.");
	}
	
	template<class MeasurmentSet>
	inline void rebuild_backward_stack(FullTensor* const * const _backwardStack, const TTTensor& _x, const MeasurmentSet& _measurments, const std::vector<std::vector<size_t>>& _backwardUpdates);
	
	template<>
	inline void rebuild_backward_stack(FullTensor* const * const _backwardStack, const TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const std::vector<std::vector<size_t>>& _backwardUpdates) {
		const size_t numMeasurments = _measurments.size();
		
		const Index r1, r2; //TODO
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
	
	template<>
	inline void rebuild_backward_stack(FullTensor* const * const _backwardStack, const TTTensor& _x, const RankOneMeasurmentSet& _measurments, const std::vector<std::vector<size_t>>& _backwardUpdates) {
		const size_t numMeasurments = _measurments.size();
		
		const Index r1, r2, i1; // TODO
		
		FullTensor reshuffledComponent;
		
		for(size_t corePosition = _x.degree()-1; corePosition > 0; --corePosition) {
			reshuffledComponent(r1, r2, i1) =  _x.get_component(corePosition)(r1, i1, r2);
			
			FullTensor mixedComponent({reshuffledComponent.dimensions[0], reshuffledComponent.dimensions[1]});
			
			// Reset the FullTensors
			for(const size_t i : _backwardUpdates[corePosition]) {
				contract(mixedComponent, reshuffledComponent, false, _measurments.positions[i][corePosition], false, 1);
				contract(*_backwardStack[i + corePosition*numMeasurments], mixedComponent, false, *_backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
			}
		}
	}
	
	/**
	 * @brief: Calculates the difference between the current and measured values at the measured positions.
	 * @note: Abuses the stack infrastructe for its caluclation. Changes only the stacks at corePosition.
	 */
	template<class MeasurmentSet>
	inline std::vector<value_t> calculate_current_differences(
									const TTTensor& _x,
									const MeasurmentSet& _measurments,
									const std::vector<std::vector<size_t>>& _forwardUpdates,
									const std::vector<std::vector<size_t>>& _backwardUpdates,
									FullTensor* const * const _forwardStack, 
									FullTensor* const * const _backwardStack, 
									const size_t _corePosition);
	
	template<>
	inline std::vector<value_t> calculate_current_differences<SinglePointMeasurmentSet>(
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
	
	template<>
	inline std::vector<value_t> calculate_current_differences<RankOneMeasurmentSet>(
									const TTTensor& _x,
									const RankOneMeasurmentSet& _measurments,
									const std::vector<std::vector<size_t>>& _forwardUpdates,
									const std::vector<std::vector<size_t>>& _backwardUpdates,
									FullTensor* const * const _forwardStack, 
									FullTensor* const * const _backwardStack, 
									const size_t _corePosition) {
		const size_t numMeasurments = _measurments.size();
		const Index r1, r2, i1;
		
		FullTensor currentValue({});
		std::vector<value_t> currentDifferences(numMeasurments);
		
		
		FullTensor reshuffledComponent;
		reshuffledComponent(r1, r2, i1) =  _x.get_component(_corePosition)(r1, i1, r2);
		FullTensor mixedComponent({reshuffledComponent.dimensions[0], reshuffledComponent.dimensions[1]});
		
		// Look which side of the stack needs less calculations
		if(_forwardUpdates[_corePosition].size() < _backwardUpdates[_corePosition].size()) {
			for(const size_t i : _forwardUpdates[_corePosition]) {
				contract(mixedComponent, reshuffledComponent, false, _measurments.positions[i][_corePosition], false, 1);
				contract(*_forwardStack[i + _corePosition*numMeasurments] , *_forwardStack[i + (_corePosition-1)*numMeasurments], false, mixedComponent, false, 1);
			}
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *_forwardStack[i + _corePosition*numMeasurments], false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
				currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
			}
		} else {
			for(const size_t i : _backwardUpdates[_corePosition]) {
				contract(mixedComponent, reshuffledComponent, false, _measurments.positions[i][_corePosition], false, 1);
				contract(*_backwardStack[i + _corePosition*numMeasurments], mixedComponent, false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
			}
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *_forwardStack[i + (_corePosition-1)*numMeasurments], false, *_backwardStack[i + _corePosition*numMeasurments], false, 1);
				currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
			}
		}
		
		return currentDifferences;
	}
	
	template<class MeasurmentSet>
	inline std::vector<FullTensor> calculate_deltas(
									const TTTensor& _x, 
									const MeasurmentSet& _measurments,
									const FullTensor* const * const _forwardStack, 
									const FullTensor* const * const _backwardStack,
									const std::vector<value_t>& _currentDifferences,
									const size_t _corePosition);
	
	template<>
	inline std::vector<FullTensor> calculate_deltas<SinglePointMeasurmentSet>(
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
	
	template<>
	inline std::vector<FullTensor> calculate_deltas<RankOneMeasurmentSet>(
									const TTTensor& _x, 
									const RankOneMeasurmentSet& _measurments,
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
			for(size_t n = 0; n < _x.dimensions[_corePosition]; ++n) {
				value_t* const deltaPtr = deltas[n].data.get();
				const value_t factor = _measurments.positions[i][_corePosition][n]*_currentDifferences[i]*_forwardStack[i + (_corePosition-1)*numMeasurments]->factor*_backwardStack[i + (_corePosition+1)*numMeasurments]->factor;
				for(size_t k = 0; k < localLeftRank; ++k) {
					for(size_t j = 0; j < localRightRank; ++j) {
						deltaPtr[k*localRightRank+j] += factor * leftPtr[k] * rightPtr[j];
					}
				}
			}
		}
		return deltas;
	}
	
	template<class MeasurmentSet>
	inline std::vector<value_t> calculate_slicewise_norm_Py(
												const std::vector<FullTensor>& _deltas,
												const MeasurmentSet& _measurments,
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
	
	
	template<class MeasurmentSet>
	inline void update_forward_stack( FullTensor* const * const _forwardStack, const TTTensor& _x, const MeasurmentSet& _measurments, const std::vector<std::vector<size_t>>& _forwardUpdates, const size_t _corePosition) {
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
	
	template<class MeasurmentSet>
	void ADFVariant::InternalSolver<MeasurmentSet>::solve_with_current_ranks() {
		double resDec1 = 1.0, resDec2 = 1.0, resDec3 = 1.0;
		
		std::vector<value_t> currentDifferences(numMeasurments);
			
		for(; iteration < maxInterations; ++iteration) {
			// Move core back to position zero
			x.move_core(0, true);
			
			rebuild_backward_stack(backwardStack, x, measurments, backwardUpdates);
			
			// Calculate the current residual
			currentDifferences = calculate_current_differences(x, measurments, forwardUpdates, backwardUpdates, forwardStack, backwardStack, 0);
			
			lastResidual = residual;
			residual = 0;
			for(size_t i = 0; i < numMeasurments; ++i) {
				residual += misc::sqr(currentDifferences[i]);
			}
			residual = std::sqrt(residual)/normMeasuredValues;
			
			perfData.add(iteration, residual, x.ranks(), 0);
// 			LOG(ADF, "Rank " << x.ranks() << " Itr: " << iteration << " Residual: " << std::scientific << residual << " Rel. Residual change: " << (lastResidual-_residual)/_lastResidual);
			
			// Check for termination criteria
			resDec3 = resDec2; resDec2 = resDec1;
			resDec1 = (lastResidual-residual)/lastResidual;
			if(residual < targetResidual || resDec1+resDec2+resDec3 < 3*minimalResidualDecrese) { break; }
			
			// Sweep from the first to the last component
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				if(corePosition > 0) { // For corePosition 0 this calculation is allready done in the calculation of the residual.
					currentDifferences = calculate_current_differences(x, measurments, forwardUpdates, backwardUpdates, forwardStack, backwardStack, corePosition);
				}
				
				const std::vector<FullTensor> deltas = calculate_deltas(x, measurments, forwardStack, backwardStack, currentDifferences, corePosition);
				
				const std::vector<value_t> PyPys = calculate_slicewise_norm_Py(deltas, measurments, forwardUpdates, backwardUpdates, forwardStack, backwardStack, corePosition);
				
				// Update each slice seperately
				for(size_t j = 0; j < x.dimensions[corePosition]; ++j) {
					// Calculate <P(y), P(X-B)> = ||deltaPlus||^2.
					const value_t PyR = misc::sqr(frob_norm(deltas[j]));
					
					// Update
					x.component(corePosition)(r1, i1, r2) = x.component(corePosition)(r1, i1, r2) + (PyR/PyPys[j])*Tensor::dirac({x.dimensions[corePosition]}, {j})(i1)*deltas[j](r1, r2);
				}
				
				// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
				if(corePosition+1 < degree) {
					x.move_core(corePosition+1, true);
					
					update_forward_stack( forwardStack, x, measurments, forwardUpdates, corePosition);
				}
			}
		}
	}
	
	template<class MeasurmentSet>
	double ADFVariant::InternalSolver<MeasurmentSet>::solve() {
		perfData.start();
		
		// We need x to be canonicalized in the sense that there is no edge with more than maximal rank.
		x.cannonicalize_left();
		
		construct_stacks(forwardStackSaveSlots, forwardUpdates, forwardStackMem, true);
		construct_stacks(backwardStackSaveSlots, backwardUpdates, backwardStackMem, false);
		
		// One inital run
		solve_with_current_ranks();
		
		// If we follow a rank increasing strategie, increase the ransk until we reach the targetResidual or the maxRanks.
		while(residual > targetResidual && x.ranks() != maxRanks) {
			
			// Increase the ranks
			x = x+((1e-6*frob_norm(x)/misc::fp_product(x.dimensions))*TTTensor::ones(x.dimensions));
			
			// Resize stacks
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				for(const size_t i : forwardUpdates[corePosition]) {
					forwardStack[i + corePosition*numMeasurments]->reset({x.rank(corePosition)}, DONT_SET_ZERO());
				}
				for(const size_t i : backwardUpdates[corePosition]) {
					backwardStack[i + corePosition*numMeasurments]->reset({x.rank(corePosition - 1)}, DONT_SET_ZERO());
				}
			}
			
			solve_with_current_ranks();
			
			// Ensure maximal ranks are not exceeded (may happen if non uniform).
			x.round(maxRanks); // TODO do not round edges that do not need it (i.e. do not even do the SVD!)
		}
		return residual;
	}
	

	double ADFVariant::operator()(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments, PerformanceData& _perfData) const {
	   InternalSolver<SinglePointMeasurmentSet> solver(_x, _x.ranks(), SinglePointMeasurmentSet(_measurments), maxInterations, targetResidual, minimalResidualDecrese, _perfData);
	   return solver.solve();
	}
	
	double ADFVariant::operator()(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData) const {
	   InternalSolver<SinglePointMeasurmentSet> solver(_x, _maxRanks, SinglePointMeasurmentSet(_measurments), maxInterations, targetResidual, minimalResidualDecrese, _perfData);
	   return solver.solve();
	}
	
	double ADFVariant::operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, PerformanceData& _perfData) const {
// 	  InternalSolver<RankOneMeasurmentSet> solver(_x, _x.ranks(), _measurments, maxInterations, targetResidual, minimalResidualDecrese, _perfData);
// 	  return solver.solve();
		return 0.0;
	}
	
	double ADFVariant::operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData) const {
// 	   InternalSolver<RankOneMeasurmentSet> solver(_x, _maxRanks, _measurments, maxInterations, targetResidual, minimalResidualDecrese, _perfData);
// 	   return solver.solve();
		return 0.0;
	}
}
