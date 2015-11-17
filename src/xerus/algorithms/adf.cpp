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
				REQUIRE(measurments.positions[_a][j].size == measurments.positions[_b][j].size, "IE: " << _a << " " << _b << " | " << measurments.positions[_a][j].size << " vs " << measurments.positions[_b][j].size);
				for(size_t k = 0; k < measurments.positions[_a][j].size; ++k) {
					if (measurments.positions[_a][j][k] < measurments.positions[_b][j][k]) { return true; }
					if (measurments.positions[_a][j][k] > measurments.positions[_b][j][k]) { return false; }
				}
			}
		} else {
			for (size_t j = degree; j > 0; --j) {
				REQUIRE(measurments.positions[_a][j-1].size == measurments.positions[_b][j-1].size, "IE: " << _a << " " << _b << " | " << measurments.positions[_a][j-1].size << " vs " << measurments.positions[_b][j-1].size);
				for(size_t k = 0; k < measurments.positions[_a][j-1].size; ++k) {
					if (measurments.positions[_a][j-1][k] < measurments.positions[_b][j-1][k]) { return true; }
					if (measurments.positions[_a][j-1][k] > measurments.positions[_b][j-1][k]) { return false; }
				}
			}
		}
		LOG(fatal, "Measurments must not appear twice.");
		return false;
	}
	
	
	template<class MeasurmentSet>
	void ADFVariant::InternalSolver<MeasurmentSet>::construct_stacks(std::unique_ptr<FullTensor[]>& _stackSaveSlot, std::vector<std::vector<size_t>>& _updates, std::unique_ptr<FullTensor*[]>& _stackMem, const bool _forward) {
		// Temporary map. For each stack entry (i.e. measurement number + corePosition) gives a measurement number of a stack entry (at same corePosition) that shall have an equal value (or its own number otherwise). 
		std::vector<size_t> calculationMap(degree*numMeasurments);
		
		// Count how many FullTensors we need for the stacks
		size_t numUniqueStackEntries = 0;
		
		// Create a reordering map
		std::vector<size_t> reorderedMeasurments(numMeasurments);
		std::iota(reorderedMeasurments.begin(), reorderedMeasurments.end(), 0);
		
		std::sort(reorderedMeasurments.begin(), reorderedMeasurments.end(), MeasurmentComparator<MeasurmentSet>(measurments, _forward) );
		
		// Create the entries for the first measurement (these are allways unqiue).
		for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
			const size_t realId = reorderedMeasurments[0];
			calculationMap[realId + corePosition*numMeasurments] = realId;
			++numUniqueStackEntries;
		}
		
		// Create the calculation map
		for(size_t i = 1; i < numMeasurments; ++i) {
			const size_t realId = reorderedMeasurments[i];
			const size_t realPreviousId = reorderedMeasurments[i-1];
			
			size_t position = 0;
			size_t corePosition = _forward ? position : degree-1-position;
			for( ; 
				position < degree && approx_equal(measurments.positions[realId][corePosition], measurments.positions[realPreviousId][corePosition]);
				++position, corePosition = _forward ? position : degree-1-position) 
			{
				if( realPreviousId < realId ) {
					calculationMap[realId + corePosition*numMeasurments] = calculationMap[realPreviousId + corePosition*numMeasurments];
				} else if(realPreviousId == calculationMap[realPreviousId + corePosition*numMeasurments]) {
					calculationMap[realPreviousId + corePosition*numMeasurments] = realId;
					calculationMap[realId + corePosition*numMeasurments] = realId;
				} else if( realId < calculationMap[realPreviousId + corePosition*numMeasurments]) {
					const size_t nextOther = calculationMap[realPreviousId + corePosition*numMeasurments];
					REQUIRE(calculationMap[nextOther + corePosition*numMeasurments] == nextOther, "IE");
					calculationMap[realPreviousId + corePosition*numMeasurments] = realId;
					calculationMap[nextOther + corePosition*numMeasurments] = realId;
					calculationMap[realId + corePosition*numMeasurments] = realId;
				} else {
					calculationMap[realId + corePosition*numMeasurments] = calculationMap[realPreviousId + corePosition*numMeasurments];
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
			_stackMem[i + 0*numMeasurments] = &_stackSaveSlot[0];
			_stackMem[i + (degree+1)*numMeasurments] = &_stackSaveSlot[0];
		}
		// TODO use stack insteasd of statMem
		for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(calculationMap[i + corePosition*numMeasurments] == i) {
					_updates[corePosition].emplace_back(i);
					_stackSaveSlot[usedSlots].reset({x.rank(corePosition - (_forward ? 0 : 1))}, DONT_SET_ZERO()); // TODO brocken -1 and degree
					_stackMem[i + (corePosition+1)*numMeasurments] = &_stackSaveSlot[usedSlots];
					usedSlots++;
				} else {
					_stackMem[i + (corePosition+1)*numMeasurments] = _stackMem[calculationMap[i + corePosition*numMeasurments] + (corePosition+1)*numMeasurments];
				}
			}
		}
		
		REQUIRE(usedSlots == numUniqueStackEntries, "Internal Error.");
		perfData << "We have " << numUniqueStackEntries << " unique stack entries. There are " << numMeasurments*degree+1 << " virtual stack entries.";
	}
	
	template<class MeasurmentSet>
	std::vector<FullTensor> ADFVariant::InternalSolver<MeasurmentSet>::get_fixed_components(const Tensor& _component) {
		std::vector<FullTensor> fixedComponents(_component.dimensions[1]);
		
		for(size_t i = 0; i < _component.dimensions[1]; ++i) {
			fixedComponents[i](r1, r2) = _component(r1, i, r2);
		}
		
		return fixedComponents;
	}
	
	template<>
	void ADFVariant::InternalSolver<SinglePointMeasurmentSet>::update_backward_stack(const size_t _corePosition, const Tensor& _currentComponent) {
		std::vector<FullTensor> fixedComponents = get_fixed_components(_currentComponent);
		
		// Reset the FullTensors
		for(const size_t i : backwardUpdates[_corePosition]) {
			contract(*backwardStack[i + _corePosition*numMeasurments], fixedComponents[measurments.positions[i][_corePosition]], false, *backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
		}
	}
	
	template<>
	void ADFVariant::InternalSolver<RankOneMeasurmentSet>::update_backward_stack(const size_t _corePosition, const Tensor& _currentComponent) {
		FullTensor reshuffledComponent;
		reshuffledComponent(i1, r1, r2) =  _currentComponent(r1, i1, r2);

		FullTensor mixedComponent({reshuffledComponent.dimensions[1], reshuffledComponent.dimensions[2]});
		
		// Reset the FullTensors
		for(const size_t i : backwardUpdates[_corePosition]) {
			contract(mixedComponent, measurments.positions[i][_corePosition], false, reshuffledComponent, false, 1);
			contract(*backwardStack[i + _corePosition*numMeasurments], mixedComponent, false, *backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
		}
	}
	
	
	template<>
	void ADFVariant::InternalSolver<SinglePointMeasurmentSet>::update_forward_stack( const size_t _corePosition, const Tensor& _currentComponent ) {
		std::vector<FullTensor> fixedComponents = get_fixed_components(_currentComponent);		
		
		// Update the stack
		for(const size_t i : forwardUpdates[_corePosition]) {
			contract(*forwardStack[i + _corePosition*numMeasurments] , *forwardStack[i + (_corePosition-1)*numMeasurments], false, fixedComponents[measurments.positions[i][_corePosition]], false, 1);
		}
	}
	
	template<>
	void ADFVariant::InternalSolver<RankOneMeasurmentSet>::update_forward_stack( const size_t _corePosition, const Tensor& _currentComponent ) {
		FullTensor reshuffledComponent;
		reshuffledComponent(i1, r1, r2) =  _currentComponent(r1, i1, r2);

		FullTensor mixedComponent({reshuffledComponent.dimensions[1], reshuffledComponent.dimensions[2]});
		
		// Reset the FullTensors
		for(const size_t i : backwardUpdates[_corePosition]) {
			contract(mixedComponent, measurments.positions[i][_corePosition], false, reshuffledComponent, false, 1);
			contract(*forwardStack[i + _corePosition*numMeasurments] , *forwardStack[i + (_corePosition-1)*numMeasurments], false, mixedComponent, false, 1);
		}
	}
	
	
	template<class MeasurmentSet>
	void ADFVariant::InternalSolver<MeasurmentSet>::calculate_residual( const size_t _corePosition ) {
		FullTensor currentValue({});
		
		// Look which side of the stack needs less calculations
		if(forwardUpdates[_corePosition].size() < backwardUpdates[_corePosition].size()) {
			update_forward_stack(_corePosition, x.get_component(_corePosition));
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *forwardStack[i + _corePosition*numMeasurments], false, *backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
				residual[i] = (measurments.measuredValues[i]-currentValue[0]);
			}
		} else {
			update_backward_stack(_corePosition, x.get_component(_corePosition));
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *forwardStack[i + (_corePosition-1)*numMeasurments], false, *backwardStack[i + _corePosition*numMeasurments], false, 1);
				residual[i] = (measurments.measuredValues[i]-currentValue[0]);
			}
		}
	}
	
	template<>
	void ADFVariant::InternalSolver<SinglePointMeasurmentSet>::calculate_projected_gradient( const size_t _corePosition) {
		const size_t localLeftRank = x.get_component(_corePosition).dimensions[0];
		const size_t localRightRank = x.get_component(_corePosition).dimensions[2];
		
		projectedGradientComponent.reset({x.dimensions[_corePosition], localLeftRank, localRightRank});
		
		for(size_t i = 0; i < numMeasurments; ++i) {
			// Interestingly writing a dyadic product on our own turns out to be faster than blas...
			const value_t* const leftPtr = forwardStack[i + (_corePosition-1)*numMeasurments]->data_pointer();
			const value_t* const rightPtr = backwardStack[i + (_corePosition+1)*numMeasurments]->data_pointer();
			value_t* const deltaPtr = projectedGradientComponent.data_pointer()+measurments.positions[i][_corePosition]*localLeftRank*localRightRank;
			REQUIRE(!forwardStack[i + (_corePosition-1)*numMeasurments]->has_factor() && !backwardStack[i + (_corePosition+1)*numMeasurments]->has_factor(), "IE");
			const value_t factor = residual[i];
			for(size_t k = 0; k < localLeftRank; ++k) {
				for(size_t j = 0; j < localRightRank; ++j) {
					deltaPtr[k*localRightRank+j] += factor * leftPtr[k] * rightPtr[j];
				}
			}
		}
		
		projectedGradientComponent(r1, i1, r2) = projectedGradientComponent(i1, r1, r2);
	}
	
	template<>
	void ADFVariant::InternalSolver<RankOneMeasurmentSet>::calculate_projected_gradient( const size_t _corePosition) {
		const size_t localLeftRank = x.get_component(_corePosition).dimensions[0];
		const size_t localRightRank = x.get_component(_corePosition).dimensions[2];
		
		projectedGradientComponent.reset({x.dimensions[_corePosition], localLeftRank, localRightRank});
		
		for(size_t i = 0; i < numMeasurments; ++i) {
			// Interestingly writing a dyadic product on our own turns out to be faster than blas...
			const value_t* const leftPtr = forwardStack[i + (_corePosition-1)*numMeasurments]->data_pointer();
			const value_t* const rightPtr = backwardStack[i + (_corePosition+1)*numMeasurments]->data_pointer();
			for(size_t n = 0; n < x.dimensions[_corePosition]; ++n) {
				value_t* const deltaPtr = projectedGradientComponent.data_pointer() + n*localLeftRank*localRightRank;
				REQUIRE(!forwardStack[i + (_corePosition-1)*numMeasurments]->has_factor() && !backwardStack[i + (_corePosition+1)*numMeasurments]->has_factor(), "IE");
				const value_t factor = measurments.positions[i][_corePosition][n]*residual[i];
				for(size_t k = 0; k < localLeftRank; ++k) {
					for(size_t j = 0; j < localRightRank; ++j) {
						deltaPtr[k*localRightRank+j] += factor * leftPtr[k] * rightPtr[j];
					}
				}
			}
		}
		
		projectedGradientComponent(r1, i1, r2) = projectedGradientComponent(i1, r1, r2);
	}
	
	template<class MeasurmentSet>
	inline size_t position_or_zero(const MeasurmentSet& _measurments, const size_t _meas, const size_t _corePosition);
	
	template<>
	inline size_t position_or_zero<SinglePointMeasurmentSet>(const SinglePointMeasurmentSet& _measurments, const size_t _meas, const size_t _corePosition) {
		return _measurments.positions[_meas][_corePosition];
	}
	
	template<>
	inline size_t position_or_zero<RankOneMeasurmentSet>(const RankOneMeasurmentSet& _measurments, const size_t _meas, const size_t _corePosition) {
		return 0;
	}
	
	template<class MeasurmentSet>
	std::vector<value_t> ADFVariant::InternalSolver<MeasurmentSet>::calculate_slicewise_norm_A_projGrad( const size_t _corePosition) {
		std::vector<value_t> normAProjGrad(x.dimensions[_corePosition], 0.0);
		
		FullTensor currentValue({});
		
		// Look which side of the stack needs less calculations
		if(forwardUpdates[_corePosition].size() < backwardUpdates[_corePosition].size()) {
			update_forward_stack(_corePosition, projectedGradientComponent);
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *forwardStack[i + _corePosition*numMeasurments], false, *backwardStack[i + (_corePosition+1)*numMeasurments], false, 1);
				normAProjGrad[position_or_zero(measurments, i, _corePosition)] += misc::sqr(currentValue[0]);
			}
		} else {
			update_backward_stack(_corePosition, projectedGradientComponent);
			
			for(size_t i = 0; i < numMeasurments; ++i) {
				contract(currentValue, *forwardStack[i + (_corePosition-1)*numMeasurments], false, *backwardStack[i + _corePosition*numMeasurments], false, 1);
				normAProjGrad[position_or_zero(measurments, i, _corePosition)] += misc::sqr(currentValue[0]);
			}
		}
		
		return normAProjGrad;
	}
	
	
// 	template<>
// 	void ADFVariant::InternalSolver<SinglePointMeasurmentSet>::update_x(const std::vector<value_t>& _normAProjGrad, const size_t _corePosition) {
// 		for(size_t j = 0; j < x.dimensions[_corePosition]; ++j) {
// 			FullTensor localDelta;
// 			localDelta(r1, r2) = projectedGradientComponent(r1, j, r2);
// 			const value_t PyR = misc::sqr(frob_norm(localDelta));
// 			
// 			// Update
// 			x.component(_corePosition)(r1, i1, r2) = x.component(_corePosition)(r1, i1, r2) + (PyR/_normAProjGrad[j])*Tensor::dirac({x.dimensions[_corePosition]}, {j})(i1)*localDelta(r1, r2);
// 		}
// 	}
	
	template<>
	void ADFVariant::InternalSolver<SinglePointMeasurmentSet>::update_x(const std::vector<value_t>& _normAProjGrad, const size_t _corePosition) {
		const value_t PyR = misc::sqr(frob_norm(projectedGradientComponent));
		
		// Update
		x.component(_corePosition)(r1, i1, r2) = x.component(_corePosition)(r1, i1, r2) + (PyR/misc::sum(_normAProjGrad))*projectedGradientComponent(r1, i1, r2);
	}
	
	template<>
	void ADFVariant::InternalSolver<RankOneMeasurmentSet>::update_x(const std::vector<value_t>& _normAProjGrad, const size_t _corePosition) {
		const value_t PyR = misc::sqr(frob_norm(projectedGradientComponent));
		
		// Update
		x.component(_corePosition)(r1, i1, r2) = x.component(_corePosition)(r1, i1, r2) + (PyR/misc::sum(_normAProjGrad))*projectedGradientComponent(r1, i1, r2);
	}
	
	template<class MeasurmentSet>
	void ADFVariant::InternalSolver<MeasurmentSet>::solve_with_current_ranks() {
		double resDec1 = 1.0, resDec2 = 1.0, resDec3 = 1.0;
			
		for(; maxIterations == 0 || iteration < maxIterations; ++iteration) {
			
			// Move core back to position zero
			x.move_core(0, true);
			
			// Rebuild backwardStack
			for(size_t corePosition = x.degree()-1; corePosition > 0; --corePosition) {
				update_backward_stack(corePosition, x.get_component(corePosition));
			}
			
			calculate_residual(0);
			
			lastResidualNorm = residualNorm;
			residualNorm = 0;
			for(size_t i = 0; i < numMeasurments; ++i) {
				residualNorm += misc::sqr(residual[i]);
			}
			residualNorm = std::sqrt(residualNorm)/normMeasuredValues;
			
			perfData.add(iteration, residualNorm, x.ranks(), 0);
			
			// Check for termination criteria
			resDec3 = resDec2; resDec2 = resDec1;
			resDec1 = (lastResidualNorm-residualNorm)/lastResidualNorm;
			if(residualNorm < targetResidualNorm || resDec1+resDec2+resDec3 < 3*minimalResidualNormDecrese) { break; }
			
			// Sweep from the first to the last component
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				if(corePosition > 0) { // For corePosition 0 this calculation is allready done in the calculation of the residual.
					calculate_residual(corePosition);
				}
				
				calculate_projected_gradient(corePosition);
				
				const std::vector<value_t> normAProjGrad = calculate_slicewise_norm_A_projGrad(corePosition);
				
				update_x(normAProjGrad, corePosition);
				
				// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
				if(corePosition+1 < degree) {
					x.move_core(corePosition+1, true);
					
					update_forward_stack(corePosition, x.get_component(corePosition));
				}
			}
		}
	}
	
	template<class MeasurmentSet>
	double ADFVariant::InternalSolver<MeasurmentSet>::solve() {
		perfData.start();
		
		// We need x to be canonicalized in the sense that there is no edge with more than maximal rank (prior to stack creation).
		x.cannonicalize_left();
		
		construct_stacks(forwardStackSaveSlots, forwardUpdates, forwardStackMem, true);
		construct_stacks(backwardStackSaveSlots, backwardUpdates, backwardStackMem, false);
		
		// One inital run
		solve_with_current_ranks();
		
		// If we follow a rank increasing strategie, increase the ransk until we reach the targetResidual, the maxRanks or the maxIterations.
		while(residualNorm > targetResidualNorm && x.ranks() != maxRanks && (maxIterations == 0 || iteration < maxIterations)) {
			
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
		return residualNorm;
	}
	

	double ADFVariant::operator()(TTTensor& _x, const SinglePointMeasurmentSet& _measurments, PerformanceData& _perfData) const {
	   InternalSolver<SinglePointMeasurmentSet> solver(_x, _x.ranks(), _measurments, maxIterations, targetResidualNorm, minimalResidualNormDecrese, _perfData);
	   return solver.solve();
	}
	
	double ADFVariant::operator()(TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData) const {
	   InternalSolver<SinglePointMeasurmentSet> solver(_x, _maxRanks, _measurments, maxIterations, targetResidualNorm, minimalResidualNormDecrese, _perfData);
	   return solver.solve();
	}
	
	double ADFVariant::operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, PerformanceData& _perfData) const {
		InternalSolver<RankOneMeasurmentSet> solver(_x, _x.ranks(), _measurments, maxIterations, targetResidualNorm, minimalResidualNormDecrese, _perfData);
		return solver.solve();
	}
	
	double ADFVariant::operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData) const {
		InternalSolver<RankOneMeasurmentSet> solver(_x, _maxRanks, _measurments, maxIterations, targetResidualNorm, minimalResidualNormDecrese, _perfData);
		return solver.solve();
	}
}
