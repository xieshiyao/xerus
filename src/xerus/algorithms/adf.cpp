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


// workaround for broken Lapack
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
extern "C"
{
    #include <cblas.h> 
}

#ifdef __has_include
    #if __has_include(<lapacke.h>)
        #include <lapacke.h>
    #elif __has_include(<lapacke/lapacke.h>)
        #include <lapacke/lapacke.h>
	#else
		#pragma error no lapacke found
    #endif
#else
    #include <lapacke.h>
#endif


namespace xerus {
	
	double calculate_norm_of_measured_values(const std::vector<value_t>& _measurments) {
		value_t normMeasuredValues = 0;
		for(const value_t measurement : _measurments) {
			normMeasuredValues += misc::sqr(measurement);
		}
		return std::sqrt(normMeasuredValues);
	}
	
	void construct_stacks(std::unique_ptr<FullTensor[]>& stackSaveSlots, std::vector<bool>& _updates, std::unique_ptr<FullTensor*[]>& _stackMem, TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const bool _forward) {
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
			_updates[realId + corePosition*numMeasurments] = true;
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
						_updates[realId + corePosition*numMeasurments] = false;
						break;
					} else if(otherId == calculationMap[otherId + corePosition*numMeasurments]) {
						REQUIRE(_updates[otherId + corePosition*numMeasurments] == true, "I.E.");
						              calculationMap[otherId + corePosition*numMeasurments] = realId;
						_updates[otherId + corePosition*numMeasurments] = false;
						              calculationMap[realId + corePosition*numMeasurments] = realId;
						_updates[realId + corePosition*numMeasurments] = true;
						break;
					} else if( realId < calculationMap[otherId + corePosition*numMeasurments]) {
						REQUIRE(_updates[otherId + corePosition*numMeasurments] == false, "I.E.");
						size_t nextOther = calculationMap[otherId + corePosition*numMeasurments];
						              calculationMap[otherId + corePosition*numMeasurments] = realId;
						otherId = nextOther;
					} else {
						otherId = calculationMap[otherId + corePosition*numMeasurments];
					}
				}
			}
			
			for( ; position < degree; ++position, corePosition = _forward ? position : degree-1-position) {
				        calculationMap[realId + corePosition*numMeasurments] = realId;
				_updates[realId + corePosition*numMeasurments] = true;
				++numUniqueStackEntries;
			}
		}
		
		// Create the stack
		numUniqueStackEntries++; // +1 for the special positions -1 and degree.
		stackSaveSlots.reset(new FullTensor[numUniqueStackEntries]); 
		size_t usedSlots = 0; 
		stackSaveSlots[usedSlots++] = Tensor::ones({1}); // Special slot reserved for the the position -1 and degree stacks
		
		// NOTE that _stackMem contains (degree+2)*numMeasurments entries and has an offset of numMeasurments (to have space for corePosition -1).
		
		// Set links for the special entries -1 and degree
		for(size_t i = 0; i < numMeasurments; ++i) {
			_stackMem[i] = &stackSaveSlots[0];
			_stackMem[i + (degree+1)*numMeasurments] = &stackSaveSlots[0];
		}
		
		for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(calculationMap[i + corePosition*numMeasurments] == i) {
					REQUIRE(_updates[i + corePosition*numMeasurments], "I.E.");
					stackSaveSlots[usedSlots].reset({_x.rank(corePosition - (_forward ? 0 : 1))}, DONT_SET_ZERO());
					_stackMem[i + (corePosition+1)*numMeasurments] = &stackSaveSlots[usedSlots++];
				} else {
					_stackMem[i + (corePosition+1)*numMeasurments] = _stackMem[calculationMap[i + corePosition*numMeasurments] + (corePosition+1)*numMeasurments];
				}
			}
		}
		
		REQUIRE(usedSlots == numUniqueStackEntries, "Internal Error.");
		
		LOG(ADF, "We have " << numUniqueStackEntries << " unique backward stack entries. There are " << numMeasurments*degree+1 << " virtual backward stack entries.");
	}
	
	
	_inline_ void rebuild_backward_stack(const TTTensor& _x, std::vector<FullTensor>& fixedComponents, const SinglePointMeasurmentSet& _measurments, const std::vector<bool>& backwardUpdates, FullTensor* const * const backwardStack) {
		const size_t numMeasurments = _measurments.size();
		
		Index r1, r2;
		
		for(size_t corePosition = _x.degree()-1; corePosition > 0; --corePosition) {
			
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
	}
	
	/**
	 * @brief: Calculates the difference between the current and measured values at the measured positions.
	 * @note: Abuses the stack infrastructe for its caluclation. Changes only the stacks at corePosition.
	 */
	_inline_ void calculate_current_differences( std::vector<value_t>& currentDifferences,
									const TTTensor& _x,
									std::vector<FullTensor>& fixedComponents,
									const SinglePointMeasurmentSet& _measurments, 
									const std::vector<bool>& forwardUpdates, 
									const std::vector<bool>& backwardUpdates, 
									FullTensor* const * const forwardStack, 
									FullTensor* const * const backwardStack, 
									const size_t corePosition) {
		const size_t numMeasurments = _measurments.size();
		
		FullTensor currentValue({});
		
		Index r1, r2;
		
		// Prepare the single slates of the current component 
		for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
			fixedComponents[i](r1, r2) = _x.get_component(corePosition)(r1, i, r2);
		}
		
		if(corePosition <= _x.degree()/2) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(forwardUpdates[i + corePosition*numMeasurments]) {
					contract(*forwardStack[i + corePosition*numMeasurments] , *forwardStack[i + (corePosition-1)*numMeasurments], false, fixedComponents[_measurments.positions[i][corePosition]], false, 1);
				}
				contract(currentValue, *forwardStack[i + corePosition*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
				        currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
			}
		} else {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(backwardUpdates[i + corePosition*numMeasurments]) {
					contract(*backwardStack[i + corePosition*numMeasurments], fixedComponents[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
				}
				contract(currentValue, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + corePosition*numMeasurments], false, 1);
				        currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
			}
		}
	}
	
	_inline_ void calculate_deltas( std::vector<FullTensor>& deltas,
											const TTTensor& _x, 
											const SinglePointMeasurmentSet& _measurments,
											const FullTensor* const * const forwardStack, 
											const FullTensor* const * const backwardStack,
										    const std::vector<value_t>& currentDifferences,
											const size_t corePosition) {
		const size_t localLeftRank = _x.get_component(corePosition).dimensions[0];
		const size_t localRightRank = _x.get_component(corePosition).dimensions[2];
		const size_t numMeasurments = _measurments.size();

		for(size_t i = 0; i < _x.dimensions[corePosition]; ++i) {
			deltas[i].ensure_own_data();
		}
		
// 		FullTensor entryAddition({localLeftRank, localRightRank}); // Non low-level
		
		for(size_t i = 0; i < numMeasurments; ++i) {
			// Non low level variant:
// 			contract(entryAddition, *_forwardStack[i + (_corePosition-1)*numMeasurments], false, *_backwardStack[i + (_corePosition+1)*numMeasurments], false, 0);
// 			deltas[_measurments.positions[i][_corePosition]] += _currentDifferences[i]*entryAddition;
			
			// Very low level calculation of dyadic product + entriewise addition (one blas call).
			cblas_dger(CblasRowMajor, (int)localLeftRank, (int)localRightRank, 
						             currentDifferences[i]*forwardStack[i + (corePosition-1)*numMeasurments]->factor*backwardStack[i + (corePosition+1)*numMeasurments]->factor,
						             forwardStack[i + (corePosition-1)*numMeasurments]->data.get(), 1, 
						             backwardStack[i + (corePosition+1)*numMeasurments]->data.get(), 1,
						deltas[_measurments.positions[i][corePosition]].data.get(), (int)localRightRank);
		}
	}
	
	_inline_ void calculate_slicewise_norm_Py(std::vector<value_t>& PyPys,
												const TTTensor& _x, 
												const SinglePointMeasurmentSet& _measurments,
												const std::vector<bool>& forwardUpdates, 
												const std::vector<bool>& backwardUpdates, 
												FullTensor* const * const forwardStack, 
												FullTensor* const * const backwardStack,
												const std::vector<FullTensor>& deltas,
												const size_t corePosition) {
		const size_t numMeasurments = _measurments.size();
		
		FullTensor currentValue({});
		
		if(corePosition <= _x.degree()/2) {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(forwardUpdates[i + corePosition*numMeasurments]) {
					contract(*forwardStack[i + corePosition*numMeasurments] , *forwardStack[i + (corePosition-1)*numMeasurments], false, deltas[_measurments.positions[i][corePosition]], false, 1);
				}
				contract(currentValue, *forwardStack[i + corePosition*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
				PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
			}
		} else {
			for(size_t i = 0; i < numMeasurments; ++i) {
				if(backwardUpdates[i + corePosition*numMeasurments]) {
					contract(*backwardStack[i + corePosition*numMeasurments], deltas[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
				}
				contract(currentValue, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + corePosition*numMeasurments], false, 1);
				PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
			}
		}
	}
	
	_inline_ void move_forward(TTTensor& _x, std::vector<FullTensor>& fixedComponents, const SinglePointMeasurmentSet& _measurments, const std::vector<bool>& forwardUpdates, FullTensor* const * const forwardStack, const size_t corePosition) {
		const size_t numMeasurments = _measurments.size();
		
		Index r1, r2;
		
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
	
	double ADFVariant::solve(TTTensor& _x, const SinglePointMeasurmentSet& _measurments) const {
		REQUIRE(_x.is_valid_tt(), "_x must be a valid TT-Tensor.");
		
		const size_t degree = _x.degree();
		const size_t numMeasurments = _measurments.size();
		
		REQUIRE(numMeasurments > 0, "Need at very least one measurment.");
		REQUIRE(_measurments.degree() == degree, "Measurment order must coincide with _x order.");
		
		const value_t normMeasuredValues = calculate_norm_of_measured_values(_measurments.measuredValues);
		
		// We want to construct 2*numMeasurments stacks of size degree+2 containing the 
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
		
		construct_stacks(forwardStackSaveSlots, forwardUpdates, forwardStackMem, _x, _measurments, true);
		construct_stacks(backwardStackSaveSlots, backwardUpdates, backwardStackMem, _x, _measurments, false);
		
		misc::TimeMeasure timer;
		
		size_t timeA=0,timeB=0,timeC=0,timeD=0, timeE=0, timeX2=0;
		
		value_t residual = 1.0, lastResidual = 1.0;
		size_t smallResidualCount = 0;
		
		std::vector<FullTensor> fixedComponents(misc::max(_x.dimensions));
		std::vector<value_t> currentDifferences(numMeasurments);
		
		for(size_t iteration = 0; iteration < maxInterations; ++iteration) {
			// Move core back to position zero
			_x.move_core(0, true);
			
			timer.step();
			rebuild_backward_stack( _x, fixedComponents,_measurments, backwardUpdates, backwardStack);
			timeA += timer.get();
			
			// Sweep from the first to the last component
			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
				const size_t localN = _x.dimensions[corePosition];
				const size_t localLeftRank = _x.get_component(corePosition).dimensions[0];
				const size_t localRightRank = _x.get_component(corePosition).dimensions[2];
				
				// Prepare the single slates of the current component 
				for(size_t i = 0; i < localN; ++i) {
					fixedComponents[i](r1, r2) = _x.get_component(corePosition)(r1, i, r2);
				}
				
				
// 				FullTensor entryAddition({localLeftRank, localRightRank});
				FullTensor currentValue({});
				
				// NOTE we abuse the stack infrastructe for local calculations, which is ok because the stacks at corePosition are set at the end of each iteration.
				
				timer.step();
				if(corePosition <= degree/2) {
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(forwardUpdates[i + corePosition*numMeasurments]) {
							contract(*forwardStack[i + corePosition*numMeasurments] , *forwardStack[i + (corePosition-1)*numMeasurments], false, fixedComponents[_measurments.positions[i][corePosition]], false, 1);
						}
						contract(currentValue, *forwardStack[i + corePosition*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
						currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
					}
				} else {
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(backwardUpdates[i + corePosition*numMeasurments]) {
							contract(*backwardStack[i + corePosition*numMeasurments], fixedComponents[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
						}
						contract(currentValue, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + corePosition*numMeasurments], false, 1);
						currentDifferences[i] = (_measurments.measuredValues[i]-currentValue[0]);
					}
				}
				timeB += timer.get();
				
				timer.step();
				std::vector<FullTensor> deltas(localN, FullTensor({localLeftRank, localRightRank}));
				for(size_t i = 0; i < localN; ++i) {
					deltas[i].ensure_own_data();
				}
				
				for(size_t i = 0; i < numMeasurments; ++i) {
					// Non low level variant:
// 					contract(entryAddition, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 0);
// 					deltas[_measurments.positions[i][corePosition]] += currentDifferences[i]*entryAddition;
					
					// Very low level calculation of dyadic product + entriewise addition (one blas call).
					cblas_dger(CblasRowMajor, (int)localLeftRank, (int)localRightRank, 
							    currentDifferences[i]*forwardStack[i + (corePosition-1)*numMeasurments]->factor*backwardStack[i + (corePosition+1)*numMeasurments]->factor,
								forwardStack[i + (corePosition-1)*numMeasurments]->data.get(), 1, 
								backwardStack[i + (corePosition+1)*numMeasurments]->data.get(), 1,
							    deltas[_measurments.positions[i][corePosition]].data.get(), (int)localRightRank);
				}
				timeC += timer.get();
				
				timer.step();
				// Calculate ||P(y)||^2 for each slice seperately, where y is the update direction.
				std::vector<value_t> PyPys(localN, 0.0);
				calculate_slicewise_norm_Py(PyPys, _x, _measurments, forwardUpdates, backwardUpdates, forwardStack, backwardStack, deltas, corePosition);
				
				
					/*
				if(corePosition <= _x.degree()/2) {
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(forwardUpdates[i + corePosition*numMeasurments]) {
							contract(*forwardStack[i + corePosition*numMeasurments] , *forwardStack[i + (corePosition-1)*numMeasurments], false, deltas[_measurments.positions[i][corePosition]], false, 1);
						}
						contract(currentValue, *forwardStack[i + corePosition*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
						PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
					}
				} else {
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(backwardUpdates[i + corePosition*numMeasurments]) {
							contract(*backwardStack[i + corePosition*numMeasurments], deltas[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
						}
						contract(currentValue, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + corePosition*numMeasurments], false, 1);
						PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
					}
				}*/
				
				/*
				if(corePosition <= degree/2) {
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(forwardUpdates[i + corePosition*numMeasurments]) {
							contract(*forwardStack[i + corePosition*numMeasurments] , *forwardStack[i + (corePosition-1)*numMeasurments], false, deltas[_measurments.positions[i][corePosition]], false, 1);
						}
						contract(currentValue, *forwardStack[i + corePosition*numMeasurments], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
						PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
					}
				} else {
					for(size_t i = 0; i < numMeasurments; ++i) {
						if(backwardUpdates[i + corePosition*numMeasurments]) {
							contract(*backwardStack[i + corePosition*numMeasurments], deltas[_measurments.positions[i][corePosition]], false, *backwardStack[i + (corePosition+1)*numMeasurments], false, 1);
						}
						contract(currentValue, *forwardStack[i + (corePosition-1)*numMeasurments], false, *backwardStack[i + corePosition*numMeasurments], false, 1);
						PyPys[_measurments.positions[i][corePosition]] += misc::sqr(currentValue[0]);
					}
				}*/
				timeD += timer.get();
				
				// Update each slice seperately
				for(size_t j = 0; j < localN; ++j) {
					// Calculate <P(y), P(X-B)> = ||deltaPlus||^2.
					const value_t PyR = misc::sqr(frob_norm(deltas[j]));
					
					// Update
					_x.component(corePosition)(r1, i1, r2) = _x.component(corePosition)(r1, i1, r2) + (PyR/PyPys[j])*Tensor::dirac({localN}, {j})(i1)*deltas[j](r1, r2);
				}
				
				timer.step();
				// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
				if(corePosition+1 < degree) {
					move_forward( _x, fixedComponents, _measurments, forwardUpdates, forwardStack, corePosition);
				}
				timeE += timer.get();
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
			
			LOG(ADF, "Rank " << _x.ranks() << " Itr: " << iteration << " Residual: " << std::scientific << residual << " Rel. Residual change: " << residual/lastResidual);
			
			if(residual <= convergenceEpsilon || smallResidualCount > 3) {
				LOG(time, " A: " << std::round(100.0*double(timeA)/double(timer.getTotal())) << " B: " << std::round(100.0*double(timeB)/double(timer.getTotal())) << " C: " << std::round(100.0*double(timeC)/double(timer.getTotal())) << " D: " << std::round(100.0*double(timeD)/double(timer.getTotal()))<< " E: " << std::round(100.0*double(timeE)/double(timer.getTotal()))<< " X2: " << std::round(100.0*double(timeX2)/double(timer.getTotal())));
				return residual;
			}
		}
				LOG(time, " A: " << std::round(100.0*double(timeA)/double(timer.getTotal())) << " B: " << std::round(100.0*double(timeB)/double(timer.getTotal())) << " C: " << std::round(100.0*double(timeC)/double(timer.getTotal())) << " D: " << std::round(100.0*double(timeD)/double(timer.getTotal()))<< " E: " << std::round(100.0*double(timeE)/double(timer.getTotal()))<< " X2: " << std::round(100.0*double(timeX2)/double(timer.getTotal())));
		return residual;
	}
}
