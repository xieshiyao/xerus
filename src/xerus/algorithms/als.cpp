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

#include <xerus/algorithms/als.h>
#include <xerus/basic.h>
#include <xerus/indexedTensorList.h>
#include <xerus/tensorNetwork.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_tensor_factorisations.h>

namespace xerus {

    void ALSVariant::lapack_solver(const TensorNetwork &_A, Tensor &_x, const Tensor &_b) {
        FullTensor A(_A);
        Index i,j;
        _x(i&0) = _b(j&0) / A(j/2, i/2);
        REQUIRE(_x.degree() <= 3, "dmrg not yet implemented in lapack_solver");// TODO split result into d-2 tensors -> choose correct tensor as core!
    }
    
    /// @brief finds the range of notes that need to be optimized and orthogonalizes @a _x properly
    /// @details finds full-rank nodes (these can wlog be set to identity and need not be optimized)
    std::pair<size_t, size_t> prepare_x_for_als(TTTensor &_x) {
		bool cannoAtTheEnd = _x.cannonicalized;
		size_t corePosAtTheEnd = _x.corePosition;
		const size_t d = _x.degree();
		Index r1,r2,n1,cr1;
		
		size_t firstOptimizedIndex = 0;
		size_t dimensionProd = 1;
		while (firstOptimizedIndex + 1 < d) {
			const size_t localDim = _x.dimensions[firstOptimizedIndex];
			size_t newDimensionProd = dimensionProd * localDim;
			if (_x.rank(firstOptimizedIndex) < newDimensionProd) {
				break;
			}
			//TODO sparse
			std::unique_ptr<Tensor> tmp(new FullTensor(static_cast<const FullTensor &>(_x.get_component(firstOptimizedIndex))));
			
			_x.set_component(firstOptimizedIndex, FullTensor(
				{dimensionProd, localDim, newDimensionProd},
				[&](const std::vector<size_t> &_idx){
					if (_idx[0]*localDim + _idx[1] == _idx[2]) {
						return 1.0;
					} else {
						return 0.0;
					}
				})
			);
			tmp->reinterpret_dimensions({tmp->dimensions[0]*tmp->dimensions[1], tmp->dimensions[2]});
			(*tmp)(r1,n1,r2) = (*tmp)(r1,cr1) * _x.get_component(firstOptimizedIndex+1)(cr1,n1,r2);
			_x.set_component(firstOptimizedIndex+1, std::move(tmp));
			REQUIRE(_x.is_in_expected_format(), "ie");
			
			firstOptimizedIndex += 1;
			dimensionProd = newDimensionProd;
		}
		size_t firstNotOptimizedIndex = d;
		dimensionProd = 1;
		while (firstNotOptimizedIndex > firstOptimizedIndex+1) {
			const size_t localDim = _x.dimensions[firstNotOptimizedIndex-1];
			size_t newDimensionProd = dimensionProd * localDim;
			if (_x.rank(firstNotOptimizedIndex-2) < newDimensionProd) {
				break;
			}
			//TODO sparse
			std::unique_ptr<Tensor> tmp(new FullTensor(static_cast<const FullTensor &>(_x.get_component(firstNotOptimizedIndex-1))));
			
			_x.set_component(firstNotOptimizedIndex-1, FullTensor(
				{newDimensionProd, localDim, dimensionProd},
				[&](const std::vector<size_t> &_idx){
					if (_idx[0] == _idx[1] + _idx[2]*localDim) {
						return 1.0;
					} else {
						return 0.0;
					}
				})
			);
			tmp->reinterpret_dimensions({tmp->dimensions[0], tmp->dimensions[1] * tmp->dimensions[2]});
			(*tmp)(r1,n1,r2) = _x.get_component(firstNotOptimizedIndex-2)(r1,n1,cr1) * (*tmp)(cr1,r2);
			_x.set_component(firstNotOptimizedIndex-2, std::move(tmp));
			REQUIRE(_x.is_in_expected_format(), "ie");
			
			firstNotOptimizedIndex -= 1;
			dimensionProd = newDimensionProd;
		}
		
		if (cannoAtTheEnd && corePosAtTheEnd < firstOptimizedIndex) {
			_x.assume_core_position(firstOptimizedIndex);
		} else {
			if (cannoAtTheEnd && corePosAtTheEnd >= firstNotOptimizedIndex) {
				_x.assume_core_position(firstNotOptimizedIndex-1);
			}
			
			_x.move_core(firstOptimizedIndex, true);
		}
		// as we will access the nodes individually in the algorithm, make sure there is no global factor
		(*_x.nodes[firstOptimizedIndex].tensorObject) *= _x.factor;
		_x.factor = 1.0;
		
		return std::pair<size_t, size_t>(firstOptimizedIndex, firstNotOptimizedIndex);
	}

    double ALSVariant::solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, value_t _convergenceEpsilon,  std::vector<value_t> *_perfData) const {
		const TTOperator &_A = *_Ap;
        LOG(ALS, "ALS("<< sites << ", " << minimumLocalResidual <<") called");
		#ifndef DISABLE_RUNTIME_CHECKS_
			REQUIRE(_x.degree() > 0, "");
			REQUIRE(_x.dimensions == _b.dimensions, "");
			REQUIRE(sites == 1, "DMRG and n-site-dmrg not yet implemented!"); // TODO
			
			if (_Ap) {
				REQUIRE(_A.dimensions.size() == _b.dimensions.size()*2, "");
				for (size_t i=0; i<_x.dimensions.size(); ++i) {
					REQUIRE(_A.dimensions[i] == _x.dimensions[i], "");
					REQUIRE(_A.dimensions[i+_A.degree()/2] == _x.dimensions[i], "");
				}
			}
        #endif
        
        const size_t d = _x.degree();
        Index cr1, cr2, cr3, r1, r2, r3, n1, n2;
        
        std::vector<FullTensor> xAxL, xAxR;
        std::vector<FullTensor> bxL, bxR;
        
		bool cannoAtTheEnd = _x.cannonicalized;
		size_t corePosAtTheEnd = _x.corePosition;
		
		std::pair<size_t, size_t> optimizedRange = prepare_x_for_als(_x);

        // Create stacks of contracted network-parts
		FullTensor tmpA({1,1,1}, [](){return 1.0;});
		FullTensor tmpB({1,1}, [](){return 1.0;});
		
		xAxL.emplace_back(tmpA);
		bxL.emplace_back(tmpB);
		// as we access all the nodes manualy we have to take care of the global factors in _A and _b outselves
		if (_Ap) {
			tmpA *= _A.factor;
		}
		tmpB *= _b.factor;
		xAxR.emplace_back(tmpA);
		bxR.emplace_back(tmpB);
		
        for (size_t i = d-1; i > optimizedRange.first + sites - 1; --i) {
			if (_Ap) {
				tmpA(r1, r2, r3) = xAxR.back()(cr1, cr2, cr3) * _x.get_component(i)(r1, n1, cr1) * _A.get_component(i)(r2, n1, n2, cr2) * _x.get_component(i)(r3, n2, cr3);
				xAxR.emplace_back(std::move(tmpA));
			}
			tmpB(r1, r2) = bxR.back()(cr1, cr2) * _b.get_component(i)(r1, n1, cr1) * _x.get_component(i)(r2, n1, cr2);
            bxR.emplace_back(std::move(tmpB));
        }
        for (size_t i = 0; i < optimizedRange.first; ++i) {
			if (_Ap) {
				tmpA(r1,r2,r3) = xAxL.back()(cr1,cr2,cr3) * _x.get_component(i)(cr1,n1,r1) * _A.get_component(i)(cr2, n1, n2, r2) * _x.get_component(i)(cr3,n2,r3);
				xAxL.emplace_back(std::move(tmpA));
			}
			tmpB(r1,r2) = bxL.back()(cr1,cr2) * _b.get_component(i)(cr1, n1, r1) * _x.get_component(i)(cr2, n1, r2);
			bxL.emplace_back(std::move(tmpB));
        }
        
        TensorNetwork ATilde;
        FullTensor BTilde;
        value_t lastEnergy2 = 1e102;
        value_t lastEnergy = 1e101;
        value_t energy = 1e100;
        bool walkingRight = true;
        bool changedSmth = false;
		bool done = false;
        size_t currIndex = optimizedRange.first;
		size_t halfSweepCount = 0;
		
		std::function<value_t()> energy_f;
		if (_Ap) {
			if (useResidualForEndCriterion) {
				energy_f = [&](){
					return frob_norm(_A(n1/2,n2/2)*_x(n2&0) - _b(n1&0));
				};
			} else {
				energy_f = [&](){
					FullTensor res;
					// 0.5*<x,Ax> - <x,b>
					res() = 0.5*xAxR.back()(cr1, cr2, cr3) 
							* _x.get_component(currIndex)(r1, n1, cr1) * _A.get_component(currIndex)(r2, n1, n2, cr2) * _x.get_component(currIndex)(r3, n2, cr3) 
							* xAxL.back()(r1, r2, r3)
							- bxR.back()(cr1, cr2) 
							* _b.get_component(currIndex)(r1, n1, cr1) * _x.get_component(currIndex)(r2, n1, cr2)
							* bxL.back()(r1, r2);
					return res[0];
				};
			}
		} else {
			if (useResidualForEndCriterion) {
				energy_f = [&](){
					return frob_norm(_x(n1&0) - _b(n1&0));
				};
			} else {
				energy_f = [&](){
					FullTensor res;
					// 0.5*<x,Ax> - <x,b> = 0.5*|x_i|^2 - <x,b>
					res() = 0.5*_x.get_component(currIndex)(r1, n1, cr1) * _x.get_component(currIndex)(r1, n1, cr1) 
							- bxR.back()(cr1, cr2) 
							* _b.get_component(currIndex)(r1, n1, cr1) * _x.get_component(currIndex)(r2, n1, cr2)
							* bxL.back()(r1, r2);
					return res[0];
				};
			}
		}
        
        // Calculate initial residual
        if (_perfData != nullptr) {
			energy = energy_f();
			_perfData->push_back(energy);
			LOG(ALS, "calculated residual for perfData: " << energy);
        }
		
        while (true) {
			LOG(ALS, "Starting to optimize index " << currIndex);
			
			// Calculate Atilde and Btilde
			if (_Ap) {
				ATilde(r1,n1,cr1, r3,n2,cr3) = xAxL.back()(r1,r2,r3) * _A.get_component(currIndex)(r2, n1, n2, cr2) * xAxR.back()(cr1, cr2, cr3);
			}
			BTilde(r2,n1,cr2) = bxL.back()(r1,r2) * _b.get_component(currIndex)(r1, n1, cr1) * bxR.back()(cr1,cr2);
			
			// Change component tensor if the local residual is large enough
			if (_Ap) {
				if (minimumLocalResidual <= 0 || frob_norm(IndexedTensorMoveable<Tensor>(ATilde(r1^3, r2^3)*_x.get_component(currIndex)(r2^3)) - BTilde(r1^3)) > minimumLocalResidual) {
					FullTensor tmpX;
					localSolver(ATilde, tmpX, BTilde);
					_x.set_component(currIndex, tmpX);
					changedSmth = true;
				}
			} else {
				_x.set_component(currIndex, std::move(BTilde));
			}
			
            if (_perfData) {
                energy = energy_f();
				_perfData->push_back(energy);
                LOG(ALS, "Calculated residual for perfData after tensor "<< currIndex <<": " << energy << " ( " << std::abs(energy - lastEnergy) << " vs " << _convergenceEpsilon << " ) ");
                if (printProgress) {
					std::cout << "                                                                      \r";
                    std::cout << "optimized tensor "<< currIndex << ": " << std::scientific << energy << " ( \t" << std::abs(energy - lastEnergy) << " vs \t" << _convergenceEpsilon << " ) \r" << std::flush;
                }
            }
			if (done && (corePosAtTheEnd == currIndex 
						|| (corePosAtTheEnd < optimizedRange.first && currIndex == optimizedRange.first)
						|| (corePosAtTheEnd >= optimizedRange.second && currIndex == optimizedRange.second-1))) {
				_x.move_core(corePosAtTheEnd);
				if (!_perfData) {
					energy = energy_f();
				}
				return energy;
			}
            
            // Are we done with the sweep?
            if ((!walkingRight && currIndex==optimizedRange.first) 
				|| (walkingRight && currIndex==optimizedRange.second-sites)) 
			{
                LOG(ALS, "Sweep Done");
				halfSweepCount += 1;
                
                if (!_perfData) {
                    LOG(ALS, "Calculating energy for loop condition");
                    energy = energy_f();
                }
                
                LOG(ALS, "Stats: " << xAxL.size() << " " << xAxR.size() << " " << bxL.size() << " " << bxR.size() << " energy: " << energy << " deltas: " << (1-lastEnergy/energy) << " " << (1-lastEnergy2/energy));
				LOG(ALS, (lastEnergy-energy) << " vs " << _convergenceEpsilon);
				
                // Conditions for loop termination
                if (!changedSmth || halfSweepCount == _numHalfSweeps || std::abs(lastEnergy-energy) < _convergenceEpsilon || std::abs(lastEnergy2-energy) < _convergenceEpsilon || (optimizedRange.second - optimizedRange.first<=sites)) {
                    // we are done! yay
                    LOG(ALS, "ALS done, " << energy << " " << lastEnergy << " " << std::abs(lastEnergy2-energy) << " " << std::abs(lastEnergy-energy) << " < " << _convergenceEpsilon);
					if (!cannoAtTheEnd || !preserveCorePosition || currIndex == corePosAtTheEnd) return energy;
					done = true;
                }
                
                lastEnergy2 = lastEnergy;
                lastEnergy = energy;
                walkingRight = !walkingRight;
                changedSmth = false;
                LOG(ALS, "Start sweep " << (walkingRight?"right":"left"));
            }
            
            
            if (walkingRight) {
				// Move core to next position
				_x.move_core(currIndex+1, true);
				
                // Move one site to the right
				if (_Ap) {
					xAxR.pop_back();
					tmpA(r1,r2,r3) = xAxL.back()(cr1,cr2,cr3) * _x.get_component(currIndex)(cr1,n1,r1) * _A.get_component(currIndex)(cr2, n1, n2, r2) * _x.get_component(currIndex)(cr3,n2,r3);
					xAxL.emplace_back(std::move(tmpA));
				}
                
				bxR.pop_back();
                tmpB(r1,r2) = bxL.back()(cr1,cr2) * _b.get_component(currIndex)(cr1, n1, r1) * _x.get_component(currIndex)(cr2, n1, r2);
                bxL.emplace_back(std::move(tmpB));
                currIndex++;
            } else {
				// Move core to next position
				_x.move_core(currIndex-1, true);
			
                // move one site to the left
				if (_Ap) {
					xAxL.pop_back();
					tmpA(r1,r2,r3) = xAxR.back()(cr1,cr2,cr3) * _x.get_component(currIndex)(r1,n1,cr1) * _A.get_component(currIndex)(r2, n1, n2, cr2) * _x.get_component(currIndex)(r3,n2,cr3);
					xAxR.emplace_back(std::move(tmpA));
				}
                
                bxL.pop_back();
				tmpB(r1,r2) = bxR.back()(cr1,cr2) * _b.get_component(currIndex)(r1, n1, cr1) * _x.get_component(currIndex)(r2, n1, cr2);
                bxR.emplace_back(std::move(tmpB));
                currIndex--;
            }
        }
    }
	
}
