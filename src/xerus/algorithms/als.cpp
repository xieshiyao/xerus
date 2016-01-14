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
 * @brief Implementation of the ALS variants.
 */

#include <xerus/algorithms/als.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/indexedTensorList.h>
#include <xerus/tensorNetwork.h>
 
#include <xerus/indexedTensorMoveable.h>
#include <xerus/indexedTensor_tensor_factorisations.h>

namespace xerus {

    void ALSVariant::lapack_solver(const TensorNetwork &_A, Tensor &_x, const Tensor &_b) {
        Tensor A(_A);
        Index i,j;
        _x(i&0) = _b(j&0) / A(j/2, i/2);
//         REQUIRE(_x.degree() <= 3, "dmrg not yet implemented in lapack_solver");// TODO split result into d-2 tensors -> choose correct tensor as core!
    }
    
    /// @brief Finds the range of notes that need to be optimized and orthogonalizes @a _x properly
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
			
			Tensor& curComponent = _x.component(firstOptimizedIndex);
			curComponent.reinterpret_dimensions({curComponent.dimensions[0]*curComponent.dimensions[1], curComponent.dimensions[2]});
			curComponent(r1,n1,r2) = curComponent(r1,cr1) * _x.get_component(firstOptimizedIndex+1)(cr1,n1,r2);
			_x.set_component(firstOptimizedIndex+1, std::move(curComponent));
			
			//TODO sparse
			_x.set_component(firstOptimizedIndex, Tensor(
				{dimensionProd, localDim, newDimensionProd},
				[&](const std::vector<size_t> &_idx){
					if (_idx[0]*localDim + _idx[1] == _idx[2]) {
						return 1.0;
					} else {
						return 0.0;
					}
				})
			);
			
			_x.require_correct_format();
			
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
			
			Tensor& curComponent = _x.component(firstNotOptimizedIndex-1);
			curComponent.reinterpret_dimensions({curComponent.dimensions[0], curComponent.dimensions[1] * curComponent.dimensions[2]});
			curComponent(r1,n1,r2) = _x.get_component(firstNotOptimizedIndex-2)(r1,n1,cr1) * curComponent(cr1,r2);
			_x.set_component(firstNotOptimizedIndex-2, std::move(curComponent));
			
			//TODO sparse
			_x.set_component(firstNotOptimizedIndex-1, Tensor(
				{newDimensionProd, localDim, dimensionProd},
				[&](const std::vector<size_t> &_idx){
					if (_idx[0] == _idx[1]*dimensionProd + _idx[2]) {
						return 1.0;
					} else {
						return 0.0;
					}
				})
			);

			_x.require_correct_format();
			
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
		
		return std::pair<size_t, size_t>(firstOptimizedIndex, firstNotOptimizedIndex);
	}

    double ALSVariant::solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, value_t _convergenceEpsilon, PerformanceData &_perfData) const {
		const TTOperator &_A = *_Ap;
        LOG(ALS, "ALS("<< sites << ", " << minimumLocalResidual <<") called");
		#ifndef DISABLE_RUNTIME_CHECKS_
			REQUIRE(_x.is_valid_tt(), "");
			REQUIRE(_b.is_valid_tt(), "");
			REQUIRE(_x.degree() > 0, "");
			REQUIRE(_x.dimensions == _b.dimensions, "");
			REQUIRE(sites == 1, "DMRG and n-site-dmrg not yet implemented!"); // TODO
			
			if (_Ap) {
				REQUIRE(_A.is_valid_tt(), "");
				REQUIRE(_A.dimensions.size() == _b.dimensions.size()*2, "");
				for (size_t i=0; i<_x.dimensions.size(); ++i) {
					REQUIRE(_A.dimensions[i] == _x.dimensions[i], "");
					REQUIRE(_A.dimensions[i+_A.degree()/2] == _x.dimensions[i], "");
				}
			}
        #endif
        
        if (_Ap) {
			_perfData << "ALS for ||A*x - b||^2, x.dimensions: " << _x.dimensions << '\n'
					<< "A.ranks: " << _A.ranks() << '\n'
					<< "x.ranks: " << _x.ranks() << '\n'
					<< "b.ranks: " << _b.ranks() << '\n'
					<< "maximum number of half sweeps: " << _numHalfSweeps << '\n'
					<< "convergence epsilon: " << _convergenceEpsilon << '\n';
		} else {
			_perfData << "ALS for ||x - b||^2, x.dimensions: " << _x.dimensions << '\n'
					<< "x.ranks: " << _x.ranks() << '\n'
					<< "b.ranks: " << _b.ranks() << '\n'
					<< "maximum number of half sweeps: " << _numHalfSweeps << '\n'
					<< "convergence epsilon: " << _convergenceEpsilon << '\n';
		}
		_perfData.start();
        
        const size_t d = _x.degree();
        Index cr1, cr2, cr3, r1, r2, r3, n1, n2;
        
        std::vector<Tensor> xAxL, xAxR;
        std::vector<Tensor> bxL, bxR;
        
		bool cannoAtTheEnd = _x.cannonicalized;
		size_t corePosAtTheEnd = _x.corePosition;
		
		std::pair<size_t, size_t> optimizedRange = prepare_x_for_als(_x);

        // Create stacks of contracted network-parts
		Tensor tmpA({1,1,1}, [](){return 1.0;});
		Tensor tmpB({1,1}, [](){return 1.0;});
		
		xAxL.emplace_back(tmpA);
		bxL.emplace_back(tmpB);
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
        Tensor BTilde;
        value_t lastEnergy2 = 1e102;
        value_t lastEnergy = 1e101;
        value_t energy = 1e100;
        bool walkingRight = true;
        bool changedSmth = false;
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
					Tensor res;
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
					Tensor res;
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
        if (_perfData) {
			_perfData.stop_timer();
			energy = energy_f();
			_perfData.continue_timer();
			_perfData.add(energy);
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
					Tensor tmpX;
					localSolver(ATilde, tmpX, BTilde);
					_x.set_component(currIndex, tmpX);
					changedSmth = true;
				}
			} else {
				_x.set_component(currIndex, std::move(BTilde));
				changedSmth = true;
			}
			
            if (_perfData) {
                _perfData.stop_timer();
				energy = energy_f();
				_perfData.continue_timer();
				_perfData.add(energy);
                LOG(ALS, "Calculated residual for perfData after tensor "<< currIndex <<": " << energy << " ( " << std::abs(energy - lastEnergy) << " vs " << _convergenceEpsilon << " ) ");
                if (printProgress) {
                    std::cout << "optimized tensor "<< currIndex << ": " << std::scientific << energy << " ( \t" << std::abs(energy - lastEnergy) << " vs \t" << _convergenceEpsilon << " ) \r" << std::flush;
					std::cout << "                                                                               \r"; // note: not flushed so it will only erase content on next output
                }
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
					if (cannoAtTheEnd && preserveCorePosition) {
						_x.move_core(corePosAtTheEnd);
					}
					return energy;
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
	
	
    const ALSVariant ALS(1, 0, EPSILON, ALSVariant::lapack_solver);
	
    const ALSVariant DMRG(2, 0, EPSILON, ALSVariant::lapack_solver);
}
