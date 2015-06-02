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

// #include <xerus.h>
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

    double ALSVariant::operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon,  std::vector<value_t> *_perfData) const {
        LOG(ALS, "ALS("<< sites << ", " << minimumLocalResidual <<") called");
		#ifndef DISABLE_RUNTIME_CHECKS_
			REQUIRE(_x.degree() > 0, "");
			REQUIRE(_x.dimensions == _b.dimensions, "");
			REQUIRE(_A.dimensions.size() == _b.dimensions.size()*2, "");
			REQUIRE(sites == 1, "DMRG and n-site-dmrg not yet implemented!"); // TODO release critical?
			
			for (size_t i=0; i<_x.dimensions.size(); ++i) {
				REQUIRE(_A.dimensions[i] == _x.dimensions[i], "");
				REQUIRE(_A.dimensions[i+_A.degree()/2] == _x.dimensions[i], "");
			}
        #endif
        
        const size_t d = _x.degree();
        Index cr1, cr2, cr3, r1, r2, r3, n1, n2;
        
        std::vector<FullTensor> xAxL, xAxR;
        std::vector<FullTensor> bxL, bxR;
        
        _x.cannonicalize_left();

        // Create stack of contracted network-parts left part
		FullTensor tmpA;
		FullTensor tmpB;
		
        tmpA(r1, r2, r3) = (*_x.nodes[0].tensorObject)(r1) * (*_A.nodes[0].tensorObject)(r2) * (*_x.nodes[0].tensorObject)(r3);
		xAxL.emplace_back(tmpA);
		
		tmpB(r1, r2) = (*_b.nodes[0].tensorObject)(r1) * (*_x.nodes[0].tensorObject)(r2);
		bxL.emplace_back(tmpB);
		
        // Create stack of contracted network-parts right part
        tmpA(r1, r2, r3) = (*_x.nodes[d+1].tensorObject)(r1) * (*_A.nodes[d+1].tensorObject)(r2) * (*_x.nodes[d+1].tensorObject)(r3);
		xAxR.emplace_back(tmpA);
		
		tmpB(r1, r2) = (*_b.nodes[d+1].tensorObject)(r1) * (*_x.nodes[d+1].tensorObject)(r2);
		bxR.emplace_back(tmpB);
		
        for (size_t i = _x.degree()-1; i > sites-1; --i) {
			tmpA(r1, r2, r3) = xAxR.back()(cr1, cr2, cr3) * _x.get_component(i)(r1, n1, cr1) * _A.get_component(i)(r2, n1, n2, cr2) * _x.get_component(i)(r3, n2, cr3);
			tmpB(r1, r2) = bxR.back()(cr1, cr2) * _b.get_component(i)(r1, n1, cr1) * _x.get_component(i)(r2, n1, cr2);
            xAxR.emplace_back(std::move(tmpA));
            bxR.emplace_back(std::move(tmpB));
        }
        
        // Calculate initial residual
        if (_perfData != nullptr) {
            _perfData->push_back(frob_norm(_A(n1/2,n2/2)*_x(n2&0) - _b(n1&0)));
            LOG(ALS, "calculated residual for perfData: " << _perfData->back());
        }
        
        TensorNetwork ATilde;
        FullTensor BTilde;
        value_t lastEnergy2 = 1e102;
        value_t lastEnergy = 1e101;
        value_t energy = 1e100;
        bool walkingRight = true;
        bool changedSmth = false;
        size_t currIndex = 0;
        while (true) {
			LOG(ALS, "Starting to optimize index " << currIndex);
			
			// Calculate Atilde and Btilde
			ATilde(r1,n1,cr1, r3,n2,cr3) = xAxL.back()(r1,r2,r3) * _A.get_component(currIndex)(r2, n1, n2, cr2) * xAxR.back()(cr1, cr2, cr3);
			BTilde(r2,n1,cr2) = bxL.back()(r1,r2) * _b.get_component(currIndex)(r1, n1, cr1) * bxR.back()(cr1,cr2);
			
			// Change component tensor if the local residual is large enough
			if (minimumLocalResidual <= 0 || frob_norm(IndexedTensorMoveable<Tensor>(ATilde(r1^3, r2^3)*_x.get_component(currIndex)(r2^3)) - BTilde(r1^3)) > minimumLocalResidual) {
				FullTensor tmpX;
				localSolver(ATilde, tmpX, BTilde);
				_x.set_component(currIndex, tmpX);
				changedSmth = true;
			}
			
            if (_perfData != nullptr) {
                _perfData->push_back(frob_norm(_A(n1/2,n2/2)*_x(n2&0) - _b(n1&0)));
                LOG(ALS, "Calculated residual for perfData after tensor "<< currIndex <<": " << _perfData->back() << " ( " << std::abs(_perfData->back() - energy) << " vs " << _convergenceEpsilon << " ) ");
                if (printProgress) {
                    std::cout << "optimized tensor "<< currIndex << ": " << _perfData->back() << " ( \t" << std::abs(_perfData->back() - energy) << " vs \t" << _convergenceEpsilon << " ) \r" << std::flush;
                }
            }
            
            // Are we done with the sweep?
            if ((!walkingRight && currIndex==0) || (walkingRight && currIndex==d-sites)) {
                LOG(ALS, "Sweep Done");
                
                lastEnergy2 = lastEnergy;
                lastEnergy = energy;
                if (_perfData != nullptr) {
                    energy = _perfData->back(); // Energy already calculated
                } else {
                    LOG(ALS, "Calculating energy for loop condition");
                    energy = frob_norm(_A(n1/2,n2/2)*_x(n2&0) - _b(n1&0));
                }
                
                LOG(ALS, "Stats: " << xAxL.size() << " " << xAxR.size() << " " << bxL.size() << " " << bxR.size() << " energy: " << energy << " deltas: " << (1-lastEnergy/energy) << " " << (1-lastEnergy2/energy));
                
                // Conditions for loop termination
                if (!changedSmth || std::abs(lastEnergy-energy) < _convergenceEpsilon || std::abs(lastEnergy2-energy) < _convergenceEpsilon || d<=sites) {
                    // we are done! yay
                    LOG(ALS, "ALS done, " << energy << " " << lastEnergy << " " << std::abs(lastEnergy2-energy) << " " << std::abs(lastEnergy-energy) << " < " << _convergenceEpsilon);
					_x.cannonicalize_left();
                    return energy; 
                }
                
                walkingRight = !walkingRight;
                changedSmth = false;
                LOG(ALS, "Start sweep " << (walkingRight?"right":"left"));
            }
            
            
            if (walkingRight) {
				// Move core to next position
				_x.move_core(currIndex+1);
				
                // Move one site to the right
                xAxR.pop_back();
                bxR.pop_back();
                
                LOG(ALS, "Calc xAxL");
				tmpA(r1,r2,r3) = xAxL.back()(cr1,cr2,cr3) * _x.get_component(currIndex)(cr1,n1,r1) * _A.get_component(currIndex)(cr2, n1, n2, r2) * _x.get_component(currIndex)(cr3,n2,r3);
                xAxL.emplace_back(std::move(tmpA));
                
                LOG(ALS, "Calc BxL");
				tmpB(r1,r2) = bxL.back()(cr1,cr2) * _b.get_component(currIndex)(cr1, n1, r1) * _x.get_component(currIndex)(cr2, n1, r2);
                bxL.emplace_back(std::move(tmpB));
                currIndex++;;
            } else {
				// Move core to next position
				_x.move_core(currIndex-1);
			
                // move one site to the left
                xAxL.pop_back();
                bxL.pop_back();
                
                LOG(ALS, "Calc xAxR");
				tmpA(r1,r2,r3) = xAxR.back()(cr1,cr2,cr3) * _x.get_component(currIndex)(r1,n1,cr1) * _A.get_component(currIndex)(r2, n1, n2, cr2) * _x.get_component(currIndex)(r3,n2,cr3);
                xAxR.emplace_back(std::move(tmpA));
                
                LOG(ALS, "Calc bxR");
				tmpB(r1,r2) = bxR.back()(cr1,cr2) * _b.get_component(currIndex)(r1, n1, cr1) * _x.get_component(currIndex)(r2, n1, cr2);
                bxR.emplace_back(std::move(tmpB));
                currIndex--;
            }
        }
    }
	
	
	
	
	double ProjectionALSVariant::solve(TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, value_t _convergenceEpsilon,  std::vector<value_t> *_perfData) const {
		LOG(ALS, "pALS("<< _numHalfSweeps << ", " << _convergenceEpsilon <<") called");
		REQUIRE(_x.degree() > 0, "");
		REQUIRE(_x.dimensions == _b.dimensions, "");
		
		const size_t d = _x.degree();
		Index cr1, cr2, r1, r2, n1;
		
		std::vector<FullTensor> bxL, bxR;
		
		bool cannoAtTheEnd = _x.cannonicalized;
		size_t corePosAtTheEnd = _x.corePosition;
		
		_x.cannonicalize_left();

		// Create stack of contracted network-parts
		FullTensor tmp({1,1}, [](){return 1.0;});
		bxL.emplace_back(tmp);
		bxR.emplace_back(tmp);
		
		for (size_t i = _x.degree()-1; i > 0; --i) {
			tmp(r1, r2) = bxR.back()(cr1, cr2) * _b.get_component(i)(r1, n1, cr1) * _x.get_component(i)(r2, n1, cr2);
			bxR.emplace_back(std::move(tmp));
		}
		
		// Calculate initial residual
		if (_perfData != nullptr) {
			_perfData->push_back(frob_norm(_x(n1&0) - _b(n1&0)));
			LOG(ALS, "calculated residual for perfData: " << _perfData->back());
		}
		
		FullTensor BTilde;
		value_t lastEnergy2 = 1e102;
		value_t lastEnergy = 1e101;
		value_t energy = 1e100;
		bool walkingRight = true;
		size_t currIndex = 0;
		size_t halfSweepCount = 0;
		while (true) {
			LOG(ALS, "Starting to optimize index " << currIndex);
			
			// Calculate Atilde and Btilde
			BTilde(r2,n1,cr2) = bxL.back()(r1,r2) * _b.get_component(currIndex)(r1, n1, cr1) * bxR.back()(cr1,cr2);
			
			_x.set_component(currIndex, std::move(BTilde));
			
			if (_perfData != nullptr) {
				_perfData->push_back(frob_norm(_x(n1&0) - _b(n1&0)));
				LOG(ALS, "Calculated residual for perfData after tensor "<< currIndex <<": " << _perfData->back() << " ( " << std::abs(_perfData->back() - energy) << " vs " << _convergenceEpsilon << " ) ");
				if (printProgress) {
					std::cout << "optimized tensor "<< currIndex << ": " << _perfData->back() << " ( \t" << std::abs(_perfData->back() - energy) << " vs \t" << _convergenceEpsilon << " ) \r" << std::flush;
				}
			}
			
			// Are we done with the sweep?
			if ((!walkingRight && currIndex==0) || (walkingRight && currIndex==d-1)) {
				LOG(ALS, "Sweep Done");
				halfSweepCount += 1;
				
				lastEnergy2 = lastEnergy;
				lastEnergy = energy;
				if (_perfData != nullptr) {
					energy = _perfData->back(); // Energy already calculated
				} else {
					LOG(ALS, "Calculating energy for loop condition");
					energy = frob_norm(_x(n1&0) - _b(n1&0));
				}
				
				// Conditions for loop termination
				if (halfSweepCount == _numHalfSweeps || std::abs(lastEnergy-energy) < _convergenceEpsilon || std::abs(lastEnergy2-energy) < _convergenceEpsilon || d<=1) {
					// we are done! yay
					LOG(ALS, "ALS done, " << energy << " " << lastEnergy << " " << std::abs(lastEnergy2-energy) << " " << std::abs(lastEnergy-energy) << " < " << _convergenceEpsilon);
					if (preserveCorePosition && cannoAtTheEnd) {
						_x.move_core(corePosAtTheEnd);
					}
					return energy; 
				}
				
				walkingRight = !walkingRight;
				LOG(ALS, "Start sweep " << (walkingRight?"right":"left"));
			}
			
			
			if (walkingRight) {
				// Move core to next position
				_x.move_core(currIndex+1);
				
				// Move one site to the right
				bxR.pop_back();
				
				LOG(ALS, "Calc BxL");
				tmp(r1,r2) = bxL.back()(cr1,cr2) * _b.get_component(currIndex)(cr1, n1, r1) * _x.get_component(currIndex)(cr2, n1, r2);
				bxL.emplace_back(std::move(tmp));
				currIndex++;;
			} else {
				// Move core to next position
				_x.move_core(currIndex-1);
			
				// move one site to the left
				bxL.pop_back();
				
				LOG(ALS, "Calc bxR");
				tmp(r1,r2) = bxR.back()(cr1,cr2) * _b.get_component(currIndex)(r1, n1, cr1) * _x.get_component(currIndex)(r2, n1, cr2);
				bxR.emplace_back(std::move(tmp));
				currIndex--;
			}
		}
	}
}
