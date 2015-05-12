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


#include "../xerus.h"

namespace xerus {

    void ALSVariant::lapack_solver(const TensorNetwork &_A, Tensor &_x, const Tensor &_b) {
        const size_t d = _x.degree();
        FullTensor A(_A);
        Index i,j;
        _x(i^d) = _b(j^d) / A(j^d, i^d);
        REQUIRE(d <= 3, "dmrg not yet implemented in lapack_solver");// TODO split result into d-2 tensors -> choose correct tensor as core!
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
        
        auto A = [&](size_t _i)->const Tensor& {
            return *_A.nodes[_i].tensorObject;
        };
        auto X = [&](size_t _i)->Tensor& {
            return *_x.nodes[_i].tensorObject;
        };
        auto B = [&](size_t _i)->const Tensor& {
            return *_b.nodes[_i].tensorObject;
        };
        
        std::vector<FullTensor> xAxL, xAxR;
        std::vector<FullTensor> bxL, bxR;
        
        _x.cannonicalize_left();

        // Create stack of contracted network-parts (from right to left)
        FullTensor tmpA, tmpB;
        for (size_t i = _x.degree()-1; i > sites-1; --i) {
            if (xAxR.empty()) {
                tmpA(r1, r2, r3) = X(i)(r1, n1) * A(i)(r2, n1, n2) * X(i)(r3, n2);
                tmpB(r1, r2) = B(i)(r1, n1) * X(i)(r2, n1);
            } else {
                tmpA(r1, r2, r3) = tmpA(cr1, cr2, cr3) * X(i)(r1, n1, cr1) * A(i)(r2, n1, n2, cr2) * X(i)(r3, n2, cr3);
                tmpB(r1, r2) = tmpB(cr1, cr2) * B(i)(r1, n1, cr1) * X(i)(r2, n1, cr2);
            }
            xAxR.emplace_back(std::move(tmpA));
            bxR.emplace_back(std::move(tmpB));
        }
        
        if (_perfData != nullptr) {
            // calculate current energy 0.5*xAx - bx
    // 		FullTensor xax, bx;
    // 		xax() = xAxL.back()(r1,r2,r3) * X[0](r1,n1,cr1) * A[0](r2,n1,n2,cr2) * X[0](r3,n2,cr3) * xAxR.back()(cr1,cr2,cr3);
    // 		bx() = bxL.back()(r1,r2) * B[0](r1,n1,cr1) * X[0](r2,n1,cr2) * bxR.back()(cr1,cr2);
    // 		_perfData->push_back(0.5*xax.data[0] - bx.data[0]);
            _perfData->push_back(frob_norm(_A(n1^d,n2^d)*_x(n2^d) - _b(n1^d)));
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
            
            // solve Atilde * Xtilde = Btilde
            size_t localD;
            if (xAxL.empty()) {
                if (xAxR.empty()) {
                    localD = 1;
                    ATilde(n1, n2) = A(currIndex)(n1, n2);
                    BTilde(n1) = B(currIndex)(n1);
                } else {
                    localD = 2;
                    ATilde(n1,cr1, n2,cr3) = A(currIndex)(n1, n2, cr2) * xAxR.back()(cr1, cr2, cr3);
                    BTilde(n1,cr2) = B(currIndex)(n1, cr1) * bxR.back()(cr1,cr2);
                }
            } else { // xAxL not empty
                if (xAxR.empty()) {
                    localD = 2;
                    ATilde(r1,n1, r3,n2) = xAxL.back()(r1,r2,r3) * A(currIndex)(r2, n1, n2);
                    BTilde(r2,n1) = bxL.back()(r1,r2) * B(currIndex)(r1, n1);
                } else { // default case:
                    localD = 3;
                    ATilde(r1,n1,cr1, r3,n2,cr3) = xAxL.back()(r1,r2,r3) * A(currIndex)(r2, n1, n2, cr2) * xAxR.back()(cr1, cr2, cr3);
                    BTilde(r2,n1,cr2) = bxL.back()(r1,r2) * B(currIndex)(r1, n1, cr1) * bxR.back()(cr1,cr2);
                }
            }
            
            // only change component tensor if the local residual is large enough
            if (minimumLocalResidual <= 0 || frob_norm(IndexedTensorMoveable<Tensor>(ATilde(r1^localD, r2^localD)*X(currIndex)(r2^localD)) - BTilde(r1^localD)) > minimumLocalResidual) {
                localSolver(ATilde, X(currIndex), BTilde);
                changedSmth = true;
            }
            
            if (_perfData != nullptr) {
                // calculate current energy 0.5*xAx - bx
    // 			FullTensor xax, bx;
    // 			xax() = xAxL.back()(r1,r2,r3) * X[currIndex](r1,n1,cr1) * A[currIndex](r2,n1,n2,cr2) * X[currIndex](r3,n2,cr3) * xAxR.back()(cr1,cr2,cr3);
    // 			bx() = bxL.back()(r1,r2) * B[currIndex](r1,n1,cr1) * X[currIndex](r2,n1,cr2) * bxR.back()(cr1,cr2);
    // 			_perfData->push_back(0.5*xax.data[0] - bx.data[0]);
                _perfData->push_back(frob_norm(_A(n1^d,n2^d)*_x(n2^d) - _b(n1^d)));
                LOG(ALS, "calculated residual for perfData after tensor "<< currIndex <<": " << _perfData->back() << " ( " << std::abs(_perfData->back() - energy) << " vs " << _convergenceEpsilon << " ) ");
                if (printProgress) {
                    std::cout << "optimized tensor "<< currIndex << ": " << _perfData->back() << " ( \t" << std::abs(_perfData->back() - energy) << " vs \t" << _convergenceEpsilon << " ) \r" << std::flush;
                }
            }
            
            // are we done with the sweep?
            if ((!walkingRight && currIndex==0) || (walkingRight && currIndex==d-sites)) {
                LOG(ALS, "Sweep Done");
                
                lastEnergy2 = lastEnergy;
                lastEnergy = energy;
                if (_perfData != nullptr) {
                    // energy already calculated
                    energy = _perfData->back();
                } else {
                    LOG(ALS, "calculating energy for loop condition");
//     				FullTensor xax(0), bx(0);
//     				xax() = xAxL.back()(r1,r2,r3) * X(currIndex)(r1,n1,cr1) * A(currIndex)(r2,n1,n2,cr2) * X(currIndex)(r3,n2,cr3) * xAxR.back()(cr1,cr2,cr3);
//     				bx() = bxL.back()(r1,r2) * B(currIndex)(r1,n1,cr1) * X(currIndex)(r2,n1,cr2) * bxR.back()(cr1,cr2);
//     				energy = 0.5*xax[0] - bx[0];
                    energy = frob_norm(_A(n1^d,n2^d)*_x(n2^d) - _b(n1^d));
                }
                
                LOG(ALS, "stats: " << xAxL.size() << " " << xAxR.size() << " " << bxL.size() << " " << bxR.size() << " energy: " << energy << " deltas: " << (1-lastEnergy/energy) << " " << (1-lastEnergy2/energy));
                
                // conditions for loop termination
                if (!changedSmth || std::abs(lastEnergy-energy) < _convergenceEpsilon || std::abs(lastEnergy2-energy) < _convergenceEpsilon || d<=sites) {
                    // we are done! yay
                    LOG(ALS, "ALS done, " << energy << " " << lastEnergy << " " << std::abs(lastEnergy2-energy) << " " << std::abs(lastEnergy-energy) << " < " << _convergenceEpsilon);
                    if (!walkingRight) {
                        _x.cannonicalize_right();
                    }
                    return energy; 
                }
//                 if (energy > lastEnergy && energy > lastEnergy2) {
//                     LOG(warning, "ALS did not converge to desired accuracy. curr:"<< energy <<" last: " << (lastEnergy) << " " << (lastEnergy2) << " eps = " << _convergenceEpsilon);
//                     if (!walkingRight) {
//                         _x.cannonicalize_right();
//                     }
//                     return energy;
//                 }
                
                walkingRight = !walkingRight;
                changedSmth = false;
                LOG(ALS, "Start sweep " << (walkingRight?"right":"left"));
            }
            
            if (walkingRight) {
                // move one site to the right
                xAxR.pop_back();
                bxR.pop_back();
                
                if (sites == 1) {
                    REQUIRE(currIndex + 1 < d, "Internal error");
                    // move core
                    FullTensor core;
                    (X(currIndex)(r1&2,r2,cr3), core(cr3,r3)) = QR(X(currIndex)(r1&2,r2,r3));
                    X(currIndex+1)(cr1, n1, r1&2) = core(cr1,cr2) * X(currIndex+1)(cr2, n1, r1&2);
                }
                
                LOG(ALS, "Calc xAxL");
                if (xAxL.empty()) {
                    tmpA(r1,r2,r3) = X(currIndex)(n1,r1) * A(currIndex)( n1, n2, r2) * X(currIndex)(n2,r3);
                } else {
                    tmpA(r1,r2,r3) = xAxL.back()(cr1,cr2,cr3) * X(currIndex)(cr1,n1,r1) * A(currIndex)(cr2, n1, n2, r2) * X(currIndex)(cr3,n2,r3);
                }
                xAxL.push_back(tmpA);
                
                LOG(ALS, "Calc BxL");
                if (bxL.empty()) {
                    tmpB(r1,r2) = B(currIndex)(n1, r1) * X(currIndex)(n1, r2);
                } else {
                    tmpB(r1,r2) = bxL.back()(cr1,cr2) * B(currIndex)(cr1, n1, r1) * X(currIndex)(cr2, n1, r2);
                }
                bxL.push_back(tmpB);
                currIndex += 1;
            } else {
                // move one site to the left
                xAxL.pop_back();
                bxL.pop_back();
                
                if (sites == 1) {
                    // move core
                    FullTensor core(2);
                    (core(r1,cr1), X(currIndex)(cr1,r2,r3&2)) = RQ(X(currIndex)(r1,r2,r3&2));
                    REQUIRE(currIndex >= 1, "internal error");
                    X(currIndex-1)(cr1&2, n1, r1) = X(currIndex-1)(cr1&2, n1, r2) * core(r2,r1);
                }
                
                LOG(ALS, "Calc xAxR");
                if (xAxR.empty()) {
                    tmpA(r1,r2,r3) = X(currIndex)(r1,n1) * A(currIndex)(r2, n1, n2) * X(currIndex)(r3,n2);
                } else {
                    tmpA(r1,r2,r3) = xAxR.back()(cr1,cr2,cr3) * X(currIndex)(r1,n1,cr1) * A(currIndex)(r2, n1, n2, cr2) * X(currIndex)(r3,n2,cr3);
                }
                xAxR.push_back(tmpA);
                
                LOG(ALS, "Calc bxR");
                if (bxR.empty()) {
                    tmpB(r1,r2) = B(currIndex)(r1, n1) * X(currIndex)(r2, n1);
                } else {
                    tmpB(r1,r2) = bxR.back()(cr1,cr2) * B(currIndex)(r1, n1, cr1) * X(currIndex)(r2, n1, cr2);
                }
                bxR.push_back(tmpB);
                currIndex -= 1;
            }
        }
    }

}
