// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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

#include <xerus/algorithms/xals.h>
#include <xerus/basic.h>
#include <xerus/misc/internal.h>
 
#include <xerus/indexedTensorMoveable.h>
#include <xerus/misc/basicArraySupport.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace xerus {
    
     class InternalSolver {
        const size_t d;
        
         std::vector<Tensor> leftStack;
         std::vector<Tensor> rightStack;
         
         std::vector<Tensor> leftBStack;
         std::vector<Tensor> rightBStack;
         
         
        TTTensor& x;
        const TTOperator& A;
        const TTTensor& b;
        
        Tensor leftStackA,rightStackA;
        Tensor leftBStackB,rightBStackB;
         
     public:
         InternalSolver(TTTensor& _x, const TTOperator& _A, const TTTensor& _b) : d(_x.degree()),leftStack(d), rightStack(d), x(_x), A(_A), b(_b) {
             
         }
         
         void calc_left_stack(const size_t _position) {
             const Tensor& xi = x.get_component(_position);
             
             if(_position == 0) {
                 Tensor tmp;
                 Tensor A0 = A.get_component(_position);
                 A0.reinterpret_dimensions({A0.dimensions[1], A0.dimensions[2], A0.dimensions[3]});
                 contract(tmp, A0, true, xi, false, 1);
                 contract(leftStack[_position], A0, true, xi, false, 1);
             } else {
                 contract(leftStack[_position], leftStack[_position-1], x.component(_position), 1);
             }
         }
         
         
         void calc_right_stack(const size_t _position) {
             
         }
         
         double calc_residual_norm() {
             Index i,j;
             Tensor tmp;
             tmp(i&0) = A(i/2, j/2)*x(j&0);
             return frob_norm(tmp);
         }
         
         void solve() {
             const double solutionsNorm = frob_norm(b);
             std::vector<double> residuals(10, 1000.0);
			const size_t maxIterations = 1000;
			
			for(size_t iteration = 0; maxIterations == 0 || iteration < maxIterations; ++iteration) {
				x.move_core(0, true);
				
				// Rebuild right stack
				for(size_t corePosition = d-1; corePosition > 0; --corePosition) {
					calc_right_stack(corePosition);
				}
				
                residuals.push_back(calc_residual_norm()/solutionsNorm);
                if(residuals.back()/residuals[residuals.size()-10] > 0.99) {
                    LOG(ALS, "Residual decrease from " << std::scientific << residuals[10] << " to " << std::scientific << residuals.back() << " in " << residuals.size()-10 << " iterations.");
                    return; // We are done!
                }
				
				// Sweep Right -> Left
				for(size_t corePosition = 0; corePosition < x.degree(); ++corePosition) {
                    // Actual Work to do
					Tensor op, rhs;
                    if(corePosition == 0) {
                        leftStackA =  A.get_component(corePosition);
                        leftStackA.reinterpret_dimensions({leftStackA.dimensions[1], leftStackA.dimensions[2], leftStackA.dimensions[3]});
                        
                        leftBStackB =  B.get_component(corePosition);
                        leftBStackB.reinterpret_dimensions({leftBStackB.dimensions[1], leftBStackB.dimensions[2]});
                    } else {
                        contract(leftStackA, leftStack[corePosition-1],  A.get_component(corePosition), 1);
                        contract(leftBStackB, leftBStack[corePosition-1],  b.get_component(corePosition), 1);
                    }
                    contract(op, leftStackA, rightStack[corePosition+1], 1);
                    contract(rhs, leftBStackB, rightBStack[corePosition+1], 1);
                    
                    solve_least_squares(x.component(corePosition), op, rhs, 0);
                    
					// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
					if(corePosition+1 < d) {
						x.move_core(corePosition+1, true);
						calc_left_stack(corePosition);
					}
				}
				
				// Sweep Right -> Left
				for(size_t corePosition = d-1; corePosition > 0; --corePosition) {
                    //Actual Work to do
                    
					// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
					if(corePosition > 0) {
						x.move_core(corePosition-1, true);
						calc_right_stack(corePosition);
					}
				}
				
			}
         }
         
     };

        void xals(TTTensor& _x, const TTOperator& _A, const TTTensor& _b)  {
            InternalSolver solver(_x, _A, _b);
            return solver.solve();
        }
}

