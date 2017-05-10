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
		const Index i1, i2, i3, j1 , j2, j3, k1, k2;
		
		const size_t d;
		
		std::vector<Tensor> leftAStack;
		std::vector<Tensor> rightAStack;
		
		std::vector<Tensor> leftBStack;
		std::vector<Tensor> rightBStack;
		
		TTTensor& x;
		const TTOperator& A;
		const TTTensor& b;
		
	public:
		
		InternalSolver(TTTensor& _x, const TTOperator& _A, const TTTensor& _b) : d(_x.degree()), leftAStack(d), rightAStack(d), leftBStack(d), rightBStack(d), x(_x), A(_A), b(_b) { }
		
		
		void calc_left_stack(const size_t _position) {
			Tensor xi = x.get_component(_position);
			Tensor Ai = A.get_component(_position);
			Tensor bi = b.get_component(_position);
			
			if(_position == 0) {
				xi.reinterpret_dimensions({xi.dimensions[1], xi.dimensions[2]});
				Ai.reinterpret_dimensions({Ai.dimensions[1], Ai.dimensions[2], Ai.dimensions[3]});
				bi.reinterpret_dimensions({bi.dimensions[1], bi.dimensions[2]});
				
				leftAStack[_position](i1, i2, i3) = xi(k1, i1)*Ai(k1, k2, i2)*xi(k2, i3);
				leftBStack[_position](i1, i2) = xi(k1, i1)*bi(k1, i2);
			} else {
				leftAStack[_position](i1, i2, i3) = leftAStack[_position-1](j1, j2, j3)*xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3);
				leftBStack[_position](i1, i2) = leftBStack[_position-1](j1, j2)*xi(j1, k1, i1)*bi(j2, k1, i2);
			}
		}
		
		
		void calc_right_stack(const size_t _position) {
			Tensor xi = x.get_component(_position);
			Tensor Ai = A.get_component(_position);
			Tensor bi = b.get_component(_position);
			
			if(_position == d-1) {
				xi.reinterpret_dimensions({xi.dimensions[0], xi.dimensions[1]});
				Ai.reinterpret_dimensions({Ai.dimensions[0], Ai.dimensions[1], Ai.dimensions[2]});
				bi.reinterpret_dimensions({bi.dimensions[0], bi.dimensions[1]});
				
				rightAStack[_position](i1, i2, i3) = xi(i1, k1)*Ai(i2, k1, k2)*xi(i3, k2);
				rightBStack[_position](i1, i2) = xi(i1, k1)*bi(i2, k1);
			} else {
				rightAStack[_position](i1, i2, i3) = xi(i1, k1, j1)*Ai(i2, k1, k2, j2)*xi(i3, k2, j3)*rightAStack[_position+1](j1, j2, j3);
				rightBStack[_position](i1, i2) = xi(i1, k1, j1)*bi(i2, k1, j2)*rightBStack[_position+1](j1, j2);
			}
		}
		
		double calc_residual_norm() {
			const Index i, j;
			TTTensor tmp;
			tmp(i&0) = A(i/2, j/2)*x(j&0)-b(i&0);
			return frob_norm(tmp);
		}
		
		
		void solve() {
			const double solutionsNorm = frob_norm(b);
			std::vector<double> residuals(10, 1000.0);
			const size_t maxIterations = 1;
			
			
			// Rebuild right stack
			x.move_core(0, true);
			for(size_t corePosition = d-1; corePosition > 0; --corePosition) {
				calc_right_stack(corePosition);
			}
			
			for(size_t iteration = 0; maxIterations == 0 || iteration < maxIterations; ++iteration) {
				// Calculate residual and check end condition
				residuals.push_back(calc_residual_norm()/solutionsNorm);
				if(residuals.back()/residuals[residuals.size()-10] > 0.99) {
					LOG(xALS, "Residual decrease from " << std::scientific << residuals[10] << " to " << std::scientific << residuals.back() << " in " << residuals.size()-10 << " iterations.");
					return; // We are done!
				} else {
// 					LOG(xALS, "Residual decrease from " << std::scientific << residuals[10] << " to " << std::scientific << residuals.back() << " in " << residuals.size()-10 << " iterations.");
				}
				
				
				// Sweep Left -> Right
				for(size_t corePosition = 0; corePosition < d; ++corePosition) {
					Tensor op, rhs;
					
					Tensor Ai = A.get_component(corePosition);
					Tensor bi = b.get_component(corePosition);
					
					if(corePosition == 0) {
						Ai.reinterpret_dimensions({Ai.dimensions[1], Ai.dimensions[2], Ai.dimensions[3]});
						bi.reinterpret_dimensions({bi.dimensions[1], bi.dimensions[2]});
						
						op(i2, i3, j2, j3) = Ai(i2, j2, k2)*rightAStack[corePosition+1](i3, k2, j3);
						rhs(i2, i3) = bi(i2, k2)*rightBStack[corePosition+1](i3, k2);
					} else if(corePosition == d-1) {
						Ai.reinterpret_dimensions({Ai.dimensions[0], Ai.dimensions[1], Ai.dimensions[2]});
						bi.reinterpret_dimensions({bi.dimensions[0], bi.dimensions[1]});
						
						op(i1, i2, j1, j2) = leftAStack[corePosition-1](i1, k1, j1)*Ai(k1, i2, j2);
						rhs(i1, i2) = leftBStack[corePosition-1](i1, k1)*bi(k1, i2);
					} else {
						op(i1, i2, i3, j1, j2, j3) = leftAStack[corePosition-1](i1, k1, j1)*Ai(k1, i2, j2, k2)*rightAStack[corePosition+1](i3, k2, j3);
						rhs(i1, i2, i3) = leftBStack[corePosition-1](i1, k1)*bi(k1, i2, k2)*rightBStack[corePosition+1](i3, k2);
					}
					
					solve_least_squares(x.component(corePosition), op, rhs, 0);
					
					if(corePosition == 0) {
						x.component(corePosition).reinterpret_dimensions({1, x.component(corePosition).dimensions[0], x.component(corePosition).dimensions[1]});
					} else if(corePosition == d-1) {
						x.component(corePosition).reinterpret_dimensions({x.component(corePosition).dimensions[0], x.component(corePosition).dimensions[1], 1});
					}
					
					// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
					if(corePosition+1 < d) {
						x.move_core(corePosition+1, true);
						calc_left_stack(corePosition);
					}
				}
				
				
				// Sweep Right -> Left
				for(size_t corePosition = d-1; corePosition > 0; --corePosition) {
// 					Tensor op, rhs;
// 					
// 					Tensor Ai = A.get_component(corePosition);
// 					Tensor bi = b.get_component(corePosition);
// 					
// 					if(corePosition == 0) {
// 						Ai.reinterpret_dimensions({Ai.dimensions[1], Ai.dimensions[2], Ai.dimensions[3]});
// 						bi.reinterpret_dimensions({bi.dimensions[1], bi.dimensions[2]});
// 						
// 						op(i2, i3, j2, j3) = Ai(i2, j2, k2)*rightAStack[corePosition+1](i3, k2, j3);
// 						rhs(i2, i3) = bi(i2, k2)*rightBStack[corePosition+1](i3, k2);
// 					} else if(corePosition == d-1) {
// 						Ai.reinterpret_dimensions({Ai.dimensions[0], Ai.dimensions[1], Ai.dimensions[2]});
// 						bi.reinterpret_dimensions({bi.dimensions[0], bi.dimensions[1]});
// 						
// 						op(i1, i2, j1, j2) = leftAStack[corePosition-1](i1, k1, j1)*Ai(k1, i2, j2);
// 						rhs(i1, i2) = leftBStack[corePosition-1](i1, k1)*bi(k1, i2);
// 					} else {
// 						op(i1, i2, i3, j1, j2, j3) = leftAStack[corePosition-1](i1, k1, j1)*Ai(k1, i2, j2, k2)*rightAStack[corePosition+1](i3, k2, j3);
// 						rhs(i1, i2, i3) = leftBStack[corePosition-1](i1, k1)*bi(k1, i2, k2)*rightBStack[corePosition+1](i3, k2);
// 					}
// 					
// 					solve_least_squares(x.component(corePosition), op, rhs, 0);
// 					
// 					if(corePosition == 0) {
// 						x.component(corePosition).reinterpret_dimensions({1, x.component(corePosition).dimensions[0], x.component(corePosition).dimensions[1]});
// 					} else if(corePosition == d-1) {
// 						x.component(corePosition).reinterpret_dimensions({x.component(corePosition).dimensions[0], x.component(corePosition).dimensions[1], 1});
// 					}
					
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

