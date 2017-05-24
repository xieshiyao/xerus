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
* @brief Implementation of a simple ALS variant. 
*/

#include <xerus/algorithms/xals.h>
#include <xerus/basic.h>
#include <xerus/misc/internal.h>

#include <xerus/indexedTensorMoveable.h>
#include <xerus/misc/basicArraySupport.h>


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
		
		InternalSolver(TTTensor& _x, const TTOperator& _A, const TTTensor& _b) : d(_x.degree()), x(_x), A(_A), b(_b) { 
			leftAStack.emplace_back(Tensor::ones(std::vector<size_t>(d, 1ul)));
			rightAStack.emplace_back(Tensor::ones(std::vector<size_t>(d, 1ul)));
			leftBStack.emplace_back(Tensor::ones(std::vector<size_t>(d, 1ul)));
			rightBStack.emplace_back(Tensor::ones(std::vector<size_t>(d, 1ul)));
		}
		
		
		void push_left_stack(const size_t _position) {
			const Tensor &xi = x.get_component(_position);
			const Tensor &Ai = A.get_component(_position);
			const Tensor &bi = b.get_component(_position);
			
			Tensor tmpA, tmpB;
			tmpA(i1, i2, i3) = leftAStack.back()(j1, j2, j3)*xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3);
			leftAStack.emplace_back(std::move(tmpA));
			tmpB(i1, i2) = leftBStack.back()(j1, j2)*xi(j1, k1, i1)*bi(j2, k1, i2);
			leftBStack.emplace_back(std::move(tmpB));
		}
		
		
		void push_right_stack(const size_t _position) {
			const Tensor &xi = x.get_component(_position);
			const Tensor &Ai = A.get_component(_position);
			const Tensor &bi = b.get_component(_position);
			
			Tensor tmpA, tmpB;
			tmpA(i1, i2, i3) = xi(i1, k1, j1)*Ai(i2, k1, k2, j2)*xi(i3, k2, j3)*rightAStack.back()(j1, j2, j3);
			rightAStack.emplace_back(std::move(tmpA));
			tmpB(i1, i2) = xi(i1, k1, j1)*bi(i2, k1, j2)*rightBStack.back()(j1, j2);
			rightBStack.emplace_back(std::move(tmpB));
		}
		
		double calc_residual_norm() {
			TTTensor tmp;
			tmp(i1&0) = A(i1/2, j1/2)*x(j1&0)-b(i1&0);
			return frob_norm(tmp);
		}
		
		
		void solve() {
			const double solutionsNorm = frob_norm(b);
			std::vector<double> residuals(10, 1000.0);
			const size_t maxIterations = 1000;
			
			// Rebuild right stack
			x.move_core(0, true);
			for(size_t pos = d-1; pos > 0; --pos) {
				push_right_stack(pos);
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
					
					const Tensor &Ai = A.get_component(corePosition);
					const Tensor &bi = b.get_component(corePosition);
					
					op(i1, i2, i3, j1, j2, j3) = leftAStack.back()(i1, k1, j1)*Ai(k1, i2, j2, k2)*rightAStack.back()(i3, k2, j3);
					rhs(i1, i2, i3) =            leftBStack.back()(i1, k1) *   bi(k1, i2, k2) *   rightBStack.back()(i3, k2);
					
					xerus::solve(x.component(corePosition), op, rhs, 0);
					
					// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
					if(corePosition+1 < d) {
						x.move_core(corePosition+1, true);
						push_left_stack(corePosition);
						rightAStack.pop_back();
						rightBStack.pop_back();
					}
				}
				
				
				// Sweep Right -> Left : only move core and update stacks
				x.move_core(0, true);
				for(size_t corePosition = d-1; corePosition > 0; --corePosition) {
					push_right_stack(corePosition);
					leftAStack.pop_back();
					leftBStack.pop_back();
				}
				
			}
		}
		
	};

	void xals(TTTensor& _x, const TTOperator& _A, const TTTensor& _b)  {
		InternalSolver solver(_x, _A, _b);
		return solver.solve();
	}
}

