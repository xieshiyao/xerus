// // Xerus - A General Purpose Tensor Library
// // Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
// // 
// // Xerus is free software: you can redistribute it and/or modify
// // it under the terms of the GNU Affero General Public License as published
// // by the Free Software Foundation, either version 3 of the License,
// // or (at your option) any later version.
// // 
// // Xerus is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// // GNU Affero General Public License for more details.
// // 
// // You should have received a copy of the GNU Affero General Public License
// // along with Xerus. If not, see <http://www.gnu.org/licenses/>.
// //
// // For further information on Xerus visit https://libXerus.org 
// // or contact us at contact@libXerus.org.
// 
// /**
//  * @file
//  * @brief Implementation of the ADF variants. 
//  */
// 
// #include <xerus/algorithms/uqAdf.h>
//  
// #include <xerus/indexedTensorMoveable.h>
// #include <xerus/misc/basicArraySupport.h>
// #include <xerus/misc/internal.h>
// 
// #ifdef _OPENMP
// 	#include <omp.h>
// #endif
// 
// namespace xerus {
// 	void UQADF::solve(TTTensor _x, const std::vector<std::vector<double>> _randomVariables, std::vector<Tensor> _solutions) {
// 		double resDec1 = 0.0, resDec2 = 0.0, resDec3 = 0.0;
// 		
// 		for(; maxIterations == 0 || iteration < maxIterations; ++iteration) {
// 			
// 			// Move core back to position zero
// 			x.move_core(0, true);
// 			
// 			// Rebuild backwardStack
// 			for(size_t corePosition = x.degree()-1; corePosition > 0; --corePosition) {
// 				update_backward_stack(corePosition, x.get_component(corePosition));
// 			}
// 			
// 			calculate_residual(0);
// 			
// 			lastResidualNorm = residualNorm;
// 			double residualNormSqr = 0;
// 			
// 			#pragma omp parallel for schedule(static) reduction(+:residualNormSqr)
// 			for(size_t i = 0; i < numMeasurments; ++i) {
// 				residualNormSqr += misc::sqr(residual[i]);
// 			}
// 			residualNorm = std::sqrt(residualNormSqr)/normMeasuredValues;
// 			
// 			perfData.add(iteration, residualNorm, x, 0);
// 			
// 			// Check for termination criteria
// 			double resDec4 = resDec3; resDec3 = resDec2; resDec2 = resDec1;
// 			resDec1 = residualNorm/lastResidualNorm;
// 			// 			LOG(wup, resDec1*resDec2*resDec3*resDec4);
// 			if(residualNorm < targetResidualNorm || resDec1*resDec2*resDec3*resDec4 > misc::pow(minimalResidualNormDecrease, 4)) { break; }
// 			
// 			
// 			// Sweep from the first to the last component
// 			for(size_t corePosition = 0; corePosition < degree; ++corePosition) {
// 				if(corePosition > 0) { // For corePosition 0 this calculation is allready done in the calculation of the residual.
// 					calculate_residual(corePosition);
// 				}
// 				
// 				calculate_projected_gradient(corePosition);
// 				
// 				const std::vector<value_t> normAProjGrad = calculate_slicewise_norm_A_projGrad(corePosition);
// 				
// 				update_x(normAProjGrad, corePosition);
// 				
// 				// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
// 				if(corePosition+1 < degree) {
// 					x.move_core(corePosition+1, true);
// 					update_forward_stack(corePosition, x.get_component(corePosition));
// 				}
// 			}
// 		}
// 	}
// 	
// } // namespace xerus
