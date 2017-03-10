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
//  * @brief Header file for the ADF algorithm and its variants.
//  */
// 
// #pragma once
// 
// #include "../index.h"
// #include "../ttNetwork.h"
// #include "../performanceData.h"
// #include "../measurments.h"
// 
// namespace xerus {
// 	
// 	
// 	
// 	class UQADF {
// 		///@brief Indices for all internal functions.
// 		const Index r1, r2, i1;
// 		
// 		///@brief Reference to the current solution (external ownership)
// 		TTTensor& x;
// 		
// 		///@brief Degree of the solution.
// 		const size_t degree;
// 		
// 		///@brief Maximally allowed ranks.
// 		const std::vector<size_t> maxRanks;
// 		
// 		///@brief Number of measurments (i.e. measurments.size())
// 		const size_t numMeasurments;
// 		
// 		///@brief The two norm of the measured values
// 		const value_t normMeasuredValues;
// 		
// 		///@brief Maximal allowed number of iterations (one iteration = one sweep)
// 		const size_t maxIterations; 
// 		
// 		///@brief The target residual norm at which the algorithm shall stop.
// 		const double targetResidualNorm;
// 		
// 		///@brief Minimal relative decrease of the residual norm ( (oldRes-newRes)/oldRes ) until either the ranks are increased (if allowed) or the algorithm stops.
// 		const double minimalResidualNormDecrease; 
// 		
// 		///@brief The current iteration.
// 		size_t iteration;
// 		
// 		///@brief Current residual norm. Updated at the beginning of each iteration.
// 		double residualNorm;
// 		
// 		///@brief The residual norm of the last iteration.
// 		double lastResidualNorm;
// 		
// 		///@brief The current residual, saved as vector (instead of a order one tensor).
// 		std::vector<value_t> residual;
// 		
// 		///@brief The current projected Gradient component. That is E(A^T(Ax-b))
// 		Tensor projectedGradientComponent;
// 		
// 		///@brief Ownership holder for a (degree+2)*numMeasurments array of Tensor pointers. (Not used directly)
// 		std::unique_ptr<Tensor*[]> forwardStackMem;
// 		
// 		/** @brief Array [numMeasurments][degree]. For positions smaller than the current corePosition and for each measurment, this array contains the pre-computed
// 		 * contraction of the first _ component tensors and the first _ components of the measurment operator. These tensors are deduplicated in the sense that for each unqiue
// 		 * part of the position only one tensor is actually stored, which is why the is an array of pointers. The Tensors at the current corePosition are used as 
// 		 * scatch space. For convinience the underlying array (forwardStackMem) is larger, wherefore also the positions -1 and degree are allow, all poining to a {1} tensor
// 		 * containing 1 as only entry. Note that position degree-1 must not be used.
// 		 **/
// 		Tensor* const * const forwardStack;
// 		
// 		/// @brief Ownership holder for the unqiue Tensors referenced in forwardStack.
// 		std::unique_ptr<Tensor[]> forwardStackSaveSlots;
// 		
// 		/// @brief Vector containing for each corePosition a vector of the smallest ids of each group of unique forwardStack entries.
// 		std::vector<std::vector<size_t>> forwardUpdates;
// 		
// 		
// 		///@brief Ownership holder for a (degree+2)*numMeasurments array of Tensor pointers. (Not used directly)
// 		std::unique_ptr<Tensor*[]> backwardStackMem;
// 		
// 		/** @brief Array [numMeasurments][degree]. For positions larger than the current corePosition and for each measurment, this array contains the pre-computed
// 		 * contraction of the last _ component tensors and the last _ components of the measurment operator. These tensors are deduplicated in the sense that for each unqiue
// 		 * part of the position only one tensor is actually stored, which is why the is an array of pointers. The Tensors at the current corePosition are used as 
// 		 * scratch space. For convinience the underlying array (forwardStackMem) is larger, wherefore also the positions -1 and degree are allow, all poining to a {1} tensor
// 		 * containing 1 as only entry. Note that position zero must not be used.
// 		 **/
// 		Tensor* const * const backwardStack;
// 		
// 		/// @brief Ownership holder for the unqiue Tensors referenced in backwardStack.
// 		std::unique_ptr<Tensor[]> backwardStackSaveSlots;
// 		
// 		/// @brief Vector containing for each corePosition a vector of the smallest ids of each group of unique backwardStack entries.
// 		std::vector<std::vector<size_t>> backwardUpdates;
// 		
// 		/// @brief: Norm of each rank one measurment operator
// 		std::unique_ptr<double[]> measurmentNorms;
// 		
// 		///@brief: Reference to the performanceData object (external ownership)
// 		PerformanceData& perfData;
// 		
// 		void solve(TTTensor _x, const std::vector<std::vector<double>> _randomVariables, std::vector<Tensor> _solutions);
// 	};
// 	
// }
