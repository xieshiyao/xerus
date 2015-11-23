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
 * @brief Header file for the ADF algorithm and its variants.
 */

#pragma once

#include "../ttNetwork.h"
#include "../performanceData.h"
#include "../measurments.h"

namespace xerus {
	/**
	 * @brief Wrapper class for all ADF variants.
	 * @details By creating a new object of this class and modifying the member variables, the behaviour of the solver can be modified.
	 */
    class ADFVariant {
		protected:
			
		template<class MeasurmentSet> 
		class InternalSolver {
		protected:
			const Index r1, r2, i1;
			
			TTTensor& x;
			const size_t degree;
			const std::vector<size_t> maxRanks;
			
			const MeasurmentSet& measurments;
			const size_t numMeasurments;
			const value_t normMeasuredValues;
			
			const size_t maxIterations; 
			const double targetResidualNorm;
			const double minimalResidualNormDecrese; 
			
			size_t iteration;
			double residualNorm;
			double lastResidualNorm;
			
			std::vector<value_t> residual;
			FullTensor projectedGradientComponent;
			
			std::unique_ptr<FullTensor*[]> forwardStackMem;
			FullTensor* const * const forwardStack;
			std::unique_ptr<FullTensor[]> forwardStackSaveSlots;
			std::vector<std::vector<size_t>> forwardUpdates;
			
			std::unique_ptr<FullTensor*[]> backwardStackMem;
			FullTensor* const * const backwardStack;
			std::unique_ptr<FullTensor[]> backwardStackSaveSlots;
			std::vector<std::vector<size_t>> backwardUpdates;
			
			PerformanceData& perfData;
			
			static double calculate_norm_of_measured_values(const MeasurmentSet& _measurments);
			
			void construct_stacks(std::unique_ptr< xerus::FullTensor[] >& _stackSaveSlot, std::vector< std::vector< size_t > >& _updates, const std::unique_ptr<FullTensor*[]>& _stackMem, const bool _forward);
			
			void resize_stack_tensors();
			
			std::vector<FullTensor> get_fixed_components(const Tensor& _component);
			
			void update_backward_stack(const size_t _corePosition, const Tensor& _currentComponent);
			
			void update_forward_stack(const size_t _corePosition, const Tensor& _currentComponent);
			
			/**
			* @brief: Calculates the current residual, i.e. b-Ax
			*/
			void calculate_residual( const size_t _corePosition );
			
			/**
			* @brief: Calculates the (component at _corePosition of the) projected gradient from the residual, i.e. E(A^T(b-Ax)), where E is the projection on the current component
			*/
			void calculate_projected_gradient( const size_t _corePosition);
			
			/**
			* @brief: Calculates ||P_n (A(E(A^T(b-Ax)))))|| = ||P_n (A(E(A^T(residual)))))|| =  ||P_n (A(E(gradient)))|| for each n, 
			* where P_n sets all entries equals zero except where the index at _corePosition is equals n. In case of RankOneMeasurments,
			* the calculation is not slicewise (only n=0 is set).
			*/
			std::vector<value_t> calculate_slicewise_norm_A_projGrad( const size_t _corePosition);
			
			void update_x(const std::vector<value_t>& _normAProjGrad, const size_t _corePosition);
			
			void solve_with_current_ranks();
			
		public:
			double solve();
			
			InternalSolver(TTTensor& _x, 
							const std::vector<size_t>& _maxRanks, 
							const MeasurmentSet& _measurments, 
							const size_t _maxIteration, 
							const double _targetResidualNorm, 
							const double _minimalResidualNormDecrese, 
							PerformanceData& _perfData ) : 
				x(_x),
				degree(_x.degree()),
				maxRanks(_maxRanks),
				
				measurments(_measurments),
				numMeasurments(_measurments.size()),
				normMeasuredValues(calculate_norm_of_measured_values(_measurments)),
				
				maxIterations(_maxIteration),
				targetResidualNorm(_targetResidualNorm),
				minimalResidualNormDecrese(_minimalResidualNormDecrese),
				
				iteration(0),
				residualNorm(std::numeric_limits<double>::max()), 
				lastResidualNorm(std::numeric_limits<double>::max()),
				
				residual(numMeasurments),
				
				forwardStackMem(new FullTensor*[numMeasurments*(degree+2)]),
				forwardStack(forwardStackMem.get()+numMeasurments),
				forwardUpdates(degree),
					
				backwardStackMem(new FullTensor*[numMeasurments*(degree+2)]),
				backwardStack(backwardStackMem.get()+numMeasurments),
				backwardUpdates(degree),
				
				perfData(_perfData) 
				{
					REQUIRE(_x.is_valid_tt(), "_x must be a valid TT-Tensor.");
					REQUIRE(numMeasurments > 0, "Need at very least one measurment.");
					REQUIRE(measurments.degree() == degree, "Measurment degree must coincide with x degree.");
				}
				
		};
		
	public:
        size_t maxIterations; ///< Maximum number of sweeps to perform. Set to 0 for infinite.
        double targetResidualNorm; ///< Target residual. The algorithm will stop upon reaching a residual smaller than this value.
        double minimalResidualNormDecrese; // The minimal relative decrese of the residual per step  ( i.e. (lastResidual-residual)/lastResidual ). If the avg. of the last three steps is smaller than this value, the algorithm stops.
        bool printProgress; ///< informs the user about the current progress via std::cout (one continuously overwritten line)
        
		/// fully defining constructor. alternatively ALSVariants can be created by copying a predefined variant and modifying it
        ADFVariant(const size_t _maxIteration, const double _targetResidual, const double _minimalResidualDecrese, const bool _printProgress) 
                : maxIterations(_maxIteration), targetResidualNorm(_targetResidual), minimalResidualNormDecrese(_minimalResidualDecrese), printProgress(_printProgress) { }
        
        /**
		* @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		* @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		* @param _measurments the available measurments.
		* @param _perfData optinal performanceData object to be used.
		* @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		*/
        double operator()(TTTensor& _x, const SinglePointMeasurmentSet& _measurments, PerformanceData& _perfData = NoPerfData) const;
		
		/**
		* @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		* @param[in,out] _x On input: an initial guess of the solution, may be of smaller rank. On output: The reconstruction found by the algorithm.
		* @param _measurments the available measurments.
		* @param _maxRanks the maximal ranks the algorithm may use to decrease the resdiual.
		* @param _perfData optinal performanceData object to be used.
		* @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		*/
        double operator()(TTTensor& _x, const SinglePointMeasurmentSet& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData = NoPerfData) const;
		
		/**
		* @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		* @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		* @param _b the available measurments.
		* @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		*/
        double operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, PerformanceData& _perfData = NoPerfData) const;
	
		/**
		* @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		* @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		* @param _b the available measurments.
		* @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		*/
        double operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData = NoPerfData) const;
    };
	
	/// @brief Default variant of the ADF algorithm
    const ADFVariant ADF(0, 1e-8, 5e-4, true);
}

