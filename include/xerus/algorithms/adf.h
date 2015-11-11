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
		void solve_with_current_ranks(TTTensor& _x, 
										size_t& _iteration,
										double& _residual,
										double& _lastResidual,
										const MeasurmentSet& _measurments,
										const value_t _normMeasuredValues,
										FullTensor* const * const _forwardStack, 
										const std::vector<std::vector<size_t>>& _forwardUpdates, 
										FullTensor* const * const _backwardStack,
										const std::vector<std::vector<size_t>>& _backwardUpdates, 
										PerformanceData& _perfData ) const;
										
		template<class MeasurmentSet>
		double solve(xerus::TTTensor& _x, const MeasurmentSet& _measurments, const std::vector< size_t >& _maxRanks, PerformanceData& _perfData) const;
		
		const Index r1, r2, i1, i2;
	public:
        size_t maxInterations; ///< Maximum number of sweeps to perform. Set to 0 for infinite.
        double targetResidual; ///< Target residual. The algorithm will stop upon reaching a residual smaller than this value.
        double minimalResidualDecrese; // The minimal relative decrese of the residual per step  ( i.e. (lastResidual-residual)/lastResidual ). If the avg. of the last three steps is smaller than this value, the algorithm stops.
        bool printProgress; ///< informs the user about the current progress via std::cout (one continuously overwritten line)
        
		/// fully defining constructor. alternatively ALSVariants can be created by copying a predefined variant and modifying it
        ADFVariant(const size_t _maxIteration, const double _targetResidual, const double _minimalResidualDecrese, const bool _printProgress) 
                : maxInterations(_maxIteration), targetResidual(_targetResidual), minimalResidualDecrese(_minimalResidualDecrese), printProgress(_printProgress) { }
        
        /**
		 * @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		 * @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		 * @param _b the available measurments.
		 * @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		 */
        double operator()(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments, PerformanceData& _perfData = NoPerfData) const;
		
		/**
		 * @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		 * @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		 * @param _b the available measurments.
		 * @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		 */
        double operator()(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments, const std::vector<size_t>& _maxRanks, PerformanceData& _perfData = NoPerfData) const;
		
		/**
		 * @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		 * @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		 * @param _b the available measurments.
		 * @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		 */
        double operator()(TTTensor& _x, const RankOneMeasurmentSet& _measurments, PerformanceData& _perfData = NoPerfData) const;
    };
	
	/// @brief Default variant of the ADF algorithm
    const ADFVariant ADF(100000, 1e-8, 5e-4, true);
}

