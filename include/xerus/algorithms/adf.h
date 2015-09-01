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
		double solve(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments) const;
    
	public:
        size_t maxInterations; ///< maximum number of sweeps to perform. set to 0 for infinite
        value_t convergenceEpsilon; ///< default value for the change in the energy functional at which the ALS assumes it is converged
        bool printProgress; ///< informs the user about the current progress via std::cout (one continuously overwritten line)
        
		/// fully defining constructor. alternatively ALSVariants can be created by copying a predefined variant and modifying it
        ADFVariant(const size_t _maxIteration, const value_t _convergenceEpsilon, const bool _printProgress) 
                : maxInterations(_maxIteration), convergenceEpsilon(_convergenceEpsilon), printProgress(_printProgress) { }
        
        /**
		 * @brief Tries to reconstruct the (low rank) tensor _x from the given measurments. 
		 * @param[in,out] _x On input: an initial guess of the solution, also defining the ranks. On output: The reconstruction found by the algorithm.
		 * @param _b the available measurments.
		 * @returns the residual @f$|P_\Omega(x-b)|_2@f$ of the final @a _x.
		 */
        double operator()(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments) const {
			return solve(_x, _measurments);
		}
    };
	
	/// default variant of the single-site ALS algorithm using the lapack solver
    const ADFVariant ADF(100, EPSILON, true);
}

