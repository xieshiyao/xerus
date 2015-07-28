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
		double solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, value_t _convergenceEpsilon, PerformanceData &_perfData = NoPerfData) const;
    
	public:
        size_t maxInterations; ///< maximum number of sweeps to perform. set to 0 for infinite
        value_t convergenceEpsilon; ///< default value for the change in the energy functional at which the ALS assumes it is converged
        bool printProgress; ///< informs the user about the current progress via std::cout (one continuously overwritten line)
        
		/// fully defining constructor. alternatively ALSVariants can be created by copying a predefined variant and modifying it
        ADFVariant(const size_t _maxIteration, const value_t _convergenceEpsilon, const bool _printProgress) 
                : maxInterations(_maxIteration), convergenceEpsilon(_convergenceEpsilon), printProgress(_printProgress) { }
        
        /**
		 * call to solve @f$ A\cdot x = b@f$ for @f$ x @f$ (in a least-squares sense)
		 * @param _A operator to solve for
		 * @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		 * @param _b right-hand side of the equation to be solved
		 * @param _convergenceEpsilon minimum change in residual / energy under which the algorithm terminates
		 * @param _perfData vector of performance data (residuals after every microiteration)
		 * @returns the residual @f$|Ax-b|@f$ of the final @a _x
		 */
        double operator()(TTTensor& _x, const std::vector<SinglePointMeasurment>& _measurments, const value_t _convergenceEpsilon) const {
			return solve(_x, _measurments);
		}
		
		/**
		 * call to solve @f$ A\cdot x = b@f$ for @f$ x @f$ (in a least-squares sense)
		 * @param _A operator to solve for
		 * @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		 * @param _b right-hand side of the equation to be solved
		 * @param _numHalfSweeps maximum number of half-sweeps to perform
		 * @param _perfData vector of performance data (residuals after every microiteration)
		 * @returns the residual @f$|Ax-b|@f$ of the final @a _x
		 */
        double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, PerformanceData &_perfData = NoPerfData) const {
			return solve(&_A, _x, _b, _numHalfSweeps, convergenceEpsilon, _perfData);
		}
		
		/**
		 * call to solve @f$ A\cdot x = b@f$ for @f$ x @f$ (in a least-squares sense)
		 * @param _A operator to solve for
		 * @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		 * @param _b right-hand side of the equation to be solved
		 * @param _perfData vector of performance data (residuals after every microiteration)
		 * @returns the residual @f$|Ax-b|@f$ of the final @a _x
		 */
        double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, PerformanceData &_perfData = NoPerfData) const {
			return solve(&_A, _x, _b, numHalfSweeps, convergenceEpsilon, _perfData);
		}
		
		/**
		 * call to minimze @f$ \|x - b\|^2 @f$ for @f$ x @f$
		 * @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		 * @param _b right-hand side of the equation to be solved
		 * @param _convergenceEpsilon minimum change in residual / energy under which the algorithm terminates
		 * @param _perfData vector of performance data (residuals after every microiteration)
		 * @returns the residual @f$|x-b|@f$ of the final @a _x
		 */
        double operator()(TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon, PerformanceData &_perfData = NoPerfData) const {
			return solve(nullptr, _x, _b, numHalfSweeps, _convergenceEpsilon, _perfData);
		}
		
		/**
		 * call to minimze @f$ \|x - b\|^2 @f$ for @f$ x @f$
		 * @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		 * @param _b right-hand side of the equation to be solved
		 * @param _numHalfSweeps maximum number of half-sweeps to perform
		 * @param _perfData vector of performance data (residuals after every microiteration)
		 * @returns the residual @f$|x-b|@f$ of the final @a _x
		 */
        double operator()(TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, PerformanceData &_perfData = NoPerfData) const {
			return solve(nullptr, _x, _b, _numHalfSweeps, convergenceEpsilon, _perfData);
		}
		
		double operator()(TTTensor &_x, const TTTensor &_b) const {
			return solve(nullptr, _x, _b, numHalfSweeps, convergenceEpsilon, NoPerfData);
		}
    };
	
	/// default variant of the single-site ALS algorithm using the lapack solver
    const ADFVariant ADF(1, 0, EPSILON);
}

