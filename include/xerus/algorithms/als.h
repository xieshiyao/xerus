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
 * Header file for the ALS algorithm and its variants.
 */

#pragma once

#include "../ttNetwork.h"

namespace xerus {

	/**
	 * @brief wrapper class for all ALS variants (dmrg etc.)
	 * 
	 * By creating a new object of this class and modifying the member variables, the behaviour of the solver can be modified.
	 * The local solver will be ignored for all calls to the x=b variants (ie. without operator A)
	 */
    class ALSVariant {
	protected:
		double solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, value_t _convergenceEpsilon,  std::vector<value_t> *_perfData) const;
    public:
        uint sites; ///< the number of sites that are simultaneously optimized
        size_t numHalfSweeps; ///< maximum number of sweeps to perform. set to 0 for infinite
        value_t convergenceEpsilon; ///< default value for the change in the energy functional at which the ALS assumes it is converged
        value_t minimumLocalResidual; ///< below this bound for the local residual, no local solver will be called.
        bool printProgress; ///< informs the user about the current progress via std::cout (one continuously overwritten line)
        bool useResidualForEndCriterion; ///< calculates the residual to decide if the ALS converged. recommended if _perfdata is given
        bool preserveCorePosition; ///< if true the core will be moved to its original position at the end
        
        // TODO printProgressToFile;
        // TODO std::function endCriterion
        // TODO all arguments should be tensorNetworks to allow efficient CG (esp. in the dmrg case)
        
        /// the algorithm that is used to solve the local problems
        std::function<void(const TensorNetwork &, Tensor &, const Tensor &)> localSolver;
        
		/// local solver that calls the corresponding lapack routines (LU solver)
        static void lapack_solver(const TensorNetwork &_A, Tensor &_x, const Tensor &_b);
		
		//TODO add local CG solver
        
		/// fully defining constructor. alternatively ALSVariants can be created by copying a predefined variant and modifying it
        ALSVariant(uint _sites, value_t _minimumLocalResidual, std::function<void(const TensorNetwork &, Tensor &, const Tensor &)> _localSolver) 
                : sites(_sites), numHalfSweeps(0), convergenceEpsilon(1e-6), 
                minimumLocalResidual(_minimumLocalResidual), printProgress(false), 
                useResidualForEndCriterion(false), preserveCorePosition(true), localSolver(_localSolver)
		{
            REQUIRE(_sites>0, "");
            REQUIRE(_minimumLocalResidual>=0, "");
        }
        
        /**
		 * call to solve @f$ A\cdot x = b@f$ for @f$ x @f$ (in a least-squares sense)
		 * @param _A operator to solve for
		 * @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		 * @param _b right-hand side of the equation to be solved
		 * @param _convergenceEpsilon minimum change in residual / energy under which the algorithm terminates
		 * @param _perfData vector of performance data (residuals after every microiteration)
		 * @returns the residual @f$|Ax-b|@f$ of the final @a _x
		 */
        double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon,  std::vector<value_t> *_perfData = nullptr) const {
			return solve(&_A, _x, _b, numHalfSweeps, _convergenceEpsilon, _perfData);
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
        double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps,  std::vector<value_t> *_perfData = nullptr) const {
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
        double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, std::vector<value_t> *_perfData = nullptr) const {
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
        double operator()(TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon,  std::vector<value_t> *_perfData = nullptr) const {
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
        double operator()(TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps,  std::vector<value_t> *_perfData = nullptr) const {
			return solve(nullptr, _x, _b, _numHalfSweeps, convergenceEpsilon, _perfData);
		}
		
		double operator()(TTTensor &_x, const TTTensor &_b) const {
			return solve(nullptr, _x, _b, numHalfSweeps, convergenceEpsilon, nullptr);
		}
    };
	
	/// default variant of the single-site ALS algorithm using the lapack solver
    const ALSVariant ALS(1, 1e-13, ALSVariant::lapack_solver);
	
	/// default variant of the two-site DMRG algorithm using the lapack solver
    const ALSVariant DMRG(2, 1e-13, ALSVariant::lapack_solver);
}

