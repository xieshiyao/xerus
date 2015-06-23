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
 * @brief Header file for the steepest descent algorithms.
 */

#pragma once

#include "../ttNetwork.h"
#include "../performanceData.h"

namespace xerus {

	/**
	* @brief Wrapper class for all steepest descent variants
	* @note only implemented for TTTensors at the moment.
	* @details By creating a new object of this class and modifying the member variables, the behaviour of the solver can be modified.
	*/
	class SteepestDescentVariant {
	protected:
		double solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numSteps, value_t _convergenceEpsilon, PerformanceData &_perfData = NoPerfData) const;
	
	public:
		size_t numSteps; ///< maximum number of steps to perform. set to 0 for infinite
		value_t convergenceEpsilon; ///< default value for the change in the residual at which the algorithm assumes it is converged
		bool printProgress; ///< informs the user about the current progress via std::cout (one continuously overwritten line)
		bool assumeSymmetricPositiveDefiniteOperator; ///< calculates the gradient as b-Ax instead of A^T(b-Ax)
		
		std::function<void(TTTensor &, const TTTensor &)> retraction; ///< the retraction to project from point + tangent vector to a new point on the manifold
		
		// TODO preconditioner
		
		/// retraction that performs a HOSVD to project back onto the Manifold
		struct HOSVDRetraction {
			bool roundByVector;
			size_t rank;
			std::vector<size_t> rankVector;
			void operator()(TTTensor &_U, const TTTensor &_change) const;
			HOSVDRetraction(size_t _rank) : roundByVector(false), rank(_rank) {}
			HOSVDRetraction(const std::vector<size_t> &_rank) : roundByVector(true), rank(~0ul), rankVector(_rank) {}
		};
		
		/// retraction that performs an ALS half-sweep to project back onto the Manifold. Automatically retains the ranks of @a _U
		static void ALSRetraction(TTTensor &_U, const TTTensor &_change);
		
		/// retraction that performs componentwise addition of @f$ U_i @f$ and @f$ W_i @f$ where @f$ W_i @f$ is the i-th component of the riemannian tangential vector representation
		static void SubmanifoldRetraction(TTTensor &_U, const TTTensor &_change);
		
		/// fully defining constructor. alternatively SteepestDescentVariant can be created by copying a predefined variant and modifying it
		SteepestDescentVariant(size_t _numSteps, value_t _convergenceEpsilon, bool _symPosOp, std::function<void(TTTensor &, const TTTensor &)> _retraction)
				: numSteps(_numSteps), convergenceEpsilon(_convergenceEpsilon),
				  assumeSymmetricPositiveDefiniteOperator(_symPosOp), retraction(_retraction)
		{ }
		
		/// definition using only the retraction. In the following an operator() including either convergenceEpsilon or numSteps must be called or the algorithm will never terminate
		SteepestDescentVariant(std::function<void(TTTensor &, const TTTensor &)> _retraction)
				: numSteps(0), convergenceEpsilon(0.0), retraction(_retraction)
		{ }
		
		/**
		* call to solve @f$ A\cdot x = b@f$ for @f$ x @f$ (in a least-squares sense)
		* @param _A operator to solve for
		* @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		* @param _b right-hand side of the equation to be solved
		* @param _convergenceEpsilon minimum change in residual / energy under which the algorithm terminates
		* @param _perfData vector of performance data (residuals after every microiteration)
		* @returns the residual @f$|Ax-b|@f$ of the final @a _x
		*/
		double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon, PerformanceData &_perfData = NoPerfData) const {
			return solve(&_A, _x, _b, numSteps, _convergenceEpsilon, _perfData);
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
		double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numSteps, PerformanceData &_perfData = NoPerfData) const {
			return solve(&_A, _x, _b, _numSteps, convergenceEpsilon, _perfData);
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
			return solve(&_A, _x, _b, numSteps, convergenceEpsilon, _perfData);
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
			return solve(nullptr, _x, _b, numSteps, _convergenceEpsilon, _perfData);
		}
		
		/**
		* call to minimze @f$ \|x - b\|^2 @f$ for @f$ x @f$
		* @param[in,out] _x in: initial guess, out: solution as found by the algorithm
		* @param _b right-hand side of the equation to be solved
		* @param _numHalfSweeps maximum number of half-sweeps to perform
		* @param _perfData vector of performance data (residuals after every microiteration)
		* @returns the residual @f$|x-b|@f$ of the final @a _x
		*/
		double operator()(TTTensor &_x, const TTTensor &_b, size_t _numSteps, PerformanceData &_perfData = NoPerfData) const {
			return solve(nullptr, _x, _b, _numSteps, convergenceEpsilon, _perfData);
		}
		
		double operator()(TTTensor &_x, const TTTensor &_b) const {
			return solve(nullptr, _x, _b, numSteps, convergenceEpsilon, NoPerfData);
		}
	};
	
	/// default variant of the steepest descent algorithm using the lapack solver
	const SteepestDescentVariant SteepestDescent(0, 1e-8, false, SteepestDescentVariant::SubmanifoldRetraction);
}

