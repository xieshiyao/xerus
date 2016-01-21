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
* @brief Header file for the ALS algorithm and its variants.
*/

#pragma once

#include "../ttNetwork.h"
#include "../performanceData.h"

namespace xerus {

	/**
	* @brief Wrapper class for all ALS variants (dmrg etc.)
	* @details By creating a new object of this class and modifying the member variables, the behaviour of the solver can be modified.
	* The local solver will be ignored for all calls to the x=b variants (ie. without operator A)
	*/
	class ALSVariant {
	public:
		enum Direction { Increasing, Decreasing };
	protected:
		double solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, value_t _convergenceEpsilon, PerformanceData &_perfData = NoPerfData) const;
	
		struct ALSAlgorithmicData {
			const ALSVariant &ALS; ///< the algorithm this data belongs to
			const TTOperator *A; ///< global operator A
			TTTensor &x; ///< current iterate x
			const TTTensor &b; ///< global right-hand-side
			std::vector<size_t> targetRank; ///< rank for the final x
			std::vector<Tensor> xAxL, xAxR; ///< stacks for the local operator
			std::vector<Tensor> bxL, bxR; ///< stacks for the right-hand-side
			std::pair<size_t, size_t> optimizedRange; ///< range of indices for the nodes of _x that need to be optimized
			bool cannonicalizeAtTheEnd; ///< whether _x should be cannonicalized at the end
			size_t corePosAtTheEnd; ///< core position that should be restored at the end of the algorithm
			size_t currIndex; ///< position that is currently being optimized
			Direction direction; ///< direction of current sweep
			std::function<value_t()> energy_f; ///< the energy functional used for this calculation
			std::function<value_t()> residual_f; ///< the functional to calculate the current residual
			
			/**
			* @brief Finds the range of notes that need to be optimized and orthogonalizes @a _x properly
			* @details finds full-rank nodes (these can wlog be set to identity and need not be optimized)
			* requires cannonicalizeAtTheEnd and corePosAtTheEnd to be set
			* sets optimizedRange
			* modifies x
			*/
			void prepare_x_for_als();
			
			/**
			* @brief prepares the initial stacks for the local operator and local right-hand-side
			* @details requires optimziedRange
			* sets xAxL, xAxR, bxL, bxR
			*/
			void prepare_stacks();
			
			/**
			* @brief chooses the fitting energy functional according to settings and whether an operator A was given
			* @details
			* sets energy_f, residual_f
			*/
			void choose_energy_functional();
			
			ALSAlgorithmicData(const ALSVariant &, const TTOperator *, TTTensor &, const TTTensor &);
			
			/**
			* @brief performs one step in @a direction, updating all stacks
			* @details requires fully initialized data
			* modifies xAxL, xAxR, bxL, bxR, currIndex
			*/
			void move_to_next_index();
		};
		
		
	public:
		uint sites; ///< the number of sites that are simultaneously optimized
		size_t numHalfSweeps; ///< maximum number of sweeps to perform. set to 0 for infinite
		value_t convergenceEpsilon; ///< default value for the change in the energy functional at which the ALS assumes it is converged
// 		value_t minimumLocalResidual; ///< below this bound for the local residual, no local solver will be called.
		bool useResidualForEndCriterion; ///< calculates the residual to decide if the ALS converged. recommended if _perfdata is given
		bool preserveCorePosition; ///< if true the core will be moved to its original position at the end
		
		// TODO std::function endCriterion
		
		/// the algorithm that is used to solve the local problems
		using LocalSolver = std::function<void(const TensorNetwork &, std::vector<Tensor> &, const TensorNetwork &, const ALSAlgorithmicData &)>;
		LocalSolver localSolver;
		
		/// local solver that calls the corresponding lapack routines (LU solver)
		static void lapack_solver(const TensorNetwork &_A, std::vector<Tensor> &_x, const TensorNetwork &_b, const ALSAlgorithmicData &_data);
		
		//TODO add local CG solver
		
		/// fully defining constructor. alternatively ALSVariants can be created by copying a predefined variant and modifying it
		ALSVariant(uint _sites, size_t _numHalfSweeps, LocalSolver _localSolver,
			bool _useResidual=false
		) 
				: sites(_sites), numHalfSweeps(_numHalfSweeps), convergenceEpsilon(1e-6), 
				useResidualForEndCriterion(_useResidual), preserveCorePosition(true), localSolver(_localSolver)
		{
			REQUIRE(_sites>0, "");
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
		double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon, PerformanceData &_perfData = NoPerfData) const {
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
		
		double operator()(TTTensor &_x, const TTTensor &_b, PerformanceData &_perfData = NoPerfData) const {
			return solve(nullptr, _x, _b, numHalfSweeps, convergenceEpsilon, _perfData);
		}
	};
	
	/// default variant of the single-site ALS algorithm using the lapack solver
	extern const ALSVariant ALS;
	
	/// default variant of the two-site DMRG algorithm using the lapack solver
	extern const ALSVariant DMRG;
}

