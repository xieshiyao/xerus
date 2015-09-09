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
 * @brief Implementation of the CG algorithm.
 */

#include <xerus/algorithms/cg.h>
#include <xerus/algorithms/als.h>
#include <xerus/basic.h>
#include <xerus/indexedTensorList.h>
#include <xerus/tensorNetwork.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_tensor_factorisations.h>

namespace xerus {
	
	value_t GeometricCGVariant::solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numSteps, value_t _convergenceEpsilon, PerformanceData &_perfData) const {
		const TTOperator &_A = *_Ap;
		static const Index i,j;
		size_t stepCount=0;
		TTTensor residual;
		value_t lastResidual=1e100;
		value_t currResidual=1e100;
		value_t normB = frob_norm(_b);
		
		if (_Ap) {
			_perfData << "Conjugated Gradients for ||A*x - b||^2, x.dimensions: " << _x.dimensions << '\n'
					<< "A.ranks: " << _A.ranks() << '\n';
			if (assumeSymmetricPositiveDefiniteOperator) {
				_perfData << " with symmetric positive definite Operator A\n";
			}
		} else {
			_perfData << "Conjugated Gradients for ||x - b||^2, x.dimensions: " << _x.dimensions << '\n';
		}
		_perfData << "x.ranks: " << _x.ranks() << '\n'
					<< "b.ranks: " << _b.ranks() << '\n'
					<< "maximum number of steps: " << _numSteps << '\n'
					<< "convergence epsilon: " << _convergenceEpsilon << '\n';
		_perfData.start();
		
		auto updateResidual = [&]() {
			if (_Ap) {
				residual(i&0) = _b(i&0) - _A(i/2,j/2)*_x(j&0);
			} else {
				residual = _b - _x;
			}
			currResidual = frob_norm(residual)/normB;
		};
		
		auto updatePerfdata = [&]() {
			_perfData.add(currResidual);
			if (printProgress && stepCount%1==0) {
				std::cout << "step \t" << stepCount << "\tresidual:" << currResidual << " (" 
						  << (lastResidual-currResidual) << ", " << std::abs(1-currResidual/lastResidual) 
						  << " vs " << _convergenceEpsilon << ")\n" << std::flush;
				std::cout << "                                                                               \r"; // note: not flushed so it will only erase content on next output
			}
		};
		updateResidual();
		updatePerfdata();
		
		TTTangentVector direction(_x, residual);
		while ((_numSteps == 0 || stepCount < _numSteps)
			&& currResidual > _convergenceEpsilon
			&& std::abs(lastResidual-currResidual) > _convergenceEpsilon
			&& std::abs(1-currResidual/lastResidual) > _convergenceEpsilon) 
		{
			if (restartInterval != 0 && (stepCount+1)%restartInterval == 0) {
				direction = TTTangentVector(_x, residual);
			}
			stepCount += 1;
			value_t alpha;
			
			if (_Ap) {
				TTTensor dirTT = TTTensor(direction);
				if (assumeSymmetricPositiveDefiniteOperator) {
					alpha = direction.scalar_product(direction)/value_t(dirTT(i&0)*_A(i/2,j/2)*dirTT(j&0));//direction(i&0)*change(i&0));
				} else {
					value_t normDir = dirTT.frob_norm();
					dirTT(i&0) = _A(i/2,j/2) * dirTT(j&0);
					alpha = normDir/frob_norm(dirTT);
				}
			} else {
				alpha = 1;
			}
			
			TTTensor oldX(_x);
			retraction(_x, direction * alpha);
			lastResidual = currResidual;
			updateResidual();
			
			// armijo backtracking
			while (alpha > 1e-20 && lastResidual < currResidual - 1e-4 * lastResidual * alpha * value_t(residual(i&0) * TTTensor(direction)(i&0))) {
// 				LOG(huch, alpha << " " << currResidual);
				alpha /= 2;
				_x = oldX;
				retraction(_x, direction * alpha);
				updateResidual();
			}
			
			updatePerfdata();
			
// 			direction(i&0) = residual(i&0) + beta * direction(i&0);
			TTTangentVector oldDirection(direction);
			value_t oldDirNorm = oldDirection.scalar_product(oldDirection);
			vectorTransport(_x, direction);
			if (assumeSymmetricPositiveDefiniteOperator || !_Ap) {
				direction = TTTangentVector(_x, residual);
			} else {
				TTTensor grad;
				grad(i&0) = _A(j/2,i/2) * residual(j&0); // grad = A^T * (b - Ax)
				direction = TTTangentVector(_x, grad);
			}
			double beta = direction.scalar_product(direction) / oldDirNorm ;//currResidual / lastResidual; // Fletcher-Reeves update
// 			LOG(ab, "\t\t\t" << alpha << " " << beta);
			direction += oldDirection * (beta);
		}
		
		return currResidual;
	}
}
