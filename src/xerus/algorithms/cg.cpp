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
	
	value_t CGVariant::solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numSteps, value_t _convergenceEpsilon, PerformanceData &_perfData) const {
		const TTOperator &_A = *_Ap;
		static const Index i,j;
		size_t stepCount=0;
		TTTensor residual;
		value_t lastResidual=1e100;
		value_t currResidual=1e100;
		
		if (_Ap) {
			_perfData << "Conjugated Gradients for ||A*x - b||^2, x.dimensions: " << _x.dimensions << '\n'
					<< "A.ranks: " << _A.ranks() << '\n';
		} else {
			_perfData << "Conjugated Gradients for ||x - b||^2, x.dimensions: " << _x.dimensions << '\n';
		}
		_perfData << "x.ranks: " << _x.ranks() << '\n'
					<< "b.ranks: " << _b.ranks() << '\n'
					<< "maximum number of steps: " << _numSteps << '\n'
					<< "convergence epsilon: " << _convergenceEpsilon << '\n';
		_perfData.start();
		
		auto updateResidualAndPerfdata = [&]() {
			lastResidual = currResidual;
			if (_Ap) {
				residual(i&0) = _b(i&0) - _A(i/2,j/2)*_x(j&0);
			} else {
				residual = _b - _x;
			}
			currResidual = frob_norm(residual);
			_perfData.add(currResidual);
			if (printProgress && stepCount%1==0) {
				std::cout << "step \t" << stepCount << "\tresidual:" << currResidual << " (" 
						  << (lastResidual-currResidual) << ", " << std::abs(1-currResidual/lastResidual) 
						  << " vs " << _convergenceEpsilon << ")\r" << std::flush;
				std::cout << "                                                                               \r"; // note: not flushed so it will only erase content on next output
			}
		};
		updateResidualAndPerfdata();
		while ((_numSteps == 0 || stepCount < _numSteps)
			   && currResidual > _convergenceEpsilon
			   && std::abs(lastResidual-currResidual) > _convergenceEpsilon
			   && std::abs(1-currResidual/lastResidual) > _convergenceEpsilon) 
		{
			TTTensor direction(residual);
			
			while ((_numSteps == 0 || stepCount < _numSteps)
				&& (restartInterval == 0 || (stepCount+1)%restartInterval != 0)
				&& currResidual > _convergenceEpsilon
				&& std::abs(lastResidual-currResidual) > _convergenceEpsilon
				&& std::abs(1-currResidual/lastResidual) > _convergenceEpsilon) 
			{
				stepCount += 1;
				TTTensor change;
				value_t alpha;
				
	// 			change(i&0) = _A(i/2,j/2) * direction(j&0);
				alpha = currResidual/value_t(direction(i&0)*_A(i/2,j/2)*direction(j&0));//direction(i&0)*change(i&0));
	// 			LOG(alpha, alpha << " " << currResidual << " / " << value_t(direction(i&0)*_A(i/2,j/2)*direction(j&0)) << " " << value_t(direction(i&0)*change(i&0)));
				retraction(_x, alpha * direction);
	// 			retraction(residual, -alpha * change);
				updateResidualAndPerfdata();
				
				double beta = currResidual / lastResidual;
	// 			direction(i&0) = residual(i&0) + beta * direction(i&0);
	// 			TTTensor oldDirection(direction);
	// 			direction = residual;
				direction *= beta;
				retraction(direction, residual);
	// 			retraction(direction, beta * oldDirection);
			}
		}
		
		return currResidual;
	}
	/*
	 * 		change(j^d) = _A(j^d,k^d) * dir(k^d);
		double alpha = fnorm/double(dir(i^d)*change(i^d));
		_x(i^d) = _x(i^d) + alpha * dir(i^d);
		res(i^d) = res(i^d) - alpha*change(i^d);
		double newfnorm = frob_norm(res(i^d));
		double beta = newfnorm / fnorm;
		dir(i^d) = res(i^d) + beta * dir(i^d);
		fnorm = newfnorm;
// 		dir.round(1e-14);
// 		_x.round(1e-14);
// 		res.round(1e-14);
		if (_perfData) _perfData->push_back(fnorm);
		std::ofstream out("test.dat", std::ios_base::app | std::ios_base::ate);
		out << fnorm << std::endl;
		out.close();
		LOG(cg, fnorm);*/
}
