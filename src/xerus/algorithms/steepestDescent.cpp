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
 * @brief Implementation of the steepest descent algorithms.
 */

#include <xerus/algorithms/steepestDescent.h>
#include <xerus/algorithms/als.h>
#include <xerus/basic.h>
#include <xerus/indexedTensorList.h>
#include <xerus/tensorNetwork.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_tensor_factorisations.h>

namespace xerus {
	void SteepestDescentVariant::HOSVDRetraction::operator()(TTTensor &_U, const TTTensor &_change) const {
		static const Index i;
		_U(i&0) = _U(i&0) + _change(i&0);
		if (roundByVector) {
			_U.round(rankVector);
		} else {
			_U.round(rank);
		}
	}
	
	void SteepestDescentVariant::ALSRetraction(TTTensor &_U, const TTTensor &_change) {
		static const ALSVariant roundingALS(1, 2, 0, ALSVariant::lapack_solver);
		static const Index i;
		TTTensor target;
		target(i&0) = _U(i&0) + _change(i&0);
		roundingALS(_U, target);
// 		_U.move_core(0);
	}
	
	//TODO this should use low-level calls and then be part of FullTensor
	static FullTensor pseudoInverse(const FullTensor &_A) {
		static const Index i1,i2,i3,i4;
		FullTensor U,S,V;
		(U(i1,i2), S(i2,i3), V(i3,i4)) = SVD(_A(i1,i4));
		S.modify_diag_elements([](value_t &a){a = 1/a;});
		FullTensor res;
		res(i1,i4) = V(i2,i1) * S(i2,i3) * U(i4,i3);
		return res;
	}
	
	void SteepestDescentVariant::SubmanifoldRetraction(TTTensor &_U, const TTTensor &_change) {
		REQUIRE(_U.cannonicalized && _U.corePosition == 0, "SubmanifoldRetraction is only implemented for core position 0 at the moment");
		REQUIRE(_U.dimensions == _change.dimensions, "");
		static const Index i1,i2,j1,j2,r,s,jx;
		std::vector<FullTensor> leftStackUV;
		std::vector<FullTensor> leftStackUU;
		FullTensor tmp({1,1}, [](){return 1.0;});
		leftStackUV.push_back(tmp);
		leftStackUU.push_back(tmp);
		for (size_t i=0; i<_U.degree()-1; ++i) {
			FullTensor newLeft;
			newLeft(j1,j2) = leftStackUV.back()(i1,i2) * _U.get_component(i)(i1,r,j1) * _change.get_component(i)(i2,r,j2);
			leftStackUV.emplace_back(std::move(newLeft));
			newLeft(j1,j2) = leftStackUU.back()(i1,i2) * _U.get_component(i)(i1,r,j1) * _U.get_component(i)(i2,r,j2);
			leftStackUU.emplace_back(std::move(newLeft));
		}
		TTTensor oldU(_U);
		FullTensor right(tmp);
		FullTensor UTV;
		for (size_t i=_U.degree(); i>0; --i) {
			const size_t currIdx = i-1;
			std::unique_ptr<FullTensor> newComponent(new FullTensor);
			const Tensor &UComp = oldU.component(currIdx);
			FullTensor V;
			V(i1,r,j1) = pseudoInverse(leftStackUU.back())(i1,s) * leftStackUV.back()(s,i2) * _change.get_component(currIdx)(i2,r,j2) * right(j1,j2);
			if (i!=_U.degree()) {
				V(i1,r,j1) = V(i1,r,j1) + UComp(i1,r,s)*UTV(s,j1);
			}
			if (currIdx!=0) {
				UTV(i1,i2) = V(i1,r,j1) * UComp(i2,r,j1);
				V(i1,r,j1) = V(i1,r,j1) - UTV(i1,s) * UComp(s,r,j1);
			}
			(*newComponent)(i1,r,j1) = (UComp(i1,r,j1) + V(i1,r,j1));
			_U.set_component(currIdx, std::move(newComponent));
			right(j1,j2) = oldU.get_component(currIdx)(j1,r,i1) * _change.get_component(currIdx)(j2,r,i2) * right(i1,i2);
			leftStackUV.pop_back();
			leftStackUU.pop_back();
		}
		REQUIRE(_U.is_valid_tt(), "ie");
// 		_U.assume_core_position(0);
		_U.move_core(0, true);
	}
	
	
	
	
	
	value_t SteepestDescentVariant::solve(const TTOperator *_Ap, TTTensor &_x, const TTTensor &_b, size_t _numSteps, value_t _convergenceEpsilon, PerformanceData &_perfData) const {
		const TTOperator &_A = *_Ap;
		static const Index i,j;
		size_t stepCount=0;
		TTTensor residual;
		value_t lastResidual=1e100;
		value_t currResidual=1e100;
		
		_perfData << "Steepest Descent for A*x = b, x.dimensions: " << _x.dimensions << '\n'
				  << "A.ranks: " << _A.ranks() << '\n'
				  << "x.ranks: " << _x.ranks() << '\n'
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
			if (printProgress && stepCount%37==0) {
				std::cout << "step \t" << stepCount << "\tresidual:" << currResidual << " (" 
				          << (lastResidual-currResidual) << ", " << std::abs(1-currResidual/lastResidual) 
						  << " vs " << _convergenceEpsilon << ")\r" << std::flush;
				std::cout << "                                                                               \r"; // note: not flushed so it will only erase content on next output
			}
		};
		updateResidualAndPerfdata();
		TTTensor y, Ay;
		value_t alpha;
		while ((_numSteps == 0 || stepCount < _numSteps)
			   && currResidual > _convergenceEpsilon
			   && (lastResidual-currResidual) > _convergenceEpsilon
			   && std::abs(1-currResidual/lastResidual) > _convergenceEpsilon) 
		{
			stepCount += 1;
			
			if (_Ap) {
				// search direction: y = A^T(b-Ax)
				y(i&0) = _A(j/2,i/2) * residual(j&0);
				// direction of change A*y
				Ay(i&0) = _A(i/2,j/2) * y(j&0);
				// "optimal" stepsize alpha = <y,y>/<Ay,Ay>
				alpha = misc::sqr(frob_norm(y)) / misc::sqr(frob_norm(Ay));
				// "optimal" change: alpha*y
				y *= alpha;
			} else {
				y = residual;
			}
			
			retraction(_x, y);
			
			updateResidualAndPerfdata();
		}
		
		return currResidual;
	}
	
}
