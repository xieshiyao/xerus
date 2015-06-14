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

#include <xerus/algorithms/decompositionAls.h>
#include <xerus/basic.h>
#include <xerus/indexedTensorList.h>
#include <xerus/tensorNetwork.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_tensor_factorisations.h>
#include <xerus/misc/namedLogger.h>

namespace xerus {

	void decomposition_als(TTTensor& _x, const Tensor& _b, const double _eps, const size_t _maxIterations) {
		Index k, iU, iX, iL, rU, rL;
		FullTensor newXi;
		
		double lastResidual = frob_norm(FullTensor(_x) - _b);
        LOG(bla, "Initial residual " << lastResidual);
		
		for(size_t iteration = 0; iteration < _maxIterations; ++iteration) {
			// Move right
			for(size_t pos = 0; pos < _x.degree(); ++pos) {
                LOG(bla, "Optimizing position " << pos);
				_x.move_core(pos);
				std::pair<TensorNetwork, TensorNetwork> split = _x.chop(pos);
				_x.component(pos)(rU, iX, rL) = split.first(iU&1, rU)*split.second(rL, iL&1)*_b(iU^pos, iX, iL&(pos+1));
			}
			
			// Move left
			for(size_t pos = _x.degree()-2; pos > 0; --pos) {
                LOG(bla, "Optimizing position " << pos);
				_x.move_core(pos);
				std::pair<TensorNetwork, TensorNetwork> split = _x.chop(pos);
				_x.component(pos)(rU, iX, rL) = split.first(iU&1, rU)*split.second(rL, iL&1)*_b(iU^pos, iX, iL&(pos+1));
			}
			
			double residual = frob_norm(FullTensor(_x) - _b);
			LOG(bla, "New residual " << residual << ", residual change " << (lastResidual-residual)/residual);
			if(residual < EPSILON || (lastResidual-residual)/residual < _eps) { return; }
			lastResidual = residual;
		}
	}
}
