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
* @brief Implementation of the indexed tensor / operator.
*/

#include <xerus/misc/internal.h>
#include <xerus/misc/containerSupport.h>

#include <xerus/indexedTensor_tensor_factorisations.h>
#include <xerus/index.h>
#include <xerus/indexedTensor.h>
#include <xerus/indexedTensorMoveable.h>
#include <xerus/tensor.h>

#include <xerus/blasLapackWrapper.h>
#include <xerus/cholmod_wrapper.h>

namespace xerus {

	void solve(internal::IndexedTensorWritable<Tensor>&& _x, internal::IndexedTensorReadOnly<Tensor>&& _a, internal::IndexedTensorReadOnly<Tensor>&& _b) {
		// x takes the dimensions of A -- also ensures that every index of x is contained in A
		_x.tensorObject->reset(_a.get_evaluated_dimensions(_x.indices), Tensor::Initialisation::None);
		
		_a.assign_indices();
		_b.assign_indices();
		
		size_t extraDims = 0;
		std::vector<Index> orderA;
		std::vector<Index> orderB;
		std::vector<Index> orderX;
		
		// If possible we don't want to reorder A, so first divide A into those shared with b and x. 
		for(const Index& idx : _a.indices) {
			if(misc::contains(_b.indices, idx)) {
				orderA.push_back(idx);
			} else {
				REQUIRE(misc::contains(_x.indices, idx), "Invalid indices");
				orderX.push_back(idx);
				for(size_t i = 0; i < idx.span; ++i) {
				}
			}
		}
		
		// This far the dimensions of A and B coincide
		orderB = orderA;
		
		// Add the second part of indices in the order obtained for X
		orderA.insert(orderA.end(), orderX.begin(), orderX.end());
		
		// Now complete indices of b and x with those not shared with A ( in order of b as we don't want to reorder b if possible).
		for(const Index& idx : _b.indices) {
			if(misc::contains(_x.indices, idx)) {
				orderB.push_back(idx);
				orderX.push_back(idx);
				for(size_t i = 0; i < idx.span; ++i) {
					extraDims++;
				}
			} else {
				REQUIRE(misc::contains(_a.indices, idx), "Invalid indices");
			}
		}
		
		Tensor reorderedA, reorderedB;
		
		reorderedA(orderA) = std::move(_a);
		reorderedB(orderB) = std::move(_b);
		
		solve_least_squares(*_x.tensorObject, reorderedA, reorderedB, extraDims);
		
		(*_x.tensorObject)(_x.indices) = (*_x.tensorObject)(orderX);
	}
	
	internal::IndexedTensorMoveable<Tensor> operator/ (internal::IndexedTensorReadOnly<Tensor>&& _b, internal::IndexedTensorReadOnly<Tensor>&& _A) {
		_A.assign_indices();
		_b.assign_indices();
		
		std::vector<Index> indicesX;
		for(const Index& idx : _A.indices) {
			if(!misc::contains(_b.indices, idx)) {
				indicesX.push_back(idx);
			}
		}
		
		internal::IndexedTensorMoveable<Tensor> tmpX(new Tensor(), std::move(indicesX));
		
		solve(std::move(tmpX), std::move(_A), std::move(_b));
		return tmpX;
	}
} // namespace xerus
