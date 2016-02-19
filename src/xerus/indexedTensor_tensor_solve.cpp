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
		
//         IF_CHECK( _x.check_indices(false); )
		
// 		for(size_t i = 0; i < bAssIndices.numIndices; ++i) {
// 			REQUIRE(!bAssIndices.indexOpen[i] || contains(AAssIndices.indices, bAssIndices.indices[i]), "Every open index of b must be contained in A.");
// 			REQUIRE(!contains(_x.indices, bAssIndices.indices[i]), "x and b must not have indices in common.");
// 		}
		
		// If possible we don't want to reorder A
		std::vector<Index> orderA;
		std::vector<Index> orderB;
		std::vector<Index> orderX;
		std::vector<size_t> dimensionsA;
		std::vector<size_t> dimensionsB;
		std::vector<size_t> dimensionsX;
		
		size_t dimensionsCount = 0;
		for(const Index& idx : _a.indices) {
			if(misc::contains(_b.indices, idx)) {
				orderA.push_back(idx);
				orderB.push_back(idx);
				for(size_t i = 0; i < idx.span; ++i) {
					dimensionsA.push_back(_a.tensorObjectReadOnly->dimensions[dimensionsCount]);
					dimensionsB.push_back(_a.tensorObjectReadOnly->dimensions[dimensionsCount]);
					dimensionsCount++;
				}
			} else {
				orderX.push_back(idx);
				for(size_t i = 0; i < idx.span; ++i) {
					dimensionsX.push_back(_a.tensorObjectReadOnly->dimensions[dimensionsCount++]);
				}
			}
		}
		orderA.insert(orderA.end(), orderX.begin(), orderX.end());
		dimensionsA.insert(dimensionsA.end(), dimensionsX.begin(), dimensionsX.end());
		
		//Save slot for eventual tmpX
		std::unique_ptr<internal::IndexedTensor<Tensor>> saveSlotX;
		internal::IndexedTensorWritable<Tensor>* usedX;
		
		bool sparseResult = (_a.tensorObjectReadOnly->is_sparse() && _b.tensorObjectReadOnly->is_sparse());
		
		if (orderX != _x.indices || (_x.tensorObject->is_sparse() != sparseResult)) {
			if (sparseResult) {
				saveSlotX.reset(new internal::IndexedTensor<Tensor>(
					new Tensor(std::move(dimensionsX), Tensor::Representation::Sparse, Tensor::Initialisation::None), orderX, true));
			} else {
				saveSlotX.reset(new internal::IndexedTensor<Tensor>(
					new Tensor(std::move(dimensionsX), Tensor::Representation::Dense, Tensor::Initialisation::None), orderX, true));
			}
			
			usedX = saveSlotX.get();
		} else {
			usedX = &_x;
		}
		
		// Assume A is an MxN matrix
		const size_t M = _b.tensorObjectReadOnly->size;
		const size_t N = usedX->tensorObjectReadOnly->size;
		
		if (_a.tensorObjectReadOnly->is_sparse()) {
			if (_b.tensorObjectReadOnly->is_sparse()) {
				internal::CholmodSparse::solve_sparse_rhs(
					usedX->tensorObject->get_unsanitized_sparse_data(), N, 
					_a.tensorObjectReadOnly->get_unsanitized_sparse_data(), false,
					_b.tensorObjectReadOnly->get_unsanitized_sparse_data(), M);
			} else {
				internal::CholmodSparse::solve_dense_rhs(
					usedX->tensorObject->get_unsanitized_dense_data(), N, 
					_a.tensorObjectReadOnly->get_unsanitized_sparse_data(), false,
					_b.tensorObjectReadOnly->get_unsanitized_dense_data(), M);
			}
			
			// Propagate the constant factor
			usedX->tensorObject->factor = _b.tensorObjectReadOnly->factor / _a.tensorObjectReadOnly->factor;
		} else {
			// dense A
			// We need tmp objects for A and b, because Lapacke wants to destroys the input
			internal::IndexedTensor<Tensor> tmpA(new Tensor(std::move(dimensionsA), Tensor::Representation::Dense, Tensor::Initialisation::None), orderA, true);
			evaluate(std::move(tmpA), std::move(_a));
			tmpA.tensorObject->ensure_own_data();
			
			internal::IndexedTensor<Tensor> tmpB(new Tensor(std::move(dimensionsB), Tensor::Representation::Dense, Tensor::Initialisation::None), orderB, true);
			evaluate(std::move(tmpB), std::move(_b));
			tmpB.tensorObject->use_dense_representation();
			tmpB.tensorObject->ensure_own_data();
			
			blasWrapper::solve_least_squares_destructive(
				usedX->tensorObject->override_dense_data(), 
				tmpA.tensorObject->get_unsanitized_dense_data(), M, N, 
				tmpB.tensorObject->get_unsanitized_dense_data());
			
			// Propagate the constant factor
			usedX->tensorObject->factor = tmpB.tensorObjectReadOnly->factor / tmpA.tensorObjectReadOnly->factor;
		}
		
		if(saveSlotX) { evaluate(std::move(_x), std::move(*usedX)); }
	}

	internal::IndexedTensorMoveable<Tensor> operator/ (internal::IndexedTensorReadOnly<Tensor>&& _b, internal::IndexedTensorReadOnly<Tensor>&& _A) {
		_A.assign_indices();
		_b.assign_indices();
		
		std::vector<Index> indicesX;
		std::vector<size_t> dimensionsX;
		
		size_t dimensionsCount = 0;
		for(const Index& idx : _A.indices) {
			if(!misc::contains(_b.indices, idx)) {
				indicesX.push_back(idx);
				for(size_t i = 0; i < idx.span; ++i) {
					dimensionsX.push_back(_A.tensorObjectReadOnly->dimensions[dimensionsCount++]);
				}
			} else {
				dimensionsCount += idx.span;
			}
		}
		internal::IndexedTensorMoveable<Tensor> tmpX(new Tensor(std::move(dimensionsX), Tensor::Representation::Dense, Tensor::Initialisation::None), std::move(indicesX));
		
		solve(std::move(tmpX), std::move(_A), std::move(_b));
		return tmpX;
	}
}
