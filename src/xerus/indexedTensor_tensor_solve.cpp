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
 * @brief Implementation of the indexed tensor / operator.
 */

#include <xerus/indexedTensor_tensor_factorisations.h>
#include <xerus/index.h>
#include <xerus/tensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/blasLapackWrapper.h>

namespace xerus {

    void solve(const IndexedTensorWritable<Tensor>& _x, const IndexedTensorReadOnly<Tensor>& _a, const IndexedTensorReadOnly<Tensor>& _b) {
        // x takes the dimensions of A -- also ensures that every index of x is contained in A
		_x.tensorObject->reset(_a.get_evaluated_dimensions(_x.indices), Tensor::Initialisation::None);
        
        const std::vector<Index> AIndices = _a.get_assigned_indices();
        const std::vector<Index> bIndices = _b.get_assigned_indices();
		
        IF_CHECK( _x.check_indices(false); )
        
        REQUIRE(!_x.tensorObjectReadOnly->is_sparse() && !_b.tensorObjectReadOnly->is_sparse(), "At the moment we only allow Tensors in solve.");
        #ifdef _CHECK
            for(size_t i = 0; i < bAssIndices.numIndices; ++i) {
                REQUIRE(!bAssIndices.indexOpen[i] || contains(AAssIndices.indices, bAssIndices.indices[i]), "Every open index of b must be contained in A.");
                REQUIRE(!contains(_x.indices, bAssIndices.indices[i]), "x and b must not have indices in common.");
            }
        #endif
        
        // If possible we don't want to reorder A
        std::vector<Index> orderA;
        std::vector<Index> orderB;
        std::vector<Index> orderX;
        std::vector<size_t> dimensionsA;
        std::vector<size_t> dimensionsB;
        std::vector<size_t> dimensionsX;
        
        size_t dimensionsCount = 0;
        for(const Index& idx : AIndices) {
            if(misc::contains(bIndices, idx)) {
                orderA.push_back(idx);
                orderB.push_back(idx);
                for(size_t i = 0; i < idx.span; ++i) {
                    dimensionsA.push_back(_a.tensorObjectReadOnly->dimensions[dimensionsCount]);
                    dimensionsB.push_back(_a.tensorObjectReadOnly->dimensions[dimensionsCount]);
                    dimensionsCount++;
                }
            } else {
                orderX.push_back(idx);
                for(size_t i=0; i< (size_t) idx.span; ++i) {
                    dimensionsX.push_back(_a.tensorObjectReadOnly->dimensions[dimensionsCount++]);
                }
            }
        }
        orderA.insert(orderA.end(), orderX.begin(), orderX.end());
        dimensionsA.insert(dimensionsA.end(), dimensionsX.begin(), dimensionsX.end());
        
        //We need tmp objects for A and b, because Lapacke wants to destroys the input
		IndexedTensor<Tensor> tmpA(new Tensor(std::move(dimensionsA), Tensor::Representation::Dense, Tensor::Initialisation::None), orderA, true);
        evaluate(tmpA, _a);
        tmpA.tensorObject->ensure_own_data();
        
		IndexedTensor<Tensor> tmpB(new Tensor(std::move(dimensionsB), Tensor::Representation::Dense, Tensor::Initialisation::None), orderB, true);
        evaluate(tmpB, _b);
        tmpB.tensorObject->ensure_own_data();
        
        //Save slot for eventual tmpX
        std::unique_ptr<IndexedTensor<Tensor>> saveSlotX;
        const IndexedTensorWritable<Tensor>* usedX;
        if(orderX != _x.indices) {
			saveSlotX.reset(new IndexedTensor<Tensor>(new Tensor(std::move(dimensionsX), Tensor::Representation::Dense, Tensor::Initialisation::None), orderX, true));
            usedX = saveSlotX.get();
        } else {
            usedX = &_x;
        }
        
        // Assume A is an MxN matrix
        const size_t M = tmpB.tensorObjectReadOnly->size;
        const size_t N = usedX->tensorObjectReadOnly->size;
        
        if(tmpA.tensorObjectReadOnly->is_sparse()) {
            LOG(fatal, "Sparse solve not yet implemented.");
        } else {
            blasWrapper::solve_least_squares_destructive(static_cast<Tensor*>(usedX->tensorObject)->get_unsanitized_dense_data(), static_cast<Tensor*>(tmpA.tensorObject)->get_unsanitized_dense_data(), M, N, static_cast<Tensor*>(tmpB.tensorObject)->get_unsanitized_dense_data());
        }
        
        if(saveSlotX) { evaluate(_x, *usedX); }
        
        // Propagate the constant factor
        _x.tensorObject->factor = tmpB.tensorObjectReadOnly->factor / tmpA.tensorObjectReadOnly->factor;
    }

    IndexedTensorMoveable<Tensor> operator/ (IndexedTensorReadOnly<Tensor> _b, IndexedTensorReadOnly<Tensor> _A) {
        const std::vector<Index> indicesA = _A.get_assigned_indices();
        const std::vector<Index> indicesB = _b.get_assigned_indices();
        
        std::vector<Index> indicesX;
        std::vector<size_t> dimensionsX;
        
        size_t dimensionsCount = 0;
        for(const Index& idx : indicesA) {
            if(!misc::contains(indicesB, idx)) {
                indicesX.push_back(idx);
                for(size_t i = 0; i < idx.span; ++i) {
                    dimensionsX.push_back(_A.tensorObjectReadOnly->dimensions[dimensionsCount++]);
                }
            } else {
                dimensionsCount += idx.span;
            }
        }
        IndexedTensorMoveable<Tensor> tmpX(new Tensor(std::move(dimensionsX), Tensor::Representation::Dense, Tensor::Initialisation::None), std::move(indicesX));
        
        solve(tmpX, _A, _b);
        return tmpX;
    }
}
