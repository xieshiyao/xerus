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
 * @brief Implementation of most indexed tensor operators.
 */

#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/tensor.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/misc/missingFunctions.h>

namespace xerus {

    template<>
    void IndexedTensorWritable<Tensor>::operator=(const IndexedTensorReadOnly<Tensor>&  _rhs)  {
        if(!_rhs.uses_tensor(tensorObject)) {
            // If LHS and RHS object don't coincide we can directly evaluate
            this->tensorObject->reset(_rhs.get_evaluated_dimensions(indices), DONT_SET_ZERO());
            evaluate(*this, _rhs);
        } else {
            // If the tensors in fact coincide we have to use a tmp object
            IndexedTensorMoveable<Tensor> tmpTensor(_rhs);
            this->tensorObject->reset(_rhs.get_evaluated_dimensions(indices), DONT_SET_ZERO());
            evaluate(*this, tmpTensor);
        }
    }


}
