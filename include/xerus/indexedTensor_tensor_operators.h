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

#pragma once

#include "indexedTensor_tensor.h"

namespace xerus {
    // Defined in the IndexedTensor classes: 
    // template<> void IndexedTensorWritable<Tensor>::operator=(const IndexedTensorReadOnly<Tensor>&  _rhs);
    
    IndexedTensorMoveable<Tensor> operator+ (const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs);
        
    IndexedTensorMoveable<Tensor> operator- (const IndexedTensorReadOnly<Tensor> & _lhs, const IndexedTensorReadOnly<Tensor> & _rhs);
    
    IndexedTensorMoveable<Tensor> operator* (const value_t _lhs, const IndexedTensorReadOnly<Tensor>& _rhs);

    IndexedTensorMoveable<Tensor> operator* (const IndexedTensorReadOnly<Tensor> & _lhs, const value_t _rhs);

    IndexedTensorMoveable<Tensor> operator/ (const IndexedTensorReadOnly<Tensor> & _lhs, const value_t _rhs);
    
    IndexedTensorMoveable<Tensor> operator/ (IndexedTensorReadOnly<Tensor> _b, IndexedTensorReadOnly<Tensor> _a);
}
