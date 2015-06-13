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

#include "indexedTensorWritable.h"

#include <type_traits>

namespace xerus {
	/**
	 * @brief Internal representation of an read and write and moveable indexed Tensor or TensorNetwork.
	 * @details This class appears inplicitly by indexing any Tensor or TensorNetwork. It is not recommended to use
	 * it explicitly or to store variables of this type (unless you really know what you are doing).
	 */
	template<class tensor_type>
    class IndexedTensorMoveable final : public IndexedTensorWritable<tensor_type> {
    public:
        /// Creates an empty indexed Tensor, should only be used internally
        IndexedTensorMoveable();
        
        /// There is no usefull copy constructor, because the handling of the tensorObject is unclear
        IndexedTensorMoveable(const IndexedTensorMoveable &_other ) = delete; // TODO actually not for 
        
        /// Move constructor
        IndexedTensorMoveable(IndexedTensorMoveable &&_other );
        
        /// Constructs an IndexedTensorMoveable with the given tensor and indices and if ordered to do so takes ownership of the tensorObject
        IndexedTensorMoveable(tensor_type* const _tensorObject, const std::vector<Index>& _indices);
        
        /// Constructs an IndexedTensorMoveable with the given tensor and indices and if ordered to do so takes ownership of the tensorObject
        IndexedTensorMoveable(tensor_type* const _tensorObject, std::vector<Index>&& _indices);
        
		/// Allow explicit cast to from IndexedTensorReadOnly.
		explicit IndexedTensorMoveable(const IndexedTensorReadOnly<tensor_type> &  _other);
		
		/// Allow explicit cast to from IndexedTensorReadOnly.
		explicit IndexedTensorMoveable(      IndexedTensorReadOnly<tensor_type>&&  _other);
		
        /// Allow conversions from indexed TensorNetworks to indexed Tensors
		template<class X = tensor_type, typename std::enable_if<std::is_base_of<Tensor, typename std::decay<X>::type>{}, int>::type = 0>
        explicit IndexedTensorMoveable(const IndexedTensorReadOnly<TensorNetwork> &  _other );
		
		template<class X = tensor_type, typename std::enable_if<std::is_base_of<Tensor, typename std::decay<X>::type>{}, int>::type = 0>
        explicit IndexedTensorMoveable(      IndexedTensorReadOnly<TensorNetwork> && _other );
        
        /// Allow implicit conversions from indexed Tensors to indexed TensorNetworks
		template<class X = tensor_type, typename std::enable_if<std::is_base_of<TensorNetwork, typename std::decay<X>::type>{}, int>::type = 0>
        IndexedTensorMoveable(const IndexedTensorReadOnly<Tensor> &  _other);
		
		template<class X = tensor_type, typename std::enable_if<std::is_base_of<TensorNetwork, typename std::decay<X>::type>{}, int>::type = 0>
        IndexedTensorMoveable(      IndexedTensorReadOnly<Tensor> && _other);
    };
}
