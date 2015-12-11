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
 * @brief Header file for the IndexedTensorReadOnly class.
 */

#pragma once

#include "basic.h"
#include <vector>

namespace xerus {
	// Necessary forward declaritons
	class Index;
	template<class tensor_type> class IndexedTensorMoveable;

	
	/**
	 * @brief Internal representation of an readable indexed Tensor or TensorNetwork.
	 * @details This class appears inplicitly by indexing any Tensor or TensorNetwork. It is not recommended to use
	 * it explicitly or to store variables of this type (unless you really know what you are doing).
	 */
	template<class tensor_type>
	class IndexedTensorReadOnly {
	public:
		/// Pointer to the associated Tensor/TensorNetwork object
		const tensor_type* tensorObjectReadOnly;
		
		/// Vector of the associates indices 
		std::vector<Index> indices;
		
		/// Flag indicating whether the indices are assinged.
		bool indicesAssigned = false;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	public:
		/// There is no usefull default constructor
		IndexedTensorReadOnly() = delete;
		
		/// There is no usefull copy constructor for IndexedTensors.
		IndexedTensorReadOnly(const IndexedTensorReadOnly& _other ) = delete;
		
		/// Move-constructor
		IndexedTensorReadOnly(IndexedTensorReadOnly<tensor_type> && _other );
		
		/// Constructs an IndexedTensorReadOnly using the given pointer and indices.
		IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, const std::vector<Index>& _indices);
		
		/// Constructs an IndexedTensorReadOnly using the given pointer and indices.
		IndexedTensorReadOnly(const tensor_type* const _tensorObjectReadOnly, std::vector<Index>&& _indices);
		
		/// Destructor must be virtual
		virtual ~IndexedTensorReadOnly();
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	public:
		
		///@brief: IndexedTensorReadOnly cannot be assigned as they are read only.
		void operator=(const IndexedTensorReadOnly&  _rhs) = delete;
		
		
		///@brief: IndexedTensorReadOnly cannot be assigned as they are read only.
		void operator=(      IndexedTensorReadOnly&& _rhs) = delete;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Aritmetic Operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		// TODO use these and implement specialized sum for Tensor
// 		IndexedTensorMoveable<tensor_type> operator+(const IndexedTensorReadOnly& _other) const;
		
// 		IndexedTensorMoveable<tensor_type> operator+(IndexedTensorMoveable<tensor_type>&& _other) const;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Others - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		///@brief Allows cast to value_t if the degree of the current object is equal to 0.
		explicit operator value_t() const;
		
		///@brief Checks whether _otherTensor is the tensorObejct of this IndexTensor.
		bool uses_tensor(const tensor_type* _otherTensor) const;
		
		///@brief Returns the degree of the associated tensorObejct
		size_t degree() const;
		
		///@brief Assignes the indices using the degree of the tensorObejct.
		void assign_indices();
		
		///@brief Assignes the indices assuming the given degree.
		void assign_indices(const size_t _degree);
		
		///@brief Assignes the indices using the current tensorObejct.
		void assign_index_dimensions();
		
		bool is_contained_and_open(const Index& idx) const;
		
		std::vector<size_t> get_evaluated_dimensions(const std::vector<Index>& _indexOrder);
		
		#ifndef DISABLE_RUNTIME_CHECKS_
			/**
			* @brief: Checks whether the indices are usefull in combination with the current degree.
			*/
			void check_indices(const bool _allowNonOpen = true) const;
			
			/**
			* @brief: Checks whether the indices are usefull in combination with the given degree.
			*/
			void check_indices(const size_t _futureDegree, const bool _allowNonOpen) const;
		#endif
	};
	
	///@brief Returns the frobenious norm of the associated tensorObejct.
	template<class tensor_type>
	value_t frob_norm(const IndexedTensorReadOnly<tensor_type>& _idxTensor);
	
	size_t get_eval_degree(const std::vector<Index>& _indices);
	
}
