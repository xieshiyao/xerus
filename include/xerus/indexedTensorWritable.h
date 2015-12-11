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
 * @brief Header file for the IndexedTensorWritable class.
 */

#pragma once

#include "indexedTensorReadOnly.h"

namespace xerus {
	// Necessary forward declaritons
	class Tensor;
	class TensorNetwork;
	
	/**
	 * @brief Abstract internal representation of an read and writeable indexed Tensor or TensorNetwork.
	 * @details This class (without specialization) should not appear or be used anywhere.
	 */
	template<class tensor_type>
	class IndexedTensorWritable : public IndexedTensorReadOnly<tensor_type> {
	public:
		/**
		 * @brief Non-const pointer to the tensor object. 
		 * @details This must always coincide with tensorObjectReadOnly.
		 * */
		tensor_type* tensorObject;
		
		///@brief Flag indicating, whether the IndexedTensorWritable is the owner of the tensorObject.
		bool deleteTensorObject;
		
	protected:
		///@brief There is no usefull default constructor;
		IndexedTensorWritable() = delete;
		
		///@brief There is no usefull copy constructor, because the handling of the tensorObject is unclear.
		IndexedTensorWritable(const IndexedTensorWritable &_other ) = delete;
		
		///@brief Move constructor.
		IndexedTensorWritable(IndexedTensorWritable &&_other );
		
		///@brief Constructs an IndexedTensorWritable with the given tensor and takes owership of the tensorObject if requested.
		IndexedTensorWritable(tensor_type* const _tensorObject, const std::vector<Index>&  _indices, const bool _takeOwnership);
		
		///@brief Constructs an IndexedTensorWritable with the given tensor and takes owership of the tensorObject if requested.
		IndexedTensorWritable(tensor_type* const _tensorObject,       std::vector<Index>&& _indices, const bool _takeOwnership);
		
	public:
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Destructor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		virtual ~IndexedTensorWritable();
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Other - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		///@brief Check whether the IndexedTensorWritable has ownership of the tensor object.
		bool is_owner() const;
		
		///@brief Deletes the tensor object if this is the owner.
		void delete_if_owner();
		
		///@brief Replacement for the original = operator
		void assign( IndexedTensorWritable&& _other);
		
		///@brief Resets the IndexedTensorWritable
		void reset(tensor_type* const _tensorObject, const std::vector<Index>&  _indices, const bool _takeOwnership = false);
		
		///@brief Resets the IndexedTensorWritable
		void reset(tensor_type* const _tensorObject,       std::vector<Index>&& _indices, const bool _takeOwnership = false);
		
		/**
		 * @brief Assignment operators -- Used for tensor assignment WITH indices.
		 * @details Note that this is NOT a classical assignment operator. In particular this and _rhs are NOT equivalent after
		 * the assignment.
		 */
		void operator=(IndexedTensorReadOnly<Tensor>&& _rhs);
		
		/**
		 * @brief Assignment operators -- Used for tensor assignment WITH indices.
		 * @details Note that this is NOT a classical assignment operator. In particular this and _rhs are NOT equivalent after
		 * the assignment.
		 */
		void operator=(IndexedTensorReadOnly<TensorNetwork>&& _rhs);
		
		///@brief The following would be deleted due to move constructor and is therefore implemented here, calls the IndexedTensorReadOnly version. 
		void operator=(IndexedTensorWritable<tensor_type>&& _rhs);
		
		/**
		 * @brief: Performes all traces induces by the current indices and therby also evaluates all fixed indices.
		 */
		void perform_traces();
	};
}
