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

#include "basic.h"
#include <string>
#include <limits>
#include "indexedTensorReadOnly.h"
#include "indexedTensor.h"

namespace xerus {
	
    /// Abstract class which defines the common functionalities of the actual tensor classes FullTensor and SparseTensor.
    class Tensor {
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Vector containing the individual dimensions of the tensor.
        std::vector<size_t> dimensions;
        
        /// Size of the Tensor -- always equal to the product of the dimensions.
        size_t size;
        
        /// Single value representing a constant scaling factor.
        value_t factor;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Empty constructor which creates an order zero tensor.
        implicit Tensor();
        
        /// Copy constructor
        implicit Tensor( const Tensor& );
        
        /// Move constructor
        implicit Tensor( Tensor&& _other );

        /// Creates a tensor with the given dimensions and (optionally) the given scaling factor.
        explicit Tensor(const std::vector<size_t>& _dimensions, const value_t _factor = 1.0);
        
        /// Creates a tensor with the given dimensions and (optionally) the given scaling factor.
        explicit Tensor(std::vector<size_t>&& _dimensions, const value_t _factor = 1.0);
        
        /// Creates a tensor with the given dimensions and (optionally) the given scaling factor.
        explicit Tensor(std::initializer_list<size_t>&& _dimensions, const value_t _factor = 1.0);

        
        /// Returns a pointer containing a copy of the tensor with same type (i.e. FullTensor or SparseTensor).
        virtual Tensor* get_copy() const = 0;
        
        /// Returns a pointer containing a moved copy of the object with same type (i.e. FullTensor or SparseTensor).
        virtual Tensor* get_moved_copy() = 0;
        
        /// Returns a pointer to a newly constructed order zero tensor of same type (i.e. FullTensor or SparseTensor) with entry equals zero.
        virtual Tensor* construct_new() const = 0;
        
        /// Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with all entries set to zero.
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions) const = 0;
        
        /// Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with all entries set to zero.
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions) const = 0;
        
        /// Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with undefined entries.
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) const = 0;
        
        /// Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with undefined entries.
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) const = 0;
        
        /// Destructor
        virtual ~Tensor();
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Ensures that this tensor is the sole owner of its data. If needed new space is allocated and all entries are copied.
        virtual void ensure_own_data() = 0;
        
        /// Ensures that this tensor is the sole owner of its data space. If needed new space is allocated with entries left undefined.
        virtual void ensure_own_data_no_copy() = 0;
        
        /// Checks whether there is a non-trivial scaling factor and applies it if nessecary.
        virtual void apply_factor() = 0;
        
        /// Checks whether there is a non-trivial factor and applies it. Even if no factor is applied ensure_own_data() is called.
        virtual void ensure_own_data_and_apply_factor() = 0;
        
        /// Resets the tensor with the given dimensions and undefined entries.
        virtual void reset(const std::vector<size_t>&  _newDim, _unused_ DONT_SET_ZERO) = 0;
        
        /// Resets the tensor with the given dimensions and undefined entries.
        virtual void reset(      std::vector<size_t>&& _newDim, _unused_ DONT_SET_ZERO) = 0;
        
        /// Resets the tensor with the given dimensions and all entries equals zero.
        virtual void reset(const std::vector<size_t>&  _newDim) = 0;
        
        /// Resets the tensor with the given dimensions and all entries equals zero.
        virtual void reset(      std::vector<size_t>&& _newDim) = 0;
    
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Allows access to the entry at _position, assuming row-major ordering.
        virtual value_t& operator[](const size_t _position) = 0;
        
        /// Allows access to the entry at _position, assuming row-major ordering.
        virtual value_t operator[](const size_t _position) const = 0;
        
        /// Allows access to the entry at _position.
        virtual value_t& operator[](const std::vector<size_t>& _positions) = 0;
        
        /// Allows access to the entry at _position.
        virtual value_t operator[](const std::vector<size_t>& _positions) const = 0;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Indexes the tensor.
        template<typename... args>
        IndexedTensor<Tensor> operator()(args... _args) {
            return IndexedTensor<Tensor>(this, std::vector<Index>({_args...}), false);
        }
        
        /// Indexes the tensor.
        template<typename... args>
        IndexedTensorReadOnly<Tensor> operator()(args... _args) const {
            return IndexedTensorReadOnly<Tensor>(this, std::vector<Index>({_args...}));
        }
        
        /// Indexes the tensor.
        IndexedTensor<Tensor> operator()(const std::vector<Index>&  _indices);
        
        /// Indexes the tensor.
        IndexedTensor<Tensor> operator()(      std::vector<Index>&& _indices);
        
        /// Indexes the tensor.
        IndexedTensorReadOnly<Tensor> operator()(const std::vector<Index>&  _indices) const;
        
        /// Indexes the tensor.
        IndexedTensorReadOnly<Tensor> operator()(      std::vector<Index>&& _indices) const;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Checks whether the object is of type SparseTensor.
        virtual bool is_sparse() const = 0;
        
        /// Returns the degree of the tensor (==dimensions.size()).
        size_t degree() const;
        
        /// Returns whether the tensor has a non-trivial factor (i.e. factor != 1.0).
        bool has_factor() const;
        
        /// Returns the number of non-zero entries.
        virtual size_t count_non_zero_entries(const value_t _eps = std::numeric_limits<value_t>::epsilon()) const = 0;
        
        /// Returns an approximation of the reorder costs
        size_t reorder_costs() const;
        
        /// Calculates the frobenious norm of the tensor.
        virtual value_t frob_norm() const = 0;
		
        /// Reinterprets the dimensions of the tensor. Opposed to change_dimensions() it is assumed that the underlying data and the size are NOT changed.
        void reinterpret_dimensions(const std::vector<size_t>& _newDimensions);
        
        /// Reinterprets the dimensions of the tensor. Opposed to change_dimensions() it is assumed that the underlying data and the size are NOT changed.
        void reinterpret_dimensions(      std::vector<size_t>&& _newDimensions);
        
        /// Reinterprets the dimensions of the tensor. Opposed to change_dimensions() it is assumed that the underlying data and the size are NOT changed.
        void reinterpret_dimensions( std::initializer_list<size_t> _newDimensions);
        
        /// Returns a string representation of the tensor.
        virtual std::string to_string() const = 0;
        
        /// Compares the Tensor entriewise to the given data.
        virtual bool compare_to_data(std::vector<value_t> _values, const double _eps = std::numeric_limits<value_t>::epsilon()) const = 0;
        
        /// Compares the Tensor entriewise to the given data.
        virtual bool compare_to_data(const value_t* _values, const double _eps = std::numeric_limits<value_t>::epsilon()) const = 0;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    protected:
        /// Assigns all member variables of tensor.
        void assign(const Tensor& _other);
        
        /// Move assigns all member variables of tensor.
        void assign(Tensor&& _other);
        
        /// Changes the dimensions of the tensor and recalculates the size of the tensor.
        void change_dimensions(const std::vector<size_t>& _newDimensions);
        
        /// Changes the dimensions of the tensor and recalculates the size of the tensor.
        void change_dimensions(      std::vector<size_t>&& _newDimensions);
    };
    
    
    /// Returns the frobenius norm of the given tensor
    _inline_ value_t frob_norm(const Tensor& _tensor) { return _tensor.frob_norm(); }
}
