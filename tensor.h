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

#include "indexedTensorList.h"

namespace xerus {
	
    class Tensor {
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Vector containing the dimensions of the tensor
        std::vector<size_t> dimensions;
        
        /// Size of the Tensor -- always equal to the product of the dimensions
        size_t size;
        
        /// Single value representing a constant factor
        value_t factor;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Empty constructor which creates an order zero tensor
        implicit Tensor();
        
        /// Copy constructor
        implicit Tensor( const Tensor& );
        
        /// Move constructor
        implicit Tensor( Tensor&& _other );

        /// Creates a tensor with the given dimensions and factor
        ALLOW_MOVE(std::vector<size_t>, Vec)
        explicit Tensor(Vec&& _dimensions, const value_t _factor = 1.0) : dimensions(std::forward<Vec>(_dimensions)), size(product(dimensions)), factor(_factor)  { }
        
        /// Creates a tensor with the given dimensions and factor
        ALLOW_MOVE(std::initializer_list<size_t>, Init)
        explicit Tensor(Init&& _dimensions, const value_t _factor = 1.0) : dimensions(std::forward<Init>(_dimensions)), size(product(dimensions)), factor(_factor)  { }
        
        /// Returns a pointer containing a copy of the object with appropriate type
        virtual Tensor* get_copy() const = 0;
        
        /// Returns a pointer containing a moved copy of the object with appropriate type
        virtual Tensor* get_moved_copy() = 0;
        
        /// Returns a pointer to a newly constructed order zero tensor of appropriate type with entry equals zero
        virtual Tensor* construct_new() const = 0;
        
        /// Returns a pointer to a newly constructed tensor of appropriate type with all entries set to zero
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions) const = 0;
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions) const = 0;
        
        /// Returns a pointer to a newly constructed tensor of appropriate type with undefined entries
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) const = 0;
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) const = 0;
        
        /// Destructor
        virtual ~Tensor();
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Ensures that this tensor is the sole owner of its data space. If needed new space is allocated with the same entries is created
        virtual void ensure_own_data() = 0;
        
        /// Ensures that this tensor is the sole owner of its data space. If needed new space is allocated with undefined entries is created
        virtual void ensure_own_data_no_copy() = 0;
        
        /// Checks whether there is a non-trivial factor and applies it if nessecary.
        virtual void apply_factor() = 0;
        
        /// Checks whether there is a non-trivial factor and applies it. Even if no factor is applied ensure_own_data() is called.
        virtual void ensure_own_data_and_apply_factor() = 0;
        
        /// Resets the tensor with the given dimensions and undefined entries
        virtual void reset(const std::vector<size_t>&  _newDim, _unused_ DONT_SET_ZERO) = 0;
        virtual void reset(      std::vector<size_t>&& _newDim, _unused_ DONT_SET_ZERO) = 0;
        
        /// Resets the tensor with the given dimensions and all entries equals zero
        virtual void reset(const std::vector<size_t>&  _newDim) = 0;
        virtual void reset(      std::vector<size_t>&& _newDim) = 0;
    
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        virtual value_t& operator[](const size_t _position) = 0;
        virtual value_t operator[](const size_t _position) const = 0;
        
        virtual value_t& operator[](const std::vector<size_t>& _positions) = 0;
        virtual value_t operator[](const std::vector<size_t>& _positions) const = 0;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        template<typename... args>
        IndexedTensor<Tensor> operator()(args... _args) {
            return IndexedTensor<Tensor>(this, std::vector<Index>({_args...}), false);
        }
        
        template<typename... args>
        IndexedTensorReadOnly<Tensor> operator()(args... _args) const {
            return IndexedTensorReadOnly<Tensor>(this, std::vector<Index>({_args...}));
        }
        
        ALLOW_MOVE(std::vector<Index>, T)
        IndexedTensor<Tensor> operator()(T&&  _indices) {
            return IndexedTensor<Tensor>(this, std::forward<T>(_indices), false);
        }
        
        ALLOW_MOVE(std::vector<Index>, T)
        IndexedTensorReadOnly<Tensor> operator()(T&&  _indices) const {
            return IndexedTensorReadOnly<Tensor>(this, std::forward<T>(_indices));
        }
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Returns whether the object is of type SparseTensor
        virtual bool is_sparse() const = 0;
        
        /// Returns the degree of the tensor (==dimensions.size() )
        _inline_ size_t degree() const {
            return dimensions.size();
        }
        
        /// Returns whether the tensor has a non-trivial factor (i.e. factor != 1.0)
        _inline_ bool has_factor() const {
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wfloat-equal"
            return (factor != 1.0);
            #pragma GCC diagnostic pop
        }
        
        /// Returns the number of non-zero entries
        virtual size_t count_non_zero_entries(const value_t _eps = 1e-14) const = 0;
        
        /// Returns the frobenious norm of the tensor
        virtual value_t frob_norm() const = 0;
		
        /// Reinterprets the dimensions. Opposed to change_dimensions() it is assumed that the underlying data and the Size is NOT changed.
        ALLOW_MOVE(std::vector<size_t>, Vec)
        void reinterpret_dimensions( Vec&& _newDimensions) {
            REQUIRE(product(_newDimensions) == product(dimensions), "New dimensions must not change the size of the tensor in reinterpretation: " << product(_newDimensions) << " != " << product(dimensions));
            dimensions = std::forward<Vec>(_newDimensions);
        }
        
        
        /// Compares the Tensor entriewise to the given data
        virtual bool compare_to_data(std::vector<value_t> _values, const double _eps = 1e-14) const = 0;
        
        /// Compares the Tensor entriewise to the given data
        virtual bool compare_to_data(const value_t* _values, const double _eps = 1e-14) const = 0;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    protected:
        /// Assigns all member variables of Tensor
        void assign(const Tensor& _other);
        
        /// Move assigns all member variables of Tensor
        void assign(Tensor&& _other);
        
        /// Changes the dimensions of the tensor. Recalculates the size of the tensor.
        ALLOW_MOVE(std::vector<size_t>, Vec)
        void change_dimensions( Vec&& _newDimensions) {
            dimensions = std::forward<Vec>(_newDimensions);
            size = product(dimensions);
        }
    };
}
