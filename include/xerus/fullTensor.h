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

#include <memory>

#include "tensor.h"
#include "misc/selectedFunctions.h"
#include "misc/sfinae.h"

namespace xerus {
    class SparseTensor;
    
    /// The main class used to represent any dense tensor.
    class FullTensor final : public Tensor {
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Shared pointer to the data array with size "size". 
        /// The data is stored such that indices increase from right to left (row-major order).
        /// If the tensor is modified and not sole owner a deep copy is performed.
        std::shared_ptr<value_t> data;
         
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Empty constructor, creates an order zero tensor with zero as single entry.
        explicit FullTensor();
        
        /// Copy constructor.
        implicit FullTensor(const FullTensor & _other);
        
        /// Move constructor.
        implicit FullTensor(      FullTensor&& _other);
        
        /// Creates a tensor of given degree, with all dimensions euqals one and zero as the single entry.
        explicit FullTensor(const size_t _degree);
        
        /// Creates a tensor with the given dimensions and undefined entries.
        explicit FullTensor(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO);
        
        /// Creates a tensor with the given dimensions and undefined entries.
        explicit FullTensor(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO);
        
        /// Creates a tensor with the given dimensions and all entries equals zero.
        explicit FullTensor(const std::vector<size_t>&  _dimensions);
        
        /// Creates a tensor with the given dimensions and all entries equals zero.
        explicit FullTensor(      std::vector<size_t>&& _dimensions);
        
        /// Creates a tensor with the given dimensions and uses the given data as entries.
        template<ADD_MOVE(std::vector<size_t>, Vec), ADD_MOVE(std::shared_ptr<value_t>, SPtr)>
        explicit FullTensor(Vec&& _dimensions, SPtr&& _data) : Tensor(std::forward<Vec>(_dimensions)), data(std::forward<SPtr>(_data)) { }
        
        /// Creates a tensor with the given dimensions and uses the given data as entries.
        explicit FullTensor(const std::vector<size_t> & _dimensions, std::unique_ptr<value_t[]>&& _data);
        
        /// Creates a tensor with the given dimensions and uses the given data as entries.
        explicit FullTensor(      std::vector<size_t>&& _dimensions, std::unique_ptr<value_t[]>&& _data);
        
        
        /// Creates a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
        ALLOW_MOVE(std::vector<size_t>, Vec)
        explicit FullTensor(Vec&& _dimensions, const std::function<value_t()>& _f) : FullTensor(std::forward<Vec>(_dimensions), DONT_SET_ZERO()) {
            value_t* realData = data.get();
            for (size_t i=0; i < size; ++i) {
                realData[i] = _f();
            }
        }
        
        /// Creates a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
        ALLOW_MOVE(std::vector<size_t>, Vec)
        explicit FullTensor(Vec&& _dimensions, const std::function<value_t(const size_t)>& _f) : FullTensor(std::forward<Vec>(_dimensions), DONT_SET_ZERO()) {
            value_t* realData = data.get();
            for (size_t i=0; i < size; ++i) {
                realData[i] = _f(i);
            }
        }
        
        /// Creates a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
        ALLOW_MOVE(std::vector<size_t>, Vec)
        explicit FullTensor(Vec&& _dimensions, const std::function<value_t(const std::vector<size_t>&)>& _f) : FullTensor(std::forward<Vec>(_dimensions), DONT_SET_ZERO()) {
            value_t* realData = data.get();
            std::vector<size_t> multIdx(degree(), 0);
            size_t idx = 0;
            while (true) {
                realData[idx] = _f(multIdx);
                // increasing indices
                idx++;
                size_t changingIndex = degree()-1;
                multIdx[changingIndex]++;
                while(multIdx[changingIndex] == dimensions[changingIndex]) {
                    multIdx[changingIndex] = 0;
                    changingIndex--;
                    // Return on overflow 
                    if(changingIndex >= degree()) { return; }
                    multIdx[changingIndex]++;
                }
            }
        }
        
        /// Creates a tensor with the given dimensions and uses the given random generator and distribution to assign the values to the entries.
        template<class generator, class distribution, ADD_MOVE(std::vector<size_t>, Vec)>
        static FullTensor construct_random(Vec&& _dimensions, generator& _rnd, distribution& _dist) {
            FullTensor result(std::forward<Vec>(_dimensions), DONT_SET_ZERO());
            for(size_t i=0; i < result.size; ++i) {
                result.data.get()[i] = _dist(_rnd);
            }
            return result;
        }
        
        
        // Unfortunaly all ALLOW_MOVE construcotrs based on vectors have to be copied for initializer_list
        
        /// Creates a tensor with the given dimensions and undefined entries.
        _inline_ explicit FullTensor(std::initializer_list<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) : FullTensor(std::vector<size_t>(_dimensions), DONT_SET_ZERO()) {}
        
        /// Creates a tensor with the given dimensions and all entries equals zero.
        _inline_ explicit FullTensor(std::initializer_list<size_t>&& _dimensions) : FullTensor(std::vector<size_t>(_dimensions)) {}
        
        /// Creates a tensor with the given dimensions and uses the given data as entries.
        ALLOW_MOVE(std::shared_ptr<value_t>, SPtr)
        _inline_ FullTensor(std::initializer_list<size_t>&& _dimensions, SPtr&& _data) : FullTensor(std::vector<size_t>(_dimensions), std::forward<SPtr>(_data)) {}
        
        /// Creates a tensor with the given dimensions and uses the given data as entries.
        _inline_ FullTensor(std::initializer_list<size_t>&& _dimensions, std::unique_ptr<value_t[]>&& _data) : FullTensor(std::vector<size_t>(_dimensions), std::move(_data)) {}
        
        /// Creates a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
        _inline_ explicit FullTensor(std::initializer_list<size_t>&& _dimensions, const std::function<value_t()>& _f)  : FullTensor(std::vector<size_t>(_dimensions), _f) {}
        
        /// Creates a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
        _inline_ explicit FullTensor(std::initializer_list<size_t>&& _dimensions, const std::function<value_t(const size_t)>& _f) : FullTensor(std::vector<size_t>(_dimensions), _f) {}
            
        /// Creates a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
        _inline_ explicit FullTensor(std::initializer_list<size_t>&& _dimensions, const std::function<value_t(const std::vector<size_t>&)>& _f)  : FullTensor(std::vector<size_t>(_dimensions), _f) {}
        
        /// Creates a tensor with the given dimensions and uses the given random generator and distribution to assign the values to the entries.
        template<class generator, class distribution>
        _inline_ static FullTensor construct_random(std::initializer_list<size_t>&& _dimensions, generator& _rnd, distribution& _dist) {
            return construct_random(std::vector<size_t>(std::move(_dimensions)), _rnd, _dist);
        }
        
        
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Virtual "Constructors" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Returns a pointer to a copy of the FullTensor.
        virtual Tensor* get_copy() const override;

        /// Returns a pointer to a moved copy of the FullTensor.
        virtual Tensor* get_moved_copy() override;
        
        /// Returns a pointer to a newly constructed order zero tensor of appropriate type with entry equals zero
        
        /// Returns a pointer to a newly constructed FullTensor.
        virtual Tensor* construct_new() const override;
        
        /// Returns a pointer to a newly constructed FullTensor with given dimensions and all entries set to zero.
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions) const override;
        
        /// Returns a pointer to a newly constructed FullTensor with given dimensions and all entries set to zero.
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions) const override;
        
        /// Returns a pointer to a newly constructed FullTensor with given dimensions and undefined entries.
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) const override;
        
        /// Returns a pointer to a newly constructed FullTensor with given dimensions and undefined entries.
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) const override;

        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Ensures that this tensor is the sole owner of its data. If needed new space is allocated and all entries are copied.
        virtual void ensure_own_data() override;
        
        /// Ensures that this tensor is the sole owner of its data space. If needed new space is allocated with entries left undefined.
        virtual void ensure_own_data_no_copy() override;
        
        /// Checks whether there is a non-trivial scaling factor and applies it if nessecary.
        virtual void apply_factor() override;
        
        /// Checks whether there is a non-trivial factor and applies it. Even if no factor is applied ensure_own_data() is called.
        virtual void ensure_own_data_and_apply_factor() override;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        FullTensor& operator=(const FullTensor&  _other);
        FullTensor& operator=(      FullTensor&& _other);
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        FullTensor& operator+=(const FullTensor& _other);
        FullTensor  operator+( const FullTensor& _other) const;
        
        FullTensor& operator-=(const FullTensor& _other);
        FullTensor  operator-( const FullTensor& _other) const;
        
        FullTensor& operator*=(const value_t _prod);
        FullTensor  operator*( const value_t _prod) const;
        
        FullTensor& operator/=(const value_t _div);
        FullTensor  operator/( const value_t _div) const;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics with SparseTensors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        FullTensor& operator+=(const SparseTensor& _other);
        FullTensor  operator+( const SparseTensor& _other) const;
        
        FullTensor& operator-=(const SparseTensor& _other);
        FullTensor  operator-( const SparseTensor& _other) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        virtual value_t& operator[](const size_t _position) override;
        virtual value_t operator[](const size_t _position) const override;
        
        virtual value_t& operator[](const std::vector<size_t>& _positions) override;
        virtual value_t operator[](const std::vector<size_t>& _positions) const override;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Modifiers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Resets the tensor with the given dimensions and undefined entries.
        virtual void reset(const std::vector<size_t>&  _newDim, _unused_ DONT_SET_ZERO) override;
        
        /// Resets the tensor with the given dimensions and undefined entries.
        virtual void reset(      std::vector<size_t>&& _newDim, _unused_ DONT_SET_ZERO) override;
        
        /// Resets the tensor with the given dimensions and all entries equals zero.
        virtual void reset(const std::vector<size_t>&  _newDim) override;
        
        /// Resets the tensor with the given dimensions and all entries equals zero.
        virtual void reset(      std::vector<size_t>&& _newDim) override;
        
        void resize_dimension(const size_t _n, const size_t _newDim, size_t _cutPos=~0ul);
        
        void remove_slate(uint _indexNb, uint _pos);
        
        void modify_diag_elements(const std::function<void(value_t&)>& _f);
            
        void modify_diag_elements(const std::function<void(value_t&, const size_t)>& _f);

        void modify_elements(const std::function<void(value_t&)>& _f);
        
        void modify_elements(const std::function<void(value_t&, const size_t)>& _f);

        void modify_elements(const std::function<void(value_t&, const std::vector<size_t>&)>& _f);
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Higher functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Checks whether the object is of type SparseTensor.
        virtual bool is_sparse() const override;
        
        
        /// Returns the number of non-zero entries.
        virtual size_t count_non_zero_entries(const value_t _eps = 1e-14) const override;
        
        /// Calculates the frobenious norm of the tensor.
        virtual value_t frob_norm() const override;
		
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Returns a string representation of the Tensor
        virtual std::string to_string() const override;
        
        /// Compares the Tensor entriewise to the given data.
        virtual bool compare_to_data(std::vector<value_t> _values, const double _eps = 1e-14) const override;
        
        /// Compares the Tensor entriewise to the given data.
        virtual bool compare_to_data(const value_t* _values, const double _eps = 1e-14) const override;
        
    };
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Other Direction arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    _inline_ FullTensor operator*(const value_t _lhs, const FullTensor& _rhs) { return _rhs*_lhs; }
    
    /// Returns the frobenius norm of the given tensor
    _inline_ value_t frob_norm(const FullTensor& _fullTensor) { return _fullTensor.frob_norm(); }
    
    /// Checks whether two FullTensor are approximately equal using either the frobenious norm of the difference (default) or by checking whether all entries are approximately equal.
    bool approx_equal(const xerus::FullTensor& _a, const xerus::FullTensor& _b, const xerus::value_t _eps = 1e-14, const bool pureDataCompare = false);
}


