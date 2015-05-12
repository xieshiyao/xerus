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

#include "fullTensor.h"

namespace xerus {
    
    /// The main class used to represent any sparse tensor.
    class SparseTensor final : public Tensor {
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        std::shared_ptr<std::map<size_t, value_t>> entries;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

        ///Empty constructor (creates an order zero Tensor)
        explicit SparseTensor();
        
        ///Copy Constructors
        implicit SparseTensor( const SparseTensor&  _other);
        
        /// Move Constructors
        implicit SparseTensor(       SparseTensor&& _other);
        
        /// Creates a SparseTensor with the given dimensions and all entries equals zero (i.e. no entries)
        ALLOW_MOVE(std::vector<size_t>, T)
        explicit SparseTensor(T&& _dimensions) : Tensor(std::forward<T>(_dimensions)), entries(new std::map<size_t, value_t>()) { }
        
        /// Creates a SparseTensor with the given dimensions and all entries equals zero (i.e. no entries)
        explicit SparseTensor(std::initializer_list<size_t>&& _dimensions);
        
        /// Creates a SparseTensor from the given FullTensor, keeps all values larger than eps
        explicit SparseTensor(const FullTensor & _full, const double _eps = 1e-14);
        
        /// Creates a SparseTensor with the given dimensions and uses the given function to create N non zero entries
        ALLOW_MOVE(std::vector<size_t>, Vec)
        explicit SparseTensor(Vec&& _dimensions, std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f, const size_t _N) : SparseTensor(std::forward<Vec>(_dimensions)) {
            REQUIRE(_N <= size, "Cannot create more non zero entries that the dimension of the Tensor.");
            for (size_t i=0; i < _N; ++i) {
                entries->insert(_f(i, size));
            }
        }
        
        explicit SparseTensor(std::initializer_list<size_t>&& _dimensions, std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f, const size_t _N) : SparseTensor(std::move(_dimensions)) {
            REQUIRE(_N <= size, "Cannot create more non zero entries that the dimension of the Tensor.");
            for (size_t i=0; i < _N; ++i) {
                std::pair<size_t, value_t> entry = _f(i, size);
                REQUIRE(entry.first < size, "Postion is out of bounds " << entry.first);
                REQUIRE(!contains(*entries, entry.first), "Allready contained " << entry.first);
                entries->insert(std::move(entry));
            } 
            
            REQUIRE(count_non_zero_entries() == _N , "Oo " << _N << " != " << count_non_zero_entries());
        }
        
        /// Creates a SparseTensor with the given dimensions and uses the given random generator and distribution to assign the values to n randomly choosen entries.
        template<ADD_MOVE(std::vector<size_t>, Vec), class generator, class distribution>
        static SparseTensor construct_random(Vec&& _dimensions, const size_t _n, generator& _rnd, distribution& _dist) {
            SparseTensor result(std::forward<Vec>(_dimensions));
            REQUIRE(_n < result.size, " Cannot create " << _n << " non zero entries in a tensor with only " << result.size << " total entries!");
            
            std::uniform_int_distribution<size_t> entryDist(0, result.size-1);
            while(result.entries->size() < _n) {
                result.entries->emplace(entryDist(_rnd), _dist(_rnd));
            }

            return result;
        }
        
        /// Creates a SparseTensor with the given dimensions and uses the given random generator and distribution to assign the values to n randomly choosen entries.
        template<class generator, class distribution>
        static SparseTensor construct_random(std::initializer_list<size_t>&& _dimensions, const size_t _n, generator& _rnd, distribution& _dist) {
            return construct_random(std::vector<size_t>(_dimensions), _n, _rnd, _dist);
        }
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Virtual "Constructors" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Returns a pointer containing a copy of the object with appropriate type
        virtual Tensor* get_copy() const;
        
        /// Returns a pointer containing a moved copy of the object with appropriate type
        virtual Tensor* get_moved_copy();
        
        /// Returns a pointer to a newly constructed order zero tensor of appropriate type with entry equals zero
        virtual Tensor* construct_new() const override;
        
        /// Returns a pointer to a newly constructed tensor of appropriate type with all entries set to zero
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions) const override;
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions) const override;
        
        /// Returns a pointer to a newly constructed tensor of appropriate type with undefined entries
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) const override;
        virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) const override;
        
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Ensures that this tensor is the sole owner of the entrie set. If needed a new set with the same entries is created
        virtual void ensure_own_data();
        
        /// Ensures that this tensor is the sole owner of the entrie set. If needed a new set with undefined entries is created
        virtual void ensure_own_data_no_copy();
        
        /// Checks whether there is a non-trivial factor and applies it.
        virtual void apply_factor();
        
        /// Checks whether there is a non-trivial factor and applies it. Even if no factor is applied ensure_own_data() is called.
        virtual void ensure_own_data_and_apply_factor();
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        SparseTensor& operator=(const SparseTensor& _other);
        SparseTensor& operator=(SparseTensor&& _other);
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Cast to FullTensor - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        operator FullTensor() const;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        virtual value_t& operator[](const size_t _position) override;
        virtual value_t  operator[](const size_t _position) const override;
        
        virtual value_t& operator[](const std::vector<size_t>& _positions) override;
        virtual value_t operator[](const std::vector<size_t>& _positions) const override;
        
        value_t at(const size_t _position) const;
        value_t at(const std::vector<size_t>& _indices) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        SparseTensor& operator+=(const SparseTensor& _other);
        SparseTensor  operator+(const SparseTensor& _other) const;
        
        SparseTensor& operator-=(const SparseTensor& _other);
        SparseTensor operator-(const SparseTensor& _other) const;
        
        SparseTensor& operator*=(const value_t _prod);
        SparseTensor operator*(const value_t _prod) const;
        
        SparseTensor& operator/=(const value_t _div);
        SparseTensor operator/(const value_t _div) const;
        
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Higher functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        /// Resets the tensor with the given dimensions and undefined entries
        virtual void reset(const std::vector<size_t>&  _newDim, _unused_ DONT_SET_ZERO);
        virtual void reset(      std::vector<size_t>&& _newDim, _unused_ DONT_SET_ZERO);
        
        /// Resets the tensor with the given dimensions and all entries equals zero
        virtual void reset(const std::vector<size_t>&  _newDim);
        virtual void reset(      std::vector<size_t>&& _newDim);
        
        /// Returns whether the object is of type SparseTensor
        virtual bool is_sparse() const;
        
        /// Returns the number of non-zero entries
        virtual size_t count_non_zero_entries(const value_t _eps = 1e-14) const override;
        
        /// Returns the frobenious norm of the tensor
        virtual value_t frob_norm() const;
        
        /// Returns a string representation of the Tensor
        virtual std::string to_string() const override;
        
        
        /// Compares the Tensor entriewise to the given data
        virtual bool compare_to_data(std::vector<value_t> _values, const double _eps = 1e-14) const override;
        
        /// Compares the Tensor entriewise to the given data
        virtual bool compare_to_data(const value_t* _values, const double _eps = 1e-14) const override;
    };
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Other Direction arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    _inline_ SparseTensor operator*(const value_t _lhs, const SparseTensor& _rhs) {
        return _rhs*_lhs;
    }
    
} 

bool approx_equal(const xerus::SparseTensor& _a, const xerus::FullTensor& _b, const xerus::value_t _eps, const bool pureDataCompare = false);
bool approx_equal(const xerus::FullTensor& _a, const xerus::SparseTensor& _b, const xerus::value_t _eps, const bool pureDataCompare = false);
