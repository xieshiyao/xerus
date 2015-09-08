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
 * @brief Header file for the SparseTensor class.
 */

#pragma once

#include "fullTensor.h"
#include "misc/check.h"
#include "misc/missingFunctions.h"
#include "misc/performanceAnalysis.h"

namespace xerus {
    
    /// @brief The xerus class used to represent all sparse tensor objects.
    class SparseTensor final : public Tensor {
    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Shared pointer to the a map containing the entries of the SparseTensor. 
		 * @details The entries are stored in a map which uses the position of each entry assuming row-major ordering as key value.
		 * If the tensor is modified and not sole owner a deep copy is performed.
		 */
        std::shared_ptr<std::map<size_t, value_t>> entries;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

         /// @brief Empty constructor, which creates an order zero tensor with no entries set (i.e. which is completely zero).
        explicit SparseTensor();
        
        /// @brief Copy Constructors
        implicit SparseTensor( const SparseTensor&  _other);
        
        /// @brief Move Constructors
        implicit SparseTensor(       SparseTensor&& _other);
        
		/** 
		 * @brief Constructs a SparseTensor with the given dimensions and no set entries (i.e. completely zero).
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit SparseTensor(const std::vector<size_t> & _dimensions);
        
		/** 
		 * @brief Constructs a SparseTensor with the given dimensions and no set entries (i.e. completely zero).
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit SparseTensor(      std::vector<size_t>&& _dimensions);
        
		/** 
		 * @brief Constructs a SparseTensor with the given dimensions and no set entries (i.e. completely zero).
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit SparseTensor(std::initializer_list<size_t>&& _dimensions);
        
		/** 
		 * @brief Constructs a SparseTensor from the given FullTensor.
		 * @details Every entry of @a _full that is larger than @a _eps is used, all other are truncated to zero.
		 * @param _full the FullTensor to be used.
		 * @param _eps (optional) the epsilon to be used to truncate the entries.
		 */
        explicit SparseTensor(const FullTensor & _full, const double _eps = EPSILON);
        
		/** 
		 * @brief Constructs a SparseTensor with the given dimensions and uses the given function @a _f to create @a _N non zero entries.
		 * @details @a _f is called with the current number of entries present and the number of possible entries (i.e. size). @a _f shall return a pair containg the position
		 * and value of the next entry. @a _f is required not to return a position twice.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to be used to create each non zero entry. 
		 * @param _N the number of non-zero entries to be created.
		 */
        ALLOW_MOVE(Vec, std::vector<size_t>)
        explicit SparseTensor(Vec&& _dimensions, std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f, const size_t _N) : SparseTensor(std::forward<Vec>(_dimensions)) {
            REQUIRE(_N <= size, "Cannot create more non zero entries that the dimension of the Tensor.");
            for (size_t i=0; i < _N; ++i) {
                std::pair<size_t, value_t> entry = _f(i, size);
                REQUIRE(entry.first < size, "Postion is out of bounds " << entry.first);
                REQUIRE(!misc::contains(*entries, entry.first), "Allready contained " << entry.first);
                entries->insert(std::move(entry));
            } 
            
            REQUIRE(count_non_zero_entries() == _N , "Oo " << _N << " != " << count_non_zero_entries());
        }
        
        /** 
		 * @brief Constructs a SparseTensor with the given dimensions and uses the given function @a _f to create @a _N non zero entries.
		 * @details @a _f is called with the current number of entries present and the number of possible entries (i.e. size). @a _f shall return a pair containg the position
		 * and value of the next entry. @a _f is required not to return a position twice.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to be used to create each non zero entry. 
		 * @param _N the number of non-zero entries to be created.
		 */
        explicit SparseTensor(std::initializer_list<size_t>&& _dimensions, std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f, const size_t _N) : SparseTensor(std::move(_dimensions)) {
            REQUIRE(_N <= size, "Cannot create more non zero entries that the dimension of the Tensor.");
            for (size_t i=0; i < _N; ++i) {
                std::pair<size_t, value_t> entry = _f(i, size);
                REQUIRE(entry.first < size, "Postion is out of bounds " << entry.first);
                REQUIRE(!misc::contains(*entries, entry.first), "Allready contained " << entry.first);
                entries->insert(std::move(entry));
            } 
            
            REQUIRE(count_non_zero_entries() == _N , "Oo " << _N << " != " << count_non_zero_entries());
        }
        
        /** 
		 * @brief Constructs a random SparseTensor with the given dimensions.
		 * @details The given random generator @a _rnd and distribution @a _dist are used to assign the values to @a _n randomly choosen entries.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _n the number of non-zero entries to be created.
		 * @param _rnd the random generator to be used.
		 * @param _dist the random distribution to be used.
		 */
        template<ADD_MOVE(Vec, std::vector<size_t>), class generator, class distribution>
        static SparseTensor random(Vec&& _dimensions, const size_t _n, generator& _rnd, distribution& _dist) {
            SparseTensor result(std::forward<Vec>(_dimensions));
            REQUIRE(_n <= result.size, " Cannot create " << _n << " non zero entries in a tensor with only " << result.size << " total entries!");
            
            std::uniform_int_distribution<size_t> entryDist(0, result.size-1);
			PA_START;
            while(result.entries->size() < _n) {
                result.entries->emplace(entryDist(_rnd), _dist(_rnd));
            }
			PA_END("Random construction", "SparseTensor", misc::to_string(_n)+"/"+misc::to_string(result.size));

            return result;
        }
        
        template<ADD_MOVE(Vec, std::vector<size_t>), class generator, class distribution>
        _deprecated_ static SparseTensor construct_random(Vec&& _dimensions, const size_t _n, generator& _rnd, distribution& _dist) {
			return SparseTensor::random(std::forward<Vec>(_dimensions), _n, _rnd, _dist);
		}
        
        /** 
		 * @brief Constructs a random SparseTensor with the given dimensions.
		 * @details The given random generator @a _rnd and distribution @a _dist are used to assign the values to @a _n randomly choosen entries.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _n the number of non-zero entries to be created.
		 * @param _rnd the random generator to be used.
		 * @param _dist the random distribution to be used.
		 */
        template<class generator, class distribution>
        static SparseTensor random(std::initializer_list<size_t>&& _dimensions, const size_t _n, generator& _rnd, distribution& _dist) {
            return SparseTensor::random(std::vector<size_t>(_dimensions), _n, _rnd, _dist);
        }
        
        template<class generator, class distribution>
        _deprecated_ static SparseTensor construct_random(std::initializer_list<size_t>&& _dimensions, const size_t _n, generator& _rnd, distribution& _dist) {
			return SparseTensor::random(std::vector<size_t>(_dimensions), _n, _rnd, _dist);
		}
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Virtual "Constructors" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        virtual Tensor* get_copy() const;
        
        virtual Tensor* get_moved_copy();
        
        virtual Tensor* construct_new() const override;
        
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions) const override;
        
		virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions) const override;
        
        virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) const override;

		virtual Tensor* construct_new(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) const override;
        
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        virtual void ensure_own_data();
        
        virtual void ensure_own_data_no_copy();
        
        virtual void apply_factor();
        
        virtual void ensure_own_data_and_apply_factor();
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Standard assignment operator.
		 * @param _other the SparseTensor to be assinged to this one.
		 * @return a reference to this SparseTensor.
		 */
        SparseTensor& operator=(const SparseTensor& _other);
		
		/** 
		 * @brief Standard move-assignment operator.
		 * @param _other the SparseTensor to be move-assinged to this one.
		 * @return a reference to this SparseTensor.
		 */
        SparseTensor& operator=(SparseTensor&& _other);
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Read/Write access a single entry.
		 * @details Note that the request for write access to this entry requires the SparseTensor to create it, i.e. to lower the sparsity. If only read access
		 * is required use at() instead.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return a reference to the selected entry.
		 */
        virtual value_t& operator[](const size_t _position) override;
		
        virtual value_t  operator[](const size_t _position) const override;
        
		/** 
		 * @brief Read/Write access a single entry.
		 * @details Note that the request for write access to this entry requires the SparseTensor to create it, i.e. to lower the sparsity. If only read access
		 * is required use at() instead.
		 * @param _position the position of the desired entry.
		 * @return a reference to the selected entry.
		 */
        virtual value_t& operator[](const std::vector<size_t>& _positions) override;
		
        virtual value_t operator[](const std::vector<size_t>& _positions) const override;
        
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return the selected entry.
		 */
        value_t at(const size_t _position) const;
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry.
		 * @return the selected entry.
		 */
        value_t at(const std::vector<size_t>& _indices) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Adds the @a _other SparseTensor entrywise to this one.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the SparseTensor to be added to this one.
		 * @return a reference to this SparseTensor.
		 */
        SparseTensor& operator+=(const SparseTensor& _other);
		
		/** 
		 * @brief Calculates the entrywise sum of this SparseTensor and @a _other.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the second summand.
		 * @return the sum as a SparseTensor.
		 */
        SparseTensor  operator+(const SparseTensor& _other) const;
        
		/** 
		 * @brief Subtracts the @a _other SparseTensor entrywise from this one.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the SparseTensor to be subtracted from this one.
		 * @return a reference to this SparseTensor.
		 */
        SparseTensor& operator-=(const SparseTensor& _other);
		
		/** 
		 * @brief Calculates the entrywise difference of this SparseTensor and @a _other.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the subtrahend.
		 * @return the difference as a SparseTensor.
		 */
        SparseTensor operator-(const SparseTensor& _other) const;
        
		/** 
		 * @brief Calculates the entrywise multiplication of this SparseTensor with a constant @a _factor.
		 * @details Internally this only results in a change in the global factor.
		 * @param _factor the factor,
		 * @return the resulting scaled SparseTensor.
		 */
        SparseTensor operator*(const value_t _factor) const;
        
		/** 
		 * @brief Calculates the entrywise divison of this SparseTensor by a constant @a _divisor.
		 * @details Internally this only results in a change in the global factor.
		 * @param _divisor the divisor,
		 * @return the resulting scaled SparseTensor.
		 */
        SparseTensor operator/(const value_t _divisor) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Higher functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		virtual void fix_slate(const size_t _dimension, const size_t _slatePosition);
        
        virtual void reset(const std::vector<size_t>&  _newDim, _unused_ DONT_SET_ZERO);
		
        virtual void reset(      std::vector<size_t>&& _newDim, _unused_ DONT_SET_ZERO);
        
        virtual void reset(const std::vector<size_t>&  _newDim);
		
        virtual void reset(      std::vector<size_t>&& _newDim);
        
        virtual bool is_sparse() const;
        
        virtual size_t count_non_zero_entries(const value_t _eps = EPSILON) const override;
        
		virtual bool all_entries_valid() const override;
		
        virtual value_t frob_norm() const;
        
        virtual std::string to_string() const override;
        
        virtual bool compare_to_data(const std::vector<value_t>& _values, const double _eps = EPSILON) const override;
        
        virtual bool compare_to_data(const value_t* _values, const double _eps = EPSILON) const override;
    };
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Other Direction arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	/** 
	* @brief Calculates the entrywise multiplication of the SparseTensor @a _lhs with a constant @a _rhs.
	* @details Internally this only results in a change in the global factor.
	* @param _lhs the SparseTensor that shall be scaled.
	* @param _rhs the factor to be used.
	* @return the resulting scaled SparseTensor.
	*/
    static _inline_ SparseTensor operator*(const value_t _lhs, const SparseTensor& _rhs) {
        return _rhs*_lhs;
    }
    
} 
