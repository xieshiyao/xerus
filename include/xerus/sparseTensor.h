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
//         std::shared_ptr<std::map<size_t, value_t>> sparseData;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

         /// @brief Empty constructor, which creates an order zero tensor with no entries set (i.e. which is completely zero).
        explicit SparseTensor();
        
        /// @brief Copy Constructors
        implicit SparseTensor( const Tensor&  _other);
        
        /// @brief Move Constructors
		implicit SparseTensor(       Tensor&& _other);
        
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
                REQUIRE(!misc::contains(*sparseData, entry.first), "Allready contained " << entry.first);
                sparseData->insert(std::move(entry));
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
                REQUIRE(!misc::contains(*sparseData, entry.first), "Allready contained " << entry.first);
                sparseData->insert(std::move(entry));
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
            while(result.sparseData->size() < _n) {
                result.sparseData->emplace(entryDist(_rnd), _dist(_rnd));
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
        
        
    };
    
    
} 
