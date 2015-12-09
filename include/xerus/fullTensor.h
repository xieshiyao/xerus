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
 * @brief Header file for the FullTensor class.
 */

#pragma once

#include "basic.h"
#include <memory>
#include "tensor.h"
#include "selectedFunctions.h"
#include "misc/sfinae.h"
#include "misc/performanceAnalysis.h"

namespace xerus {
    class SparseTensor;
    class TensorNetwork;
    
    /// @brief The xerus class used to represent all dense tensor objects.
    class FullTensor final : public Tensor {
    protected:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /** 
		 * @brief Shared pointer to the data array with size "size". 
		 * @details The data is stored such that indices increase from right to left (row-major order). 
		 * If the tensor is modified and not sole owner a deep copy is performed.
		 */
//         std::shared_ptr<value_t> denseData;
	public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// @brief Empty constructor, which creates an order zero tensor with zero as single entry.
        explicit FullTensor();
        
        /// @brief Copy constructor.
        implicit FullTensor(const FullTensor & _other);
        
        /// @brief Move constructor.
        implicit FullTensor(      FullTensor&& _other);
        
		
		/** 
		 * @brief Constructs a FullTensor from another Tensor.
		 * @param _other the Tensor which shall be used.
		 */
         FullTensor(const Tensor&  _other);
		
		/** 
		 * @brief Move-constructs a FullTensor from another Tensor.
		 * @param _other the Tensor which shall be used.
		 */
         FullTensor(      Tensor&& _other);
		
		
		/** 
		 * @brief Constructs a FullTensor by completely contracting a Tensor network.
		 * @param _other the TensorNetwork which shall be used.
		 */
        explicit FullTensor(const TensorNetwork&  _other);
		
        
		
		/** 
		 * @brief Constructs a FullTensor with the given degree, with all dimensions euqals one and zero as the single entry.
		 * @param _degree the future degree of the Tensor.
		 */
        explicit FullTensor(const size_t _degree);
        
		/** 
		 * @brief Constructs a FullTensor with the given dimensions and undefined entries.
		 * @details The second parameter is a DONT_SET_ZERO helper object that is only used to provide the function overload.
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit FullTensor(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO);
        
		/** 
		 * @brief Constructs a FullTensor with the given dimensions and undefined entries.
		 * @details The second parameter is a DONT_SET_ZERO helper object that is only used to provide the function overload.
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit FullTensor(      std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO);
        
		/** 
		 * @brief Constructs a FullTensor with the given dimensions and all entries equals zero.
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit FullTensor(const std::vector<size_t>&  _dimensions);
        
		/** 
		 * @brief Constructs a FullTensor with the given dimensions and all entries equals zero.
		 * @param _dimensions the future dimensions of the Tensor.
		 */
        explicit FullTensor(      std::vector<size_t>&& _dimensions);
        
		/** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the provided data as entries.
		 * @details The FullTensor will NOT modify the provided data in any way, unsless it is the sole owner (i.e. all other shared pointers are deleted).
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _data a shared_ptr to the data the FullTensor shall use. This must be (at least) of the size determined by the dimensions.
		 */
        template<ADD_MOVE(Vec, std::vector<size_t>), ADD_MOVE(SPtr, std::shared_ptr<value_t>)>
        explicit FullTensor(Vec&& _dimensions, SPtr&& _data) : Tensor(std::forward<Vec>(_dimensions), _data) { }
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the provided data as entries.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _data a unique_ptr to the data the FullTensor shall use. This must be (at least) of the size determined by the dimensions.
		 */
        explicit FullTensor(const std::vector<size_t> & _dimensions, std::unique_ptr<value_t[]>&& _data);
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the provided data as entries.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _data a unique_ptr to the data the FullTensor shall use. This must be (at least) of the size determined by the dimensions.
		 */
		explicit FullTensor(      std::vector<size_t>&& _dimensions, std::unique_ptr<value_t[]>&& _data);
        
        
		/** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details In this overload no value is passed to _f, i.e. _f must determine the values of the entries independend of their position,
		 * or keep track of the position itself. _f may assume that it is called for the entries in the order they are stored (i.e. row-major order)
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to use to set the entries of the FullTensor. 
		 */
        ALLOW_MOVE(Vec, std::vector<size_t>)
        explicit FullTensor(Vec&& _dimensions, const std::function<value_t()>& _f) : FullTensor(std::forward<Vec>(_dimensions), DONT_SET_ZERO()) {
            value_t* realData = denseData.get();
            for (size_t i=0; i < size; ++i) {
                realData[i] = _f();
            }
        }
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details In this overload the position of each entry assuming row-major order is passed to  _f.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to use to set the entries of the FullTensor. 
		 */
        ALLOW_MOVE(Vec, std::vector<size_t>)
        explicit FullTensor(Vec&& _dimensions, const std::function<value_t(const size_t)>& _f) : FullTensor(std::forward<Vec>(_dimensions), DONT_SET_ZERO()) {
            value_t* realData = denseData.get();
            for (size_t i=0; i < size; ++i) {
                realData[i] = _f(i);
            }
        }
        
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details In this overload the complete position of each entry is passed to  _f.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to use to set the entries of the FullTensor. 
		 */
        ALLOW_MOVE(Vec, std::vector<size_t>)
        explicit FullTensor(Vec&& _dimensions, const std::function<value_t(const std::vector<size_t>&)>& _f) : FullTensor(std::forward<Vec>(_dimensions), DONT_SET_ZERO()) {
            value_t* realData = denseData.get();
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
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given random generator and distribution to assign the values to the entries.
		 * @details The entries are assigned in the order they are stored (i.e. row-major order). Each assigned is a seperate call to the random distribution.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _rnd the random generator to be use.
		 * @param _dist the random distribution to be used.
		 */
        template<class generator, class distribution, ADD_MOVE(Vec, std::vector<size_t>)>
        static FullTensor random(Vec&& _dimensions, generator& _rnd, distribution& _dist) {
            FullTensor result(std::forward<Vec>(_dimensions), DONT_SET_ZERO());
			PA_START;
            for(size_t i=0; i < result.size; ++i) {
                result.denseData.get()[i] = _dist(_rnd);
            }
			PA_END("Random construction", "FullTensor", misc::to_string(result.size));
            return result;
        }
        
        template<class generator, class distribution, ADD_MOVE(Vec, std::vector<size_t>)>
        _deprecated_ static FullTensor construct_random(Vec&& _dimensions, generator& _rnd, distribution& _dist) {
            return FullTensor::random(std::forward<Vec>(_dimensions), _rnd, _dist);
        }
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given random generator and distribution to assign the values to the entries.
		 * @details See the std::vector variant for details.
		 */
        template<class generator, class distribution>
        _inline_ static FullTensor random(std::initializer_list<size_t>&& _dimensions, generator& _rnd, distribution& _dist) {
            return FullTensor::random(std::vector<size_t>(std::move(_dimensions)), _rnd, _dist);
        }
        
        template<class generator, class distribution>
        _deprecated_ _inline_ static FullTensor construct_random(std::initializer_list<size_t>&& _dimensions, generator& _rnd, distribution& _dist) {
            return FullTensor::random(std::vector<size_t>(std::move(_dimensions)), _rnd, _dist);
        }
        
        /** 
		 * @brief Creates a tensor with the given dimensions and undefined entries.
		 * @details See the std::vector variant for details.
		 */
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) : FullTensor(std::vector<size_t>(_dimensions), DONT_SET_ZERO()) {}
        
        
        /// Creates a tensor with the given dimensions and all entries equals zero. See the std::vector variant for more details.
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and all entries equals zero.
		 * @details See the std::vector variant for details.
		 */
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions) : FullTensor(std::vector<size_t>(_dimensions)) {}
        
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the provided data as entries.
		 * @details See the std::vector variant for details.
		 */
        ALLOW_MOVE(SPtr, std::shared_ptr<value_t>)
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions, SPtr&& _data) : FullTensor(std::vector<size_t>(_dimensions), std::forward<SPtr>(_data)) {}
        
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the provided data as entries.
		 * @details See the std::vector variant for details.
		 */
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions, std::unique_ptr<value_t[]>&& _data) : FullTensor(std::vector<size_t>(_dimensions), std::move(_data)) {}
        
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details See the std::vector variant for details.
		 */
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions, const std::function<value_t()>& _f)  : FullTensor(std::vector<size_t>(_dimensions), _f) {}
        
        
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details See the std::vector variant for details.
		 */
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions, const std::function<value_t(const size_t)>& _f) : FullTensor(std::vector<size_t>(_dimensions), _f) {}
            
            
        /** 
		 * @brief Constructs a FullTensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details See the std::vector variant for details.
		 */
        explicit FullTensor(std::initializer_list<size_t>&& _dimensions, const std::function<value_t(const std::vector<size_t>&)>& _f)  : FullTensor(std::vector<size_t>(_dimensions), _f) {}
        
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /** 
		 * @brief Standard assignment operator.
		 * @param _other the FullTensor to be assinged to this one.
		 * @return a reference to this FullTensor.
		 */
        FullTensor& operator=(const FullTensor&  _other);
        
		/** 
		 * @brief Standard move-assignment operator.
		 * @param _other the FullTensor to be move-assinged to this one.
		 * @return a reference to this FullTensor.
		 */
        FullTensor& operator=(      FullTensor&& _other);
		
		/** 
		 * @brief Assignment operator for a general tensor.
		 * @param _other the Tensor to be assinged to this one.
		 * @return a reference to this FullTensor.
		 */
        FullTensor& operator=(const Tensor&  _other);
        
		/** 
		 * @brief Move-assignment operator for a general tensor.
		 * @param _other the Tensor to be move-assinged to this one.
		 * @return a reference to this FullTensor.
		 */
        FullTensor& operator=(      Tensor&& _other);
        
		
    };
    
    
    
	
}
