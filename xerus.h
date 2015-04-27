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

#define MISC_NAMESPACE xerus

// The misc stuff needed by xerus
#include "misc/standard.h"
#include "misc/namedLogger.h"
#include "misc/missingFunctions.h"
#include "misc/stringUtilities.h"
#include "misc/sfinae.h"
#include "misc/blasLapackWrapper.h"
#include "misc/selectedFunctions.h"
#include "misc/callStack.h"
#include "misc/simpleNumerics.h"

// Additionall std stuff needed by xerus
#include <atomic>


// File which sets the custom log levels
#include "tensorLogger.h"

   
namespace xerus {
    /// The type of values to be used by xerus. In future versions this should be allowed to be float, double, or complex.
    typedef double value_t;
    
    namespace internal {
        /// Internal deleter function, needed because std::shared_ptr misses an array overload.
        void array_deleter_vt(value_t* const _toDelete);
        
        /// Internal deleter functions, needed because std::shared_ptr misses an array overload.
        void array_deleter_st(size_t* const _toDelete);
    }
    
    /// Helper class to provide possible overloads of several Tensor constructors
    class DONT_SET_ZERO {};
}


// All the xerus headers
#include "index.h"
#include "assignedIndices.h"
#include "indexedTensorReadOnly.h"
#include "indexedTensorWritable.h"
#include "indexedTensor.h"
#include "indexedTensorMoveable.h"
#include "indexedTensorList.h"
#include "tensor.h"
#include "fullTensor.h"
#include "sparseTensor.h"
#include "sparseTensor_contraction.h"
#include "indexedTensor_tensor.h"
#include "indexedTensor_tensor_operators.h"
#include "indexedTensor_tensor_factorisations.h"
#include "tensorNode.h"
#include "tensorNetwork.h"
#include "indexedTensor_TN.h"
#include "indexedTensor_TN_operators.h"
#include "contractionHeuristic.h"
#include "ttTensor.h"
// #include "ttTensor_specializations.h"
#include "algorithm/als.h"

