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

// This is the main .cpp file used for the creation of the xerus library. All pre compiled classes are included in this file, apart from the ones coming from includes.*

#include "xerus.h"

#include <initializer_list>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <memory>
#include <climits>

// Provide an implementation of the internal deleter functions
namespace xerus {
    namespace internal {
        void array_deleter_vt(value_t* const _toDelete) { delete[] _toDelete; }
        void array_deleter_st(size_t* const _toDelete) { delete[] _toDelete; }
    }
}

#include "index.hpp"
#include "tensor.hpp"
#include "tensor_specializations.hpp"
#include "fullTensor.hpp"
#include "sparseTensor.hpp"
#include "sparseTensor_contraction.hpp"
#include "indexedTensor_tensor_evaluate.hpp"
#include "indexedTensor_tensor_contraction.hpp"
#include "indexedTensor_tensor_solve.hpp"
#include "indexedTensor_tensor_operators.hpp"
#include "indexedTensor_tensor_factorisations.hpp"
#include "tensorNetwork_specializations.hpp"
#include "tensorNetwork.hpp"
#include "tensorNode.hpp"
#include "ttTensor_specializations.h"
#include "indexedTensor_TN_operators.hpp"
#include "contractionHeuristic.hpp"
#include "algorithm/als.hpp"

namespace xerus {
    template class TTNetwork<false>;
    template class TTNetwork<true>;
}
