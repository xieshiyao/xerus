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

#include "unknownContraction.h"

namespace xerus {

    std::shared_ptr<Tensor> unknown_contraction(
            std::shared_ptr<Tensor> _lhs, std::vector<Index> _lhsIndices,
            std::shared_ptr<Tensor> _rhs, std::vector<Index> _rhsIndices,
            std::vector<Index> _resultIndices) 
    {
            // only fullTensor format implemented so far
            std::shared_ptr<FullTensor> lhsFull = std::dynamic_pointer_cast<FullTensor>(_lhs);
            std::shared_ptr<FullTensor> rhsFull = std::dynamic_pointer_cast<FullTensor>(_rhs);
            if (lhsFull && rhsFull) {
                    std::shared_ptr<FullTensor> result(new FullTensor((uint)_resultIndices.size()));
                    (*result)(_resultIndices) = contract((*lhsFull)(_lhsIndices), (*rhsFull)(_rhsIndices));
                    return result;
            }
            LOG(fatal, "unknown_contraction for given types not yet implemented");
            exit(1);
    }
}
