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

#include <vector>
#include <functional>

namespace xerus {
    // Necessary forward declaritons
    template<class tensor_type> class IndexedTensorWritable;
    class Tensor;
        
    class IndexedTensorList {
    public:
        std::vector<const IndexedTensorWritable<Tensor>*> tensors;
        
        IndexedTensorList() = delete;
        IndexedTensorList(const IndexedTensorList& _old) = delete;
        
        IndexedTensorList(IndexedTensorList&& _old);
        
        IndexedTensorList(const IndexedTensorWritable<Tensor>& _first, const IndexedTensorWritable<Tensor>& _second);
        
        // Generic =operator 
        void operator=(std::function<void(const std::vector<const IndexedTensorWritable<Tensor>*>&)> _f) const;
    };

    IndexedTensorList operator,(const IndexedTensorWritable<Tensor>& _first, const IndexedTensorWritable<Tensor>& _second);

    IndexedTensorList operator,(IndexedTensorList &&_first, const IndexedTensorWritable<Tensor> &_second);
}
