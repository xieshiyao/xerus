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

#include <xerus/indexedTensorList.h>

namespace xerus {
    IndexedTensorList::IndexedTensorList(IndexedTensorList&& _old) : tensors(std::move(_old.tensors)) { } 
    
    IndexedTensorList::IndexedTensorList(const IndexedTensorWritable<Tensor>& _first, const IndexedTensorWritable<Tensor>& _second) {
        tensors.emplace_back(&_first);
        tensors.emplace_back(&_second);
    }
    
    void IndexedTensorList::operator=(std::function<void(const std::vector<const IndexedTensorWritable<Tensor>*>&)> _f) const {
        _f(tensors);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - External functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    IndexedTensorList operator,(const IndexedTensorWritable<Tensor>& _first, const IndexedTensorWritable<Tensor>& _second) {
        return IndexedTensorList(_first, _second);
    }

    IndexedTensorList operator,(IndexedTensorList &&_first, const IndexedTensorWritable<Tensor> &_second) {
        _first.tensors.emplace_back(&_second); // Hope this is standardconform. maybe we have to move-construct a new object
        return std::move(_first);
    }
}
