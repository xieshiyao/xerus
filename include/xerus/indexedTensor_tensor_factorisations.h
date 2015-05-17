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

#include "indexedTensorWritable.h"

namespace xerus {

    // Helper class for IndexedTensorLists
    class SVD {
    public:
        const IndexedTensorReadOnly<Tensor>& input;
        const double epsilon;
        SVD(const IndexedTensorReadOnly<Tensor>& _input, const double _epsilon = 1e-14) : input(_input), epsilon(_epsilon) { }
        
        void operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const ;
    };

    class QR {
    public:
        const IndexedTensorReadOnly<Tensor>* input;
        QR(const IndexedTensorReadOnly<Tensor>& _input) : input(&_input) { }
        
        void operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const;
    };

    class RQ {
    public:
        const IndexedTensorReadOnly<Tensor>* input;
        RQ(const IndexedTensorReadOnly<Tensor>& _input) : input(&_input) { }
        
        void operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const;
    };
    
    // TODO Split Core
}