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

#include "xerus.h"

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    Tensor::Tensor() : size(1), factor(1.0) {}

    Tensor::Tensor(const  Tensor&  _other ) : dimensions(_other.dimensions), size(_other.size), factor(_other.factor) { }
    
    Tensor::Tensor(       Tensor&& _other ) : dimensions(_other.dimensions), size(_other.size), factor(_other.factor) { }
    
    Tensor::~Tensor() {}
    
    void Tensor::assign(const Tensor& _other) {
        dimensions = _other.dimensions;
        size = _other.size;
        factor = _other.factor;
    }
    
    void Tensor::assign(Tensor&& _other) {
        dimensions = std::move(_other.dimensions);
        size = _other.size;
        factor = _other.factor;
    }
}
