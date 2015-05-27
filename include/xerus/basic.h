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

#include "misc/standard.h"
#include "tensorLogger.h"

namespace xerus {
    /// The type of values to be used by xerus. In future versions this should be allowed to be float, double, or complex.
    typedef double value_t;
    
    namespace internal {
        /// Internal deleter function, needed because std::shared_ptr misses an array overload.
        void array_deleter_vt(value_t* const _toDelete);
        
        /// Internal deleter function, needed because std::shared_ptr misses an array overload.
        void array_deleter_st(size_t* const _toDelete);
    }
    
    /// Helper class to provide possible overloads of several Tensor constructors.
    class DONT_SET_ZERO {};
}
