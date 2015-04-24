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
    
    /// The Index class is used to write indexed tensor expressen, e.g. A(i,j)*B(j,k). Here i,j,k are of type xerus::Index. 
    class Index {
    private:
        static std::atomic<long> idThreadInitCounter;
        static thread_local long idCounter;
        
    public:
        /// If negative the unique Id of the index. If positive, the index is fixed (e.g. = 5).
        long valueId;
        
        /// The span allows a single index to cover more than one dimension of a tensor.
        size_t span;
        
        /// Flag indicating that this index covers all but span dimensions of the tensor.
        bool inverseSpan;
        
        /// Empty constructor that creates a new Index with new ID. Use this to create indices.
        Index();
        
        Index(const Index&) = default;
        Index(Index&&) = default;
        
        /// Integers are implicitly allowed to be casted to Index, to allow expression as A(i) = B(3,i), i.e. A is the third row of B.
        implicit Index(const long _i);
        
        /// Internal constructor to create indices with larger span. Do not use this unless you know what you are doing.
        explicit Index(const long _valueId, const size_t _span, const bool _inverseSpan = false);
        
    public:
        /// Indices are default assignable
        Index& operator=(const Index&) = default;
        
        /// Indices are default moveable
        Index& operator=(Index&&) = default;
        
        /// Checks whether the Index represents a fixed number
        _inline_ bool is_fixed() const {
            return valueId >= 0;
        }
        
        /// Allow the creation of Indices covering more than one dimension using the power operator. E.g. A(i^2) = B(i^2) + C(i^2), defines A as the entriewise sum of the matrices B and C.
        _inline_ Index operator^(const size_t _span) const {
            return Index(valueId, _span);
        }
        
        /** Allow the creation of Indices covering all but x dimensions using the and operator. 
         *  E.g. A() = B(i&0) * C(i&0), defines A as the full contraction between B and C, indifferent of the order of B and C (which have to coincide however). 
         */
        _inline_ Index operator&(const size_t _span) const {
            return Index(valueId, _span, true);
        }
        
        /// Two Indices are equal if their valueId coincides.
        _inline_ bool operator==(const Index& _other) const {
            return valueId == _other.valueId;
        }
        
        /// Two Indices are equal if their valueId coincides.
        _inline_ bool operator!=(const Index& _other) const {
            return valueId != _other.valueId;
        }
        
        /// The Comparision operator is needed for indices to be orderable in std::set, the valueId is used.
        _inline_ bool operator<(const Index& _other) const {
            return valueId < _other.valueId;
        }
    };
    
    /// Allows to pretty print Indices, giving the valueId and span.
    std::ostream& operator<<(std::ostream& _out, const xerus::Index& _idx); 
}
