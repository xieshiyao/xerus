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
    public:
        #ifndef DISABLE_RUNTIME_CHECKS_
            enum Flag {FIXED, INVERSE_SPAN, FRACTIONAL_SPAN, OPEN, OPEN_CHECKED, NUM_FLAGS};
        #else
            enum Flag {FIXED, INVERSE_SPAN, FRACTIONAL_SPAN, OPEN, NUM_FLAGS};
        #endif
        
    private:
        static std::atomic<size_t> idThreadInitCounter;
        static thread_local size_t idCounter;
        
    public:
        /// Unqiue id of the index. In case the fixed flag is set, this is the fixed position.
        size_t valueId;
        
        /// The span states how many dimensions are covered by the index.
        size_t span;
        
        /// Bitset of all possible flags the index may possess.
        std::bitset<NUM_FLAGS> flags;
        
        /// Empty constructor that creates a new Index with new ID. Use this to create indices.
        Index();
        
        /// Indices are default copy constructable.
        Index(const Index&) = default;
        
        /// Indices are default move constructable.
        Index(Index&&) = default;
        
        /// Integers are implicitly allowed to be casted to Index, to allow expression as A(i) = B(3,i), i.e. A is the third row of B.
        implicit Index(const long _i);
        
        /// Internal constructor, do not use this unless you know what you are doing.
        explicit Index(const size_t _valueId, const size_t _span);
        
        /// Internal constructor, do not use this unless you know what you are doing.
        explicit Index(const size_t _valueId, const size_t _span, const Flag _flag1, const bool _flagValue1 = true);
        
        /// Internal constructor, do not use this unless you know what you are doing.
        explicit Index(const size_t _valueId, const size_t _span, const Flag _flag1, const Flag _flag2, const bool _flagValue1 = true, const bool _flagValue2 = true);
        
        /// Internal constructor, do not use this unless you know what you are doing.
        explicit Index(const size_t _valueId, const size_t _span, const Flag _flag1, const Flag _flag2, const Flag _flag3, const bool _flagValue1 = true, const bool _flagValue2 = true, const bool _flagValue3 = true);
        
    public:
        /// Indices are default assignable.
        Index& operator=(const Index&) = default;
        
        /// Indices are default moveable.
        Index& operator=(Index&&) = default;
        
        /// Checks whether the Index represents a fixed number.
        _inline_ bool fixed() const {
            return flags[Index::Flag::FIXED];
        }
        
        /// Checks whether the index is open.
        _inline_ bool open() const {
            REQUIRE(flags[Index::Flag::OPEN_CHECKED], "Check for index openness only allowed if the openness was checked before.");
            return flags[Index::Flag::OPEN];
        }
        
        /// Allow the creation of Indices covering more than one dimension using the power operator. E.g. A(i^2) = B(i^2) + C(i^2), defines A as the entriewise sum of the matrices B and C.
        _inline_ Index operator^(const size_t _span) const {
            REQUIRE(flags.none(), "Cannot apply ^ operator to an index that has any flag set.");
            return Index(valueId, _span);
        }
        
        /** Allow the creation of Indices covering all but x dimensions using the and operator. 
         *  E.g. A() = B(i&0) * C(i&0), defines A as the full contraction between B and C, indifferent of the order of B and C (which have to coincide however). 
         */
        _inline_ Index operator&(const size_t _span) const {
            REQUIRE(flags.none(), "Cannot apply & operator to an index that has any flag set.");
            return Index(valueId, _span, Flag::INVERSE_SPAN);
        }
        
        /// Two Indices are equal if their valueId coincides. Fixed indices are never equal.
        _inline_ bool operator==(const Index& _other) const {
            return valueId == _other.valueId && !fixed() && !_other.fixed();
        }
        
        /// Two Indices are equal if their valueId coincides. Fixed indices are never equal.
        _inline_ bool operator!=(const Index& _other) const {
            return valueId != _other.valueId || fixed() || _other.fixed();
        }
        
        /// The Comparision operator is needed for indices to be orderable in std::set, the valueId is used.
        _inline_ bool operator<(const Index& _other) const {
            return valueId < _other.valueId;
        }
    };
    
    /// Allows to pretty print indices, giving the valueId and span.
    std::ostream& operator<<(std::ostream& _out, const xerus::Index& _idx); 
}
