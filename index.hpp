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
#include "index.h"
#include "tensor.h"

namespace xerus {
    std::atomic<size_t> Index::idThreadInitCounter(0);
    thread_local size_t Index::idCounter = (idThreadInitCounter++)<<54;
    
    Index::Index() : valueId(idCounter++), span(1) { REQUIRE(idCounter < 1ul<<54, "Index ID counter left thread safe range."); }
    
    Index::Index(const long _i) : valueId(_i), span(1) {
        REQUIRE(_i >= 0, "Negative valueId= " <<_i<< " given");
        flags[Flag::FIXED] = true;
    }
    
    
    Index::Index(const size_t _valueId, const size_t _span) : valueId(_valueId), span(_span) {}
    
    Index::Index(const size_t _valueId, const size_t _span, const Flag _flag1, const bool _flagValue1) : valueId(_valueId), span(_span) {
        IF_CHECK(
        if(_flag1 == Flag::OPEN) {
            flags[Flag::OPEN_CHECKED] = true;
        })
        
        flags[_flag1] = _flagValue1;
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const Flag _flag1, const Flag _flag2, const bool _flagValue1, const bool _flagValue2) : valueId(_valueId), span(_span) {
        IF_CHECK(
        if(_flag1 == Flag::OPEN || _flag2 == Flag::OPEN) {
            flags[Flag::OPEN_CHECKED] = true;
        })
        
        flags[_flag1] = _flagValue1;
        flags[_flag2] = _flagValue2;
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const Flag _flag1, const Flag _flag2, const Flag _flag3, const bool _flagValue1, const bool _flagValue2, const bool _flagValue3) : valueId(_valueId), span(_span) {
        IF_CHECK(
        if(_flag1 == Flag::OPEN || _flag2 == Flag::OPEN || _flag3 == Flag::OPEN) {
            flags[Flag::OPEN_CHECKED] = true;
        })
        
        flags[_flag1] = _flagValue1;
        flags[_flag2] = _flagValue2;
        flags[_flag3] = _flagValue3;
    }
    
    std::ostream& operator<<(std::ostream& _out, const xerus::Index& _idx) {
        _out << _idx.valueId << "(" << (_idx.flags[Index::Flag::INVERSE_SPAN] ? "-" : "") << _idx.span << ")";
        return _out;
    }
}

