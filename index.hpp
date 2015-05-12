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
    
    bool Index::all_open(const std::vector<Index>& _indices) {
        for(const Index& idx : _indices) {
            if(!idx.open()) { return false; }
        }
        return true;
    }
            
    std::atomic<size_t> Index::idThreadInitCounter(0);
    thread_local size_t Index::idCounter = (idThreadInitCounter++)<<54;
    
    Index::Index() : valueId(idCounter++), span(1) { REQUIRE(idCounter < 1ul<<54, "Index ID counter left thread safe range."); }
    
    Index::Index(const long _i) : valueId(_i), span(1) {
        REQUIRE(_i >= 0, "Negative valueId= " <<_i<< " given");
        flags[Flag::FIXED] = true;
    }
    
    
    Index::Index(const size_t _valueId, const size_t _span) : valueId(_valueId), span(_span) {}
    
    Index::Index(const size_t _valueId, const size_t _span, const Flag _flag1, const bool _flagValue1) : valueId(_valueId), span(_span) {
        flags[_flag1] = _flagValue1;
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const size_t _dimension) : valueId(_valueId), span(_span), assingedDimension(_dimension) {}
    
    Index::Index(const size_t _valueId, const size_t _span, const size_t _dimension, const Flag _flag1, const bool _flagValue1) : valueId(_valueId), span(_span), assingedDimension(_dimension) {
        IF_CHECK( flags[Flag::ASSINGED] = true; )
        flags[_flag1] = _flagValue1;
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const size_t _dimension, const Flag _flag1, const Flag _flag2, const bool _flagValue1, const bool _flagValue2) : valueId(_valueId), span(_span), assingedDimension(_dimension) {
        IF_CHECK( flags[Flag::ASSINGED] = true; )
        flags[_flag1] = _flagValue1;
        flags[_flag2] = _flagValue2;
    }
    
    std::ostream& operator<<(std::ostream& _out, const xerus::Index& _idx) {
        _out << _idx.valueId << "(" << (_idx.flags[Index::Flag::INVERSE_SPAN] ? "-" : "") << _idx.span << ")";
        return _out;
    }
}

