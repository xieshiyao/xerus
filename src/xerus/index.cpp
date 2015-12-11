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

/**
 * @file
 * @brief Implementation of the Index class.
 */

#include <xerus/index.h>

#include <xerus/misc/standard.h>
#include <xerus/misc/check.h>

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
    
    Index::Index(const size_t _i) : valueId(_i), span(1) {
        flags[Flag::FIXED] = true;
    }
    
    Index::Index(const int _i) : Index((size_t) _i) {
        REQUIRE(_i >= 0, "Negative valueId= " <<_i<< " given");
    }
    
    Index::Index(const size_t _valueId, const size_t _span) : valueId(_valueId), span(_span) {}
    
    Index::Index(const size_t _valueId, const size_t _span, const Flag _flag1, const bool _flagValue1) : valueId(_valueId), span(_span) {
        flags[_flag1] = _flagValue1;
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const size_t _dimension) : valueId(_valueId), span(_span), assingedDimension(_dimension) {
        IF_CHECK( flags[Flag::ASSINGED] = true; )
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const size_t _dimension, const Flag _flag1, const bool _flagValue1) : valueId(_valueId), span(_span), assingedDimension(_dimension) {
        IF_CHECK( flags[Flag::ASSINGED] = true; )
        flags[_flag1] = _flagValue1;
    }
    
    Index::Index(const size_t _valueId, const size_t _span, const size_t _dimension, const Flag _flag1, const Flag _flag2, const bool _flagValue1, const bool _flagValue2) : valueId(_valueId), span(_span), assingedDimension(_dimension) {
        IF_CHECK( flags[Flag::ASSINGED] = true; )
        flags[_flag1] = _flagValue1;
        flags[_flag2] = _flagValue2;
    }
	
	void Index::set_span(const size_t _degree) {
		REQUIRE(!flags[Flag::FIXED] || span == 1, "Fixed indices must have span one.");
		if(flags[Flag::INVERSE_SPAN]) {
			REQUIRE(!flags[Flag::FIXED], "Fixed indices must not have inverse span."); 
			REQUIRE(span <= _degree, "Index with inverse span would have negative actual span. Tensor degree: " << _degree << ", inverse span " << span);
			span = _degree-span;
			flags.reset();
		} else if( flags[Flag::FRACTIONAL_SPAN] ) {
			REQUIRE(!flags[Flag::FIXED], "Fixed indices must not have fractional span.");
			REQUIRE(_degree%span == 0, "Fractional span must divide the tensor degree. Here tensor degree = " << _degree << ", span = " << span);
			span = _degree/span;
			flags.reset();
		}
	}
	
	size_t Index::actual_span(const size_t _degree) const {
		if(flags[Flag::INVERSE_SPAN]) {
			REQUIRE(!flags[Flag::FIXED], "Fixed indices must not have inverse span."); 
			REQUIRE(span <= _degree, "Index with inverse span would have negative actual span. Tensor degree: " << _degree << ", inverse span " << span);
			return _degree-span;
		} else if( flags[Flag::FRACTIONAL_SPAN] ) {
			REQUIRE(!flags[Flag::FIXED], "Fixed indices must not have fractional span.");
			REQUIRE(_degree%span == 0, "Fractional span must divide the tensor degree. Here tensor degree = " << _degree << ", span = " << span);
			return _degree/span;
		} else {
			REQUIRE(!flags[Flag::FIXED] || span == 1, "Fixed indices must have span one.");
			return span;
		}
	}
	
	bool Index::fixed() const {
		#ifndef DISABLE_RUNTIME_CHECKS_
			if(flags[Index::Flag::FIXED]) {
				REQUIRE(!flags[Flag::INVERSE_SPAN], "Fixed indices must not have inverse span."); 
				REQUIRE(!flags[Flag::FRACTIONAL_SPAN], "Fixed indices must not have fractional span.");
				REQUIRE(span == 1, "Fixed indices must have span one.");
				return true;
			} else {
				return false;
			}
		#else
			return flags[Index::Flag::FIXED];
		#endif
	}
    
    bool Index::open() const {
        REQUIRE(flags[Index::Flag::ASSINGED], "Check for index openness only allowed if the index is assinged.");
        return flags[Index::Flag::OPEN];
    }
    
    void Index::open(const bool _open) {
        flags[Flag::OPEN] = _open;
        IF_CHECK( flags[Flag::ASSINGED] = true; )
    }
    
    size_t Index::dimension() const {
        REQUIRE(flags[Index::Flag::ASSINGED], "Check for index dimension only allowed if the index is assinged.");
        return assingedDimension;
    }
    
    Index Index::operator^(const size_t _span) const {
        REQUIRE(flags.none(), "Cannot apply ^ operator to an index that has any flag set.");
        return Index(valueId, _span);
    }
    
    Index Index::operator&(const size_t _span) const {
        REQUIRE(flags.none(), "Cannot apply & operator to an index that has any flag set.");
        return Index(valueId, _span, Flag::INVERSE_SPAN);
    }
    
    Index Index::operator/(const size_t _span) const {
        REQUIRE(flags.none(), "Cannot apply & operator to an index that has any flag set.");
        return Index(valueId, _span, Flag::FRACTIONAL_SPAN);
    }
    
    bool Index::operator==(const Index& _other) const {
        return valueId == _other.valueId && !fixed() && !_other.fixed();
    }
    
    bool Index::operator!=(const Index& _other) const {
        return valueId != _other.valueId || fixed() || _other.fixed();
    }
    
    bool Index::operator<(const Index& _other) const {
        return valueId < _other.valueId;
    }
    
    
    std::ostream& operator<<(std::ostream& _out, const xerus::Index& _idx) {
        _out << _idx.valueId << "(" << (_idx.flags[Index::Flag::INVERSE_SPAN] ? "-" : "") << _idx.span << ")";
        return _out;
    }
}

