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
    std::atomic<long> Index::idThreadInitCounter(0);
    thread_local long Index::idCounter = -((idThreadInitCounter++)<<55);
    
    Index::Index() : valueId(--idCounter), span(1), inverseSpan(false) { REQUIRE(idCounter < 0, "Index ID counter overflowed"); }
    
    Index::Index(const long _i) : valueId(_i), span(1), inverseSpan(false) {
        CHECK(_i >= 0, fatal, "Negative valueId= " <<_i<< " given");
    }
    
    Index::Index(const long _valueId, const size_t _span, const bool _inverseSpan) : valueId(_valueId), span(_span), inverseSpan(_inverseSpan) { }
    
    std::ostream& operator<<(std::ostream& _out, const xerus::Index& _idx) {
        _out << _idx.valueId << "(" << (_idx.inverseSpan? "-" : "") << _idx.span << ")";
        return _out;
    }
}



