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

// Only do anything if the check macros if DISABLE_RUNTIME_CHECKS_ is inaktive
#ifndef DISABLE_RUNTIME_CHECKS_
	#include "namedLogger.h"

	#ifdef TEST_COVERAGE_
		#include "test.h"
		
		#define CHECK(condition, level, message) REQUIRE_TEST; if(IS_LOGGING(level) && !(condition)) { LOG(level, #condition << " failed msg: " << message); } else void(0)
	#else
		#define CHECK(condition, level, message) if(IS_LOGGING(level) && !(condition)) { LOG(level, #condition << " failed msg: " << message); } else void(0)
	#endif

	#define REQUIRE(condition, message) CHECK(condition, fatal, message)
	
    #define IF_CHECK(expression) expression
#else
	#define CHECK(condition, level, message) void(0)
	#define REQUIRE(condition, message) void(0)
    #define IF_CHECK(expression)  
#endif
