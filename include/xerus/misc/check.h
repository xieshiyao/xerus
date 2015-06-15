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
 * @brief Header file for CHECK and REQUIRE macros.
 */

#pragma once

/**
 * @def CHECK(condition, level, message)
 * @brief Checks whether @a condition is true and calls LOG(level, message) otherwise.
 * @details The check is omitted if level is not logged anyways.
 */

/**
 * @def REQUIRE(condition, message)
 * @brief Checks whether @a condition is true. logs a fatal error otherwise via LOG(fatal, message).
 */

/**
 * @def IF_CHECK(expression)
 * @brief Executes @a expression only if the compilation was without DISABLE_RUNTIME_CHECKS_ set.
 */

#ifdef TEST_COVERAGE_
	#include "test.h"
#else
	#define REQUIRE_TEST (void)0
#endif

// Only do anything if the check macros if DISABLE_RUNTIME_CHECKS_ is inaktive
#ifndef DISABLE_RUNTIME_CHECKS_
	#include "namedLogger.h"

	#ifdef TEST_COVERAGE_
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
