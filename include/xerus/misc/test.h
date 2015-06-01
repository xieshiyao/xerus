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

#include "namedLogger.h"
#include "testManager.h"

// Only define the check macros if DISABLE_RUNTIME_CHECKS_ is inaktive
#ifndef DISABLE_RUNTIME_CHECKS_

	#define PCHECK(precondition, condition, level, message) if(IS_LOGGING(level)) { precondition; if(condition) { LOG(level, #condition << " failed msg: " << message); }} else void(0)

	#ifdef TEST_
		#define CHECK(condition, level, message) REQUIRE_TEST; if(IS_LOGGING(level) && !(condition)) { LOG(level, #condition << " failed msg: " << message); } else void(0)
	#else
		#define CHECK(condition, level, message) if(IS_LOGGING(level) && !(condition)) { LOG(level, #condition << " failed msg: " << message); } else void(0)
	#endif

	#define REQUIRE(condition, message) CHECK(condition, fatal, message)
	
    #define IF_CHECK(expression) expression

#else

	#define PCHECK(precondition, condition, level, message) void(0)

	#define CHECK(condition, level, message) void(0)

	#define REQUIRE(condition, message) void(0)
	
    #define IF_CHECK(expression)  

#endif
 
#ifdef TEST_
    #define PRINTCHECK std::cout << u8"\033[1;32m\u2713 \033[0m" << std::flush
    #define PRINTFAIL  std::cout << u8"\033[1;31m\u2717 \033[0m" << std::flush

	#define PASTE2( a, b) a##b
	#define PASTE( a, b) PASTE2( a, b)

	#define TEST(...) if (!(__VA_ARGS__)) {PRINTFAIL; LOG(error, #__VA_ARGS__ << " failed"); passed = false;} else {PRINTCHECK;} void(0)
	
	#ifdef TEST_COVERAGE_
		#define REQUIRE_TEST \
			do { \
				static const char * xerus_test_fname = __PRETTY_FUNCTION__;\
				struct xerus_test_a{ static void rt() {\
					xerus::misc::internal::RequiredTest::register_test(xerus_test_fname, __FILE__, __LINE__);\
				} };\
				typedef void (*xerus_test_t)();\
				static xerus_test_t xerus_test_rtp __attribute__((section(".init_array"))) = &xerus_test_a::rt; \
				(void)xerus_test_rtp; \
				xerus::misc::internal::RequiredTest::increase_counter(xerus_test_fname, __FILE__, __LINE__); \
			} while(false)
	#else
		#define REQUIRE_TEST (void)0
	#endif
	
	#ifndef DISABLE_RUNTIME_CHECKS_
		#define FAILTEST(test) \
            xerus::misc::internal::logFilePrefix = "failtest/"; xerus::misc::internal::silenced = true; \
            try { test; PRINTFAIL; LOG(error, #test << " returned without error"); passed = false; } catch (...) {PRINTCHECK;} \
            xerus::misc::internal::logFilePrefix.clear();  xerus::misc::internal::silenced = false;  \
            void(0)
	#else
		#define FAILTEST(test) LOG(warning, "Failtest is not useful with flag DISABLE_RUNTIME_CHECKS_")
	#endif
	
	#define UNIT_TEST(grp, name, ...) \
		xerus::misc::internal::UnitTest PASTE(grp,name)(#grp, #name, []()->bool{\
				bool passed = true;\
				__VA_ARGS__\
				return passed;\
			});
			
	#define main(...) original_main_function_that_was_disabled_by_xerus_unit_test_enviroment_horst( __VA_ARGS__ )

#else
	#define REQUIRE_TEST (void)0
	#define UNIT_TEST(grp, name, ...)

#endif

