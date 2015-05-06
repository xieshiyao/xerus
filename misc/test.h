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

// Only define the check macros if CHECK_ is aktive
#ifdef CHECK_

	#define PCHECK(precondition, condition, level, message) if(IS_LOGGING(level)) { precondition; if(condition) { LOG(level, #condition << " failed msg: " << message); }} else void(0)

	#ifdef TEST_
		#define CHECK(condition, level, message) REQUIRE_TEST; if(IS_LOGGING(level) && !(condition)) { LOG(level, #condition << " failed msg: " << message); } else void(0)
	#else
		#define CHECK(condition, level, message) if(IS_LOGGING(level) && !(condition)) { LOG(level, #condition << " failed msg: " << message); } else void(0)
	#endif

	#define REQUIRE(condition, message) CHECK(condition, fatal, message)

#else

	#define PCHECK(precondition, condition, level, message) void(0)

	#define CHECK(condition, level, message) void(0)

	#define REQUIRE(condition, message) void(0)

#endif
 
#ifdef TEST_
    #define PRINTCHECK std::cout << u8"\033[1;32m\u2713 \033[0m" << std::flush
    #define PRINTFAIL  std::cout << u8"\033[1;31m\u2717 \033[0m" << std::flush

	#define PASTE2( a, b) a##b
	#define PASTE( a, b) PASTE2( a, b)

	#define TEST(cond) if (!(cond)) {PRINTFAIL; LOG(error, #cond << " failed"); passed = false;} else {PRINTCHECK;} void(0)
	
	#define REQUIRE_TEST \
		do { \
			static const char * ___fname = __PRETTY_FUNCTION__;\
			struct ___a{ static void rt() {\
				___RequiredTest::register_test(___fname, __FILE__, __LINE__);\
			} };\
			typedef void (*___t)();\
			static ___t ___rtp __attribute__((section(".init_array"))) = &___a::rt; \
			(void)___rtp; \
			___RequiredTest::increase_counter(___fname, __FILE__, __LINE__); \
		} while(false)
	
	#ifdef CHECK_
		#define FAILTEST(test) \
			err::logFilePrefix = "failtest/"; err::silenced = true; \
			try { test; PRINTFAIL; LOG(error, #test << " returned without error"); passed = false; } catch (...) {PRINTCHECK;} \
			err::logFilePrefix.clear();  err::silenced = false;  \
			void(0)
	#else
		#define FAILTEST(test) LOG(warning, "Failtest is not useful without flag CHECK_")
	#endif
	
	#define UNIT_TEST(grp, name, ...) \
		___UnitTest PASTE(grp,name)(#grp, #name, []()->bool{\
				bool passed = true;\
				__VA_ARGS__\
				return passed;\
			});
			
	#define main(...) ___horst_main_will_not_be_called( __VA_ARGS__ )

#else

	#define UNIT_TEST(grp, name, ...)

#endif

