// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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
 * @brief Header file for unit test and code-coverage helper classes.
 */

#pragma once

/**
 * @def UNIT_TEST(grp, name, code)
 * @brief Defines a unit test of name @a name in the group @a grp dat performs the source code @a code.
 * @details Will automatically be performed if it is part of the executable and the application was built with -D TEST_.
 */

/**
 * @def REQUIRE_TEST
 * @brief Marked position for the code-coverage test. Any not passed REQUIRE_TEST macro will result in a warning.
 * REQUIRE_TEST is implied by any REQUIRE(...) macro if test coverage is enabled.
 */

#ifdef TEST_
	#include "namedLogger.h"
	
	#include <map>
    #include <string>
    #include <functional>

    #define PRINTCHECK std::cout << u8"\033[1;32m\u2713 \033[0m" << std::flush
    #define PRINTFAIL  std::cout << u8"\033[1;31m\u2717 \033[0m" << std::flush

	#define PASTE2( a, b) a##b
	#define PASTE( a, b) PASTE2( a, b)

	#define TEST(...) if (!(__VA_ARGS__)) {PRINTFAIL; LOG(error, #__VA_ARGS__ << " failed"); passed = false;} else {PRINTCHECK;} void(0)
	#define MTEST(cond, ...) if (!(cond)) {PRINTFAIL; LOG(error, #cond << " failed, msg: " << __VA_ARGS__); passed = false;} else {PRINTCHECK;} void(0)
	
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
			{\
				bool failtestFailed = false; \
				xerus::misc::internal::logFilePrefix = "failtest/"; xerus::misc::internal::silenced = true; \
				try { test; } catch (...) {failtestFailed = true; PRINTCHECK;} \
				xerus::misc::internal::logFilePrefix.clear();  xerus::misc::internal::silenced = false;  \
				if(!failtestFailed) { PRINTFAIL; LOG(error, #test << " returned without error"); passed = false; } \
			}\
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
			
			
	#define UNIT_TEST2(grp, name) \
		xerus::misc::internal::UnitTest PASTE(grp,name)(#grp, #name, [](bool& passed) {\
				std::mt19937_64 rnd;\
				rnd.seed(0xC0CAC01A);\
				\
				std::normal_distribution<value_t> normalDist (0.0, 1.0);\
				std::uniform_real_distribution<value_t> uniformDist (-1.0, 1.0);\
			
	
	#define main(...) original_main_function_that_was_disabled_by_xerus_unit_test_enviroment_horst( __VA_ARGS__ )

    namespace xerus { namespace misc { namespace internal {
		struct UnitTest final {
			// order of contruction of global objects is random so the first one has to create the map
			static std::map<std::string, std::map<std::string, std::function<bool ()>>> *tests;
			
			UnitTest(std::string _group, std::string _name, std::function<bool ()> _f);
			UnitTest(std::string _group, std::string _name, std::function<void(bool&)> _f);
		};

		#ifdef TEST_COVERAGE_
			struct RequiredTest {
				struct Identifier {
					std::string functionName;
					std::string filename;
					size_t lineNumber;
					Identifier(std::string _func, std::string _file, size_t _line);
					bool operator<(const Identifier &_rhs) const;
				};
				
				// order of construction of global objects is random so the first one has to create the map
				static std::map<Identifier, size_t> *tests;
				
				static void register_test(std::string _functionName, std::string _fileName, size_t _lineNb);
				
				static void increase_counter(std::string _functionName, std::string _fileName, size_t _lineNb);
			};
		#endif
	}}}
#else
	#define UNIT_TEST(grp, name, ...)
	#define UNIT_TEST2(grp, name, ...)
	#define REQUIRE_TEST (void)0
#endif
