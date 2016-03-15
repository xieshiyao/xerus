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
 * @brief Header file for all logging macros and log-buffer functionality.
 */

#pragma once

#include <mutex>
#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "callStack.h"

#ifdef LOGFILE_
    #include <fstream>
    #define XERUS_LOGSTREAM xerus::misc::internal::fileStream
#else
    #include <iostream>
    #define XERUS_LOGSTREAM std::cerr
#endif

namespace xerus {
    namespace misc {
        namespace internal {
			/**
			 * @brief Hashes a given c-string using the FNV-1a standard hash.
			 * @details This is used eg. to use strings as template arguments.
			 */
            constexpr uint64_t log_namehash(const char* x) {
                return *x ? (uint64_t(*x) ^ xerus::misc::internal::log_namehash(x+1))*1099511628211ul : 14695981039346656037ul;
            }

			// If the LOG_BUFFER is aktive there is the additional option only to print the log if an error occours.
			#ifdef LOG_BUFFER_
				enum {
					NOT_LOGGING = 0,
					LOGGING_ON_ERROR = 1,
					LOGGING_FULL = 2
				};
			#else
				enum {
					NOT_LOGGING = 0,
					LOGGING_FULL = 2
				};
				static const auto LOGGING_ON_ERROR = NOT_LOGGING;
			#endif

            extern std::mutex namedLoggerMutex;
            extern std::string logFilePrefix;
            extern bool silenced;
			extern std::chrono::system_clock::time_point programStartTime;

			extern std::ofstream fileStream;

			namespace buffer {
				extern std::stringstream current;
				extern std::stringstream old;
				
				void checkSwitch();
				
				void dump_log(std::string _comment);
			}
        }
    }
}

#define STRINGIFY2( x) #x
#define STRINGIFY(x) STRINGIFY2(x)
#define PASTE2( a, b) a##b
#define PASTE( a, b) PASTE2( a, b)

/**
 * @def SET_LOGGING(lvl, value)
 * @brief set the logging behaviour of severity level @a lvl to @a value (either NOT_LOGGING, LOGGING_ON_ERROR or LOGGING_FULL)
 * @details this definition must not be repeated and must be define din a global header that is included before any msg is logged with that lvl
 */
#define SET_LOGGING(lvl, value) \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>{ static const int flag = value; };

/**
 * @def SET_LOGGING_DEFAULT(value)
 * @brief sets the logging behaviour of all levels that are not otherwise specified to @a value
 */
#define SET_LOGGING_DEFAULT(value) \
    template<uint64_t lvl> struct XERUS_logFlag { static const int flag = value; }; \
    SET_DEFAULT_LOG_LEVELS
    

// Default log levels  
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; }; 
    
#ifdef FATAL_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; }; 
#endif

#ifdef CRITICAL_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef ERROR_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef WARNING_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef INFO_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef DEBUG_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("debug")>{ static const int flag = xerus::misc::internal::LOGGING_FULL; }; 
#endif


// No exit if inside unit-test 
#if defined(TEST_) && defined(NO_XERUS_EXCEPTIONS)
    #pragma warning "Tried to compile with TEST but without exceptions... failtests will exit the program."
#endif

    
#ifndef NO_XERUS_EXCEPTIONS
	#include "exceptions.h"
#endif

/**
 * @def COMPILE_TIME_EVAL(e)
 * @brief forces the compiler to evaluate @a e during compilation
 */
#define COMPILE_TIME_EVAL(e) (std::integral_constant<decltype(e), e>::value)


// in the following note that:
// - everything is piped into a stringstream for two reasons: 1. so that functions appearing in msg will only be executed once
//                                                            2. so that buffer(pointers) can be piped as usual even though they empty the buffer while being piped
// - in performance critical programs logging should be disabled completely, the mutex in the following code will block the whole function. a read/write mutex could
//   increase performance a bit (by only making sure that the delete inside checkSwitch is not called while another thread is piping something) and the mutex
//   in the LOG_NO_BUFFER version could be disabled completely (the streams are thread"safe" (only the output will be jumbled a bit)

#ifdef NO_XERUS_EXCEPTIONS
    #define NAMED_LOGGER_ON_FATAL \
        if (!xerus::misc::internal::silenced) { \
            if (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag != xerus::misc::internal::LOGGING_FULL) { \
                XERUS_LOGSTREAM << tmpStream.str() << "callstack:\n" << xerus::misc::get_call_stack());\
            } else { \
                XERUS_LOGSTREAM << "callstack:\n" << xerus::misc::get_call_stack());\
            } \
        } \
        exit(1);
#else // with xerus exceptions
    #define NAMED_LOGGER_ON_FATAL XERUS_THROW(xerus::misc::generic_error() << tmpStream.str() << "callstack:\n" << xerus::misc::get_call_stack());
#endif


#ifdef LOG_BUFFER_
    #define NAMED_LOGGER_LOGBUFFER \
        xerus::misc::internal::buffer::current << tmpStream.str(); \
        xerus::misc::internal::buffer::checkSwitch(); \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("error"))) {xerus::misc::internal::buffer::dump_log(std::string("error invoked:\n")+tmpStream.str()); }; \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("critical"))) {xerus::misc::internal::buffer::dump_log(std::string("critical error invoked:\n")+tmpStream.str()); }; \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("fatal"))) {xerus::misc::internal::buffer::dump_log(std::string("fatal error invoked:\n")+tmpStream.str()); }; 
#else // no log buffer
    #define NAMED_LOGGER_LOGBUFFER
#endif

#ifdef LOG_ABSOLUTE_TIME
	#define XERUS_LOGGER_TIMESTAMP \
		std::time_t xerus_err_t=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); \
		tmpStream \
				<< std::put_time(std::localtime(&xerus_err_t), "%F %T") << " "
#else
	#define XERUS_LOGGER_TIMESTAMP \
		std::chrono::system_clock::time_point xerus_err_t = std::chrono::system_clock::now();\
		auto xerus_err_timediff = std::chrono::duration_cast<std::chrono::milliseconds>(xerus_err_t - xerus::misc::internal::programStartTime).count(); \
		tmpStream \
				<< std::right << '+' << std::setw(2) << std::setfill('0') << (xerus_err_timediff/3600000) \
				<< ':' << std::setw(2) << std::setfill('0') << ((xerus_err_timediff/60000)%60) \
				<< ':' << std::setw(2) << std::setfill('0') <<  ((xerus_err_timediff/1000)%60) \
				<< ',' << std::setw(3) << std::setfill('0') <<  (xerus_err_timediff%1000) << " "
#endif


/**
 * @def LOG(lvl, msg)
 * @brief logs the message @a msg with severity level @a lvl
 * @details the exact behaviour is modified by the SET_DEFAULT_LOG_LEVELS and SET_LOGGING macros. In case @a lvl is not being logged with the
 *   current configuration, this macro evaluates to an `if (false) {}` expression and is fully removed by the compiler.
 */
#define LOG(lvl, ...) \
    if (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag != xerus::misc::internal::NOT_LOGGING) { \
        std::stringstream tmpStream; \
        XERUS_LOGGER_TIMESTAMP \
				<< std::setfill(' ') << std::setw(20) << std::left << xerus::misc::explode(__FILE__, '/').back() << ":" \
				<< std::right << std::setfill(' ') << std::setw(4) <<__LINE__ << " : " \
				<< std::setfill(' ') << std::setw(12) << std::left \
				<< std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << std::endl; \
        xerus::misc::internal::namedLoggerMutex.lock(); \
        if (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag == xerus::misc::internal::LOGGING_FULL && !xerus::misc::internal::silenced) { \
            XERUS_LOGSTREAM << tmpStream.str() << std::flush; \
        } \
        NAMED_LOGGER_LOGBUFFER \
        xerus::misc::internal::namedLoggerMutex.unlock(); \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("fatal"))) { \
            NAMED_LOGGER_ON_FATAL \
        } \
    } else \
        (void)0

/**
 * @def LOG_SHORT(lvl, msg)
 * @brief logs the message @a msg with severity level @a lvl, omits the current file name and line number
 * @details the exact behaviour is modified by the SET_DEFAULT_LOG_LEVELS and SET_LOGGING macros. In case @a lvl is not being logged with the
 *   current configuration, this macro evaluates to an `if (false) {}` expression and is fully removed by the compiler.
 */
#define LOG_SHORT(lvl, ...) \
    if (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag != xerus::misc::internal::NOT_LOGGING) { \
        std::stringstream tmpStream; \
        XERUS_LOGGER_TIMESTAMP \
				<< std::setfill(' ') << std::setw(12) << std::left \
				<< std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << std::endl; \
        xerus::misc::internal::namedLoggerMutex.lock(); \
        if (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag == xerus::misc::internal::LOGGING_FULL && !xerus::misc::internal::silenced) { \
            XERUS_LOGSTREAM << tmpStream.str() << std::flush; \
        } \
        NAMED_LOGGER_LOGBUFFER \
        xerus::misc::internal::namedLoggerMutex.unlock(); \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("fatal"))) { \
            NAMED_LOGGER_ON_FATAL \
        } \
    } else \
        (void)0

/**
 * @def LOG_ONCE(lvl, msg)
 * @brief logs the message @a msg with severity level @a lvl at most once per program execution
 * @details the exact behaviour is modified by the SET_DEFAULT_LOG_LEVELS and SET_LOGGING macros. In case @a lvl is not being logged with the
 *   current configuration, this macro evaluates to an `if (false) {}` expression and is fully removed by the compiler.
 */
#define LOG_ONCE(lvl, ...) \
	{\
		static bool logged = false;\
		if (!logged) {\
			LOG(lvl, __VA_ARGS__);\
			logged = true;\
		}\
	}

/**
 * @def IS_LOGGING(lvl)
 * @brief evaluates to true if @a lvl is begin logged (either to cerr or into a file on error) in the current configuration 
 */
#define IS_LOGGING(lvl) \
    (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag != xerus::misc::internal::NOT_LOGGING)


SET_LOGGING_DEFAULT(xerus::misc::internal::LOGGING_FULL)
