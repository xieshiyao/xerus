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

#include <cmath>
#include <type_traits>
#include <iomanip>
#include <mutex>
#include <time.h>

#include "stringUtilities.h"
#include "callStack.h"

// ------------------------------------------------------------------------------------------------------------------------------
// Hash the given name (so that we can use an int as template argument, strings are not possible!)
// using FNV-1a standard hash
constexpr uint64_t ___log_namehash(const char* x) {
    return *x ? (uint64_t(*x) ^ ___log_namehash(x+1))*1099511628211ul : 14695981039346656037ul;
}


// Define the loglevels
#define STRINGIFY2( x) #x
#define STRINGIFY(x) STRINGIFY2(x)
#define PASTE2( a, b) a##b
#define PASTE( a, b) PASTE2( a, b)

namespace err {
    enum {
        NOT_LOGGING = 0,
        LOGGING_ON_ERROR,
        LOGGING_FULL
    };
}

#define SET_LOGGING(name, value) \
    template<> struct ___logFlag<___log_namehash(STRINGIFY(name))>{ static const int flag = value; };
    
#define SET_LOGGING_DEFAULT(value) \
    template<uint64_t lvl> struct ___logFlag { static const int flag = value; }; \
    SET_DEFAULT_LOG_LEVELS
    

// default log levels  
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct ___logFlag<___log_namehash("fatal")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("critical")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("error")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("warning")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("info")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("debug")>{ static const int flag = err::LOGGING_ON_ERROR; }; 
    
#ifdef CRITICAL_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct ___logFlag<___log_namehash("fatal")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("critical")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("error")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("warning")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("info")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("debug")>{ static const int flag = err::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef ERROR_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct ___logFlag<___log_namehash("fatal")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("critical")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("error")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("warning")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("info")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("debug")>{ static const int flag = err::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef WARNING_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct ___logFlag<___log_namehash("fatal")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("critical")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("error")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("warning")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("info")>{ static const int flag = err::LOGGING_ON_ERROR; };\
    template<> struct ___logFlag<___log_namehash("debug")>{ static const int flag = err::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef INFO_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct ___logFlag<___log_namehash("fatal")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("critical")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("error")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("warning")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("info")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("debug")>{ static const int flag = err::LOGGING_ON_ERROR; }; 
#endif
    
#ifdef DEBUG_
#undef SET_DEFAULT_LOG_LEVELS
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct ___logFlag<___log_namehash("fatal")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("critical")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("error")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("warning")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("info")>{ static const int flag = err::LOGGING_FULL; };\
    template<> struct ___logFlag<___log_namehash("debug")>{ static const int flag = err::LOGGING_FULL; }; 
#endif


    
    
#ifdef LOGFILE_
    #include <fstream>
    #define ___LOGSTREAM err::fileStream
    namespace err {
        extern std::ofstream fileStream;
    }
#else
    #include <iostream>
    #define ___LOGSTREAM std::cerr
#endif  


#define ___LOGBUFFER (*err::bufferStreams.curr)
#include <sstream>
#include <limits>
namespace err {
    struct BufferStreams {
        std::stringstream *curr;
        std::stringstream *old;
        BufferStreams() {
            old = new std::stringstream;
            curr = new std::stringstream;
        }
        void checkSwitch() {
            if (curr->str().size() > 1024*1024) {
                delete old; 
                old = curr;
                curr = new std::stringstream;
            }
        }
        ~BufferStreams() {
            delete old;
            delete curr;
        }
    };
    extern BufferStreams bufferStreams;
    void dump_log_buffer(std::string _comment);
    extern std::mutex namedLoggerMutex;
    extern std::string logFilePrefix;
}
 
    
// no exit if inside unit-test
#ifdef TEST_
	#include "exceptions.h"
	#define exit(a) XERUS_THROW(::MISC::generic_error() << "exit(a) called (ignored by unit-testing)\ncallstack:\n" << ::MISC::get_call_stack());
#endif
    
#define ___LOG_TIME std::right << (1900+__ltm->tm_year) << '-' <<std::setw(2) << std::setfill('0') <<  __ltm->tm_mon \
				<< '-' <<std::setw(2) << std::setfill('0') <<  __ltm->tm_mday \
				<< ' ' << std::setw(2) << std::setfill('0') << __ltm->tm_hour \
				<< ':' <<std::setw(2) << std::setfill('0') <<  __ltm->tm_min \
				<< ':' <<std::setw(2) << std::setfill('0') <<  __ltm->tm_sec << " " \
				<< std::setfill(' ') << std::setw(20) << std::left << ::MISC::explode(__FILE__, '/').back() << ":" \
				<< std::right << std::setfill(' ') << std::setw(4) <<__LINE__ << " : " \
				<< std::setfill(' ') << std::setw(12) << std::left
    
#define COMPILE_TIME_EVAL(e) (std::integral_constant<decltype(e), e>::value)
    
#ifndef TEST_
#ifndef LOG_BUFFER_
    #define LOG(lvl, ...) \
        if (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag == err::LOGGING_FULL) { \
			time_t __t=std::time(0); tm *__ltm = localtime(&__t);\
			err::namedLoggerMutex.lock(); \
            ___LOGSTREAM << ___LOG_TIME << std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << std::endl; \
			err::namedLoggerMutex.unlock(); \
            if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("fatal"))) { exit(1); }; \
        } else \
            (void)0
            
#else
    #define LOG(lvl, ...) \
        if (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag != err::NOT_LOGGING) { \
			std::stringstream __tmp; \
			time_t __t=std::time(0); tm *__ltm = localtime(&__t);\
			__tmp << ___LOG_TIME << std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << std::endl; \
			err::namedLoggerMutex.lock(); \
            if (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag == err::LOGGING_FULL) { \
                ___LOGSTREAM << __tmp.str() << std::flush; \
            } \
            ___LOGBUFFER << __tmp.str(); \
            err::bufferStreams.checkSwitch(); \
            if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("error"))) {err::dump_log_buffer(std::string("error invoked:\n")+__tmp.str()); }; \
            if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("critical"))) {err::dump_log_buffer(std::string("critical error invoked:\n")+__tmp.str()); }; \
            if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("fatal"))) {err::dump_log_buffer(std::string("fatal error invoked:\n")+__tmp.str()); exit(1); }; \
			err::namedLoggerMutex.unlock(); \
        } else \
            (void)0
// in the above note the following:
// - everything is piped into a stringstream for two reasons: 1. so that functions appearing in msg will only be executed once
//                                                            2. so that buffer(pointers) can be piped as usual even though they empty the buffer while being piped
// - in performance critical programs logging should be disabled completely, the mutex in the above code will block the whole function. a read/write mutex could
//   increase performance a bit (by only making sure that the delete inside checkSwitch is not called why another thread is piping something) and the mutex
//   in the LOG_NO_BUFFER version could be disabled completely (the streams are thread"safe" (only the output will be jumbled a bit)
            
#endif

#else 
// ie TEST_ is defined
#ifndef LOG_BUFFER_
    #define LOG(lvl, ...) \
        if (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag == err::LOGGING_FULL) { \
			time_t __t=std::time(0); tm *__ltm = localtime(&__t);\
			if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("fatal"))) { \
				std::stringstream ss; \
				ss << std::endl << "\033[36m" << ___LOG_TIME << std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << "\033[0m";\
				XERUS_THROW(::MISC::generic_error() << ss.str() << "\ncallstack:\n" << ::MISC::get_call_stack()); \
			} else { \
				err::namedLoggerMutex.lock(); \
				___LOGSTREAM << ___LOG_TIME << std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << std::endl; \
				err::namedLoggerMutex.unlock(); \
			} \
        } else \
            (void)0
            
#else
    #define LOG(lvl, ...) \
        if (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag != err::NOT_LOGGING) { \
			std::stringstream __tmp; \
			time_t __t=std::time(0); tm *__ltm = localtime(&__t); \
			__tmp << ___LOG_TIME << std::string(STRINGIFY(lvl) ": ") << __VA_ARGS__ << std::endl; \
			if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("fatal"))) { \
				err::namedLoggerMutex.lock(); \
				___LOGBUFFER << __tmp.str(); \
				err::bufferStreams.checkSwitch(); \
				err::dump_log_buffer(std::string("fatal error invoked (might be expected):\n")+__tmp.str()); \
				err::namedLoggerMutex.unlock(); \
				XERUS_THROW(::MISC::generic_error() << __tmp.str() << "\ncallstack:\n" << MISC::get_call_stack()); \
			} else { \
				err::namedLoggerMutex.lock(); \
				if (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag == err::LOGGING_FULL) { \
					___LOGSTREAM << __tmp.str() << std::flush; \
				} \
				___LOGBUFFER << __tmp.str(); \
				err::bufferStreams.checkSwitch(); \
				if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("error"))) {err::dump_log_buffer(std::string("error invoked:\n")+__tmp.str()); }; \
				if (COMPILE_TIME_EVAL(___log_namehash(STRINGIFY(lvl))==___log_namehash("critical"))) {err::dump_log_buffer(std::string("critical error invoked:\n")+__tmp.str()); }; \
				err::namedLoggerMutex.unlock(); \
			} \
        } else \
            (void)0


#endif
#endif

#define IS_LOGGING(lvl) \
    (___logFlag<___log_namehash(STRINGIFY(lvl))>::flag != err::NOT_LOGGING)

    
SET_LOGGING_DEFAULT(err::LOGGING_FULL)

/*---------------------- Get name of current object -----------------*/
#define OBJECT_NAME GET_OBJECT_NAME(this)

template<typename T>
std::string GET_OBJECT_NAME(const T* const) {
	std::string name = typeid(T).name();
	return ::MISC::demangle_cxa(name);
}

