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
#include <sstream>

#include "stringUtilities.h"
#include "callStack.h"

#ifdef LOGFILE_
    #include <fstream>
    #define XERUS_LOGSTREAM xerus::misc::internal::fileStream
#else
    #include <iostream>
    #define XERUS_LOGSTREAM std::cerr
#endif

#ifdef LOG_BUFFER_
	#include <fstream>
#endif

namespace xerus {
    namespace misc {
        namespace internal {
            // Hash the given name (so that we can use an int as template argument, strings are not possible!)
            // using FNV-1a standard hash
            constexpr uint64_t log_namehash(const char* x) {
                return *x ? (uint64_t(*x) ^ xerus::misc::internal::log_namehash(x+1))*1099511628211ul : 14695981039346656037ul;
            }
            
#ifdef LOG_BUFFER_
            enum {
                NOT_LOGGING = 0,
                LOGGING_ON_ERROR,
                LOGGING_FULL
            };
#else
            enum {
                NOT_LOGGING = 0,
                LOGGING_FULL
            };
            static const auto LOGGING_ON_ERROR = NOT_LOGGING;
#endif
            
            extern std::mutex namedLoggerMutex;
            extern std::string logFilePrefix;
            extern bool silenced;
            
            #ifdef LOGFILE_
                extern std::ofstream fileStream;
            #endif
                
            #ifdef LOG_BUFFER_
                // TODO This class should be replaced by two global stringstream objects and one global function (both in this namespace of course).
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
            #endif
        }
    }
}

#define STRINGIFY2( x) #x
#define STRINGIFY(x) STRINGIFY2(x)
#define PASTE2( a, b) a##b
#define PASTE( a, b) PASTE2( a, b)

#define SET_LOGGING(name, value) \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(name))>{ static const int flag = value; };
    
#define SET_LOGGING_DEFAULT(value) \
    template<uint64_t lvl> struct XERUS_logFlag { static const int flag = value; }; \
    SET_DEFAULT_LOG_LEVELS
    

// Default log levels  
#define SET_DEFAULT_LOG_LEVELS \
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("fatal")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("critical")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("error")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("warning")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
    template<> struct XERUS_logFlag<xerus::misc::internal::log_namehash("info")>{ static const int flag = xerus::misc::internal::LOGGING_ON_ERROR; };\
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
    #pragma warning "Tried to compile with TEST but without exceptions... failtests would exit the program." // TODO would?
    #undef NO_XERUS_EXCEPTIONS
#endif

    
#ifndef NO_XERUS_EXCEPTIONS
	#include "exceptions.h"
#endif

    
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
        (*xerus::misc::internal::bufferStreams.curr) << tmpStream.str(); \
        xerus::misc::internal::bufferStreams.checkSwitch(); \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("error"))) {xerus::misc::internal::dump_log_buffer(std::string("error invoked:\n")+tmpStream.str()); }; \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("critical"))) {xerus::misc::internal::dump_log_buffer(std::string("critical error invoked:\n")+tmpStream.str()); }; \
        if (COMPILE_TIME_EVAL(xerus::misc::internal::log_namehash(STRINGIFY(lvl))==xerus::misc::internal::log_namehash("fatal"))) {xerus::misc::internal::dump_log_buffer(std::string("fatal error invoked:\n")+tmpStream.str()); }; 
#else // no log buffer
    #define NAMED_LOGGER_LOGBUFFER
#endif


#define LOG(lvl, ...) \
    if (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag != xerus::misc::internal::NOT_LOGGING) { \
        std::stringstream tmpStream; \
        time_t xerus_err_t=std::time(0); tm *xerus_err_ltm = localtime(&xerus_err_t); \
        tmpStream \
				<< std::right << (1900+xerus_err_ltm->tm_year) << '-' <<std::setw(2) << std::setfill('0') <<  xerus_err_ltm->tm_mon \
				<< '-' <<std::setw(2) << std::setfill('0') <<  xerus_err_ltm->tm_mday \
				<< ' ' << std::setw(2) << std::setfill('0') << xerus_err_ltm->tm_hour \
				<< ':' <<std::setw(2) << std::setfill('0') <<  xerus_err_ltm->tm_min \
				<< ':' <<std::setw(2) << std::setfill('0') <<  xerus_err_ltm->tm_sec << " " \
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

        
#define IS_LOGGING(lvl) \
    (XERUS_logFlag<xerus::misc::internal::log_namehash(STRINGIFY(lvl))>::flag != xerus::misc::internal::NOT_LOGGING)

    
SET_LOGGING_DEFAULT(xerus::misc::internal::LOGGING_FULL)
