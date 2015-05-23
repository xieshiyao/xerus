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

#include <cstdint>
#include <cstddef>

#ifdef MISC_NAMESPACE
    #define START_MISC_NAMESPACE namespace MISC_NAMESPACE {
    #define END_MISC_NAMESPACE }
    #define MISC MISC_NAMESPACE
#else
    #define START_MISC_NAMESPACE
    #define END_MISC_NAMESPACE 
    #define MISC 
#endif

START_MISC_NAMESPACE

// Shorter names for unsigned types
typedef uint8_t byte;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

// Shorter names for fixed width types
typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

END_MISC_NAMESPACE

//Make likely & unlikely paramters useable. (Probably totaly useless xD)
#ifdef __GNUC__
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

// For formatting it is somethimes nice to declare constructors implicit 
#define implicit      

// very thorough always_inline version
// all of these can be writen like [[gnu::unused]] but kdevelop will not recognize them then
#define _inline_  		__attribute__((always_inline, gnu_inline)) inline

#define _noinline_ 		__attribute__((noinline))
#define _flatten_ 		__attribute__((flatten))
#define _const_ 		__attribute__((const, pure))
#define _pure_ 			__attribute__((pure))
#define _deprecated_ 	__attribute__((deprecated))
#define _noreturn_ 		__attribute__((noreturn))
#define _hot_ 			__attribute__((hot))
#define _cold_ 			__attribute__((cold))
#define _unused_ 		__attribute__((unused))

