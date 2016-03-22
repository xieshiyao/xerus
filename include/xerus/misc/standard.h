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
 * @brief Header file for global shorthand notations of elementary integer types and attribute lists.
 */

#pragma once

#include <cstdint>
#include <cstddef>

namespace xerus {
	/**
	 * @brief Collection of classes and functions that provide elementary functionality that is not special to xerus as a tensor library.
	 */
	namespace misc {
		/// @brief Namespace for function and classes designated only for internal use.
		namespace internal {}
	}

    // Shorter names for unsigned types
    typedef uint8_t byte; ///< unsigned int type of exactly 8 bit
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

}

//Make likely & unlikely paramters useable. (Probably totaly useless xD)
#ifdef __GNUC__
    #define likely(x)       __builtin_expect((x),1)
    #define unlikely(x)     __builtin_expect((x),0)
#else
    #define likely(x)       (x)
    #define unlikely(x)     (x)
#endif

// counterpart to the explicit keyword for constructors
/**
 * @def implicit
 * @brief Counterpart to explicit keyword. has no effect as 'implicit' is always implied as per c++11.
 */
#define implicit      

// all of these can be writen like [[gnu::unused]] but kdevelop will not recognize them then
/** 
 * @def _inline_ 
 * @brief Collection of attributes to force gcc to inline a specific function.
 */
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
#define _warn_unused_	__attribute__((warn_unused_result))
