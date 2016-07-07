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
 * @brief Header file for some elementary string manipulation routines.
 */

#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <utility>

#include "standard.h"
#include "sfinae.h"
#include "containerOutput.h"

namespace xerus {
    namespace misc {

        /// @brief Demangles the function and class names created by gcc into a more readable format.
        std::string demangle_cxa(const std::string &_cxa);

        /// @brief Resolves 'folder/..' occurences in pathnames.
        __attribute__((const, pure)) std::string normalize_pathname(const std::string &_name);

        ///@brief: Explodes a string at positions indicated by _delim.
        __attribute__((const, pure)) std::vector<std::string> explode(const std::string& _string, const char _delim);

        ///@brief: Replaces all occurences of _search in _string by _replace.
        void replace(std::string& _string, const std::string& _search, const std::string& _replace);

        ///@brief: Removes all leading and trailing whitespaces from _string.
        void trim(std::string& _string, const std::string& whitespace = " \t\n\r\v");
		
		///@brief: Removes all leading and trailing whitespaces from _string.
		std::string trim(const std::string& _string, const std::string& whitespace = " \t\n\r\v");

        ///@brief: Removes all leading and trailing whitespaces from _string, and reduces all double whitespaces to one.
        void reduce(std::string& _string, const std::string& whitespace = " \t\n\r\v", const std::string& fill = " ");


        // Enable all standard to string functions also in xerus misc namespace
        using std::to_string;
        
        static XERUS_force_inline std::string to_string(const bool obj) {
            return obj ? std::string("TRUE") : std::string("FALSE");
        }
        
        static XERUS_force_inline std::string to_string(const std::string& obj) {
            return obj;
        }

        static XERUS_force_inline std::string to_string(const char * const obj) {
            return std::string(obj);
        }
        
//         #if __GNUC__ > 4 || defined(__clang__)
        
            XERUS_GENERATE_EXISTS_FUNCTION(to_string)
            
            ///@brief: Converts any object (which has no specialized function for it) to string.
            template<typename T, typename std::enable_if<!sfinae::exists_to_string<T>::value, int>::type = 0>
            std::string to_string(const T& obj) {
                std::ostringstream stream;
                stream << obj;
                return stream.str();
            }
            
//         #else
//             
//             ///@brief: Converts any object (which has no specialized function for it) to string.
//             template<typename T>
//             std::string to_string(const T& obj) {
//                 std::ostringstream stream;
//                 stream << obj;
//                 return stream.str();
//             }
//             
//         #endif
        
        
        ///@brief: Converts an arbitary Object to string with fixed precision
        template<typename T>
        std::string to_string(const T& obj, const long _precision) {
            std::ostringstream stream;
            stream.precision(_precision);
            stream << std::fixed << obj;
            return stream.str();
        }

        ///@brief: Creates an arbitary Object from string
        template<typename T>
        T from_string(const std::string& _s){
            std::stringstream stream(_s);
            T tmp;
            stream >> tmp;
            return tmp;
        }

    }
}
