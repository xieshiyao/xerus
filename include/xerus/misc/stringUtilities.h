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
 * @brief Header file for some elementary string manipulation routines.
 */

#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <utility>

#include "standard.h"

namespace xerus {
    namespace misc {

        /// @brief Demangles the function and class names created by gcc into a more readable format.
        std::string demangle_cxa(const std::string &_cxa);

        /// @brief Resolves 'folder/..' occurences in pathnames.
        _const_ std::string normalize_pathname(const std::string &_name);

        ///@brief: Explodes a string at positions indicated by _delim.
        _const_ std::vector<std::string> explode(const std::string& _string, const char _delim);

        ///@brief: Replaces all occurences of _search in _string by _replace.
        void replace(std::string& _string, const std::string& _search, const std::string& _replace);

        ///@brief: Removes all leading and trailing whitespaces from _string.
        void trim(std::string& _string, const std::string& whitespace = " \t\n\r\v");

        ///@brief: Removes all leading and trailing whitespaces from _string, and reduces all double whitespaces to one.
        void reduce(std::string& _string, const std::string& whitespace = " \t\n\r\v", const std::string& fill = " ");


        ///@brief: Converts an arbitary Object to string
        template<typename T>
        static std::string to_string(const T& obj) {
            std::stringstream stream;
            stream.precision(8);
            stream << std::fixed << obj;
            return stream.str();
        }

        // Avoid stringstream if possible
        //TODO following with SFINAE (if std::to_string exists...)
        template<>
        _inline_ std::string to_string<int>(const int& obj) {
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<unsigned>(const unsigned& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<long>(const long& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<unsigned long>(const unsigned long& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<long long>(const long long& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<unsigned long long>(const unsigned long long& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<float>(const float& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<double>(const double& obj) {
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<long double>(const long double& obj){
            return std::to_string(obj);
        }

        template<>
        _inline_ std::string to_string<std::string>(const std::string& obj){
            return obj;
        }

        static _inline_ std::string to_string(const char * obj) {
            return std::string(obj);
        }

        ///@brief: Converts an arbitary Object to string with fixed precision
        template<typename T>
        std::string to_string(const T& obj, const size_t _precision){
            std::stringstream stream;
            stream.precision((long) _precision);
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
