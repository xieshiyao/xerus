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

#include <exception>
#include <sstream>

namespace xerus {
    namespace misc {
        class generic_error : public std::exception {
        protected:
            std::stringstream error_info;
            
        public:
            generic_error();
            
            generic_error(const generic_error &_other) noexcept;
            
            const char* what() const noexcept override;
            
            template<class T>
            generic_error& operator<<(const T &_info) noexcept {
                error_info << _info;
                return *this;
            }
        };
    }
}

#define XERUS_THROW(...) throw (__VA_ARGS__ << "\nexception thrown in function: " << (__func__) << " (" << (__FILE__) <<" : " << (__LINE__) << ")\n")
