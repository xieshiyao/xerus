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

#include <map>
#include <string>
#include <functional>
#include "namedLogger.h"

#ifdef TEST_
struct ___UnitTest {
	static uint uses; // indirection over 'nifty counter' because elementary data types are created ar compile time and 
					  // order of initialization of global objects is random
	static std::map<std::string, std::map<std::string, std::function<bool ()>>> *tests;
	
	___UnitTest(std::string _group, std::string _name, std::function<bool ()> _f) {
		if (0 == uses++) {
			tests = new std::map<std::string, std::map<std::string, std::function<bool ()>>>();
		}
		if (tests->count(_group) > 0 && (*tests)[_group].count(_name) > 0) {
			LOG(error, "Unit test '" << _group << "::" << _name << "' defined multiple times!");
		}
		(*tests)[_group][_name] = _f;
	}
};
#endif
