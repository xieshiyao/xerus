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


#ifdef TEST_

#include <map>
#include <string>
#include <functional>
#include "namedLogger.h"
#include "stringUtilities.h"


struct ___UnitTest {
	// order of contruction of global objects is random so the first one has to create the map
	static std::map<std::string, std::map<std::string, std::function<bool ()>>> *tests;
	
	___UnitTest(std::string _group, std::string _name, std::function<bool ()> _f) {
		if (!tests) {
			tests = new std::map<std::string, std::map<std::string, std::function<bool ()>>>();
		}
		if (tests->count(_group) > 0 && (*tests)[_group].count(_name) > 0) {
			LOG(error, "Unit test '" << _group << "::" << _name << "' defined multiple times!");
		}
		(*tests)[_group][_name] = _f;
	}
};

struct ___RequiredTest {
	struct identifier {
		std::string functionName;
		std::string filename;
		size_t lineNumber;
		identifier(std::string _func, std::string _file, size_t _line) 
			: functionName(_func), filename(_file), lineNumber(_line) {}
		bool operator<(const identifier &_rhs) const {
			if (functionName < _rhs.functionName) {
				return true;
			} else if (functionName == _rhs.functionName && lineNumber < _rhs.lineNumber) {
				return true;
			} else if (functionName == _rhs.functionName && lineNumber == _rhs.lineNumber && filename < _rhs.filename) {
				return true;
			} else {
				return false;
			}
		}
	};
	
	// order of contruction of global objects is random so the first one has to create the map
	static std::map<identifier, size_t> *tests;
	
	___RequiredTest(std::string _functionName, std::string _fileName, size_t _lineNb) 
	{
		if (!tests) {
			tests = new std::map<identifier, size_t>();
		}
		identifier key = identifier(_functionName, _fileName, _lineNb);
		if (tests->count(key) > 0) {
			LOG(warning, "test of function " << MISC::demangle_cxa(_functionName) << " file " << _fileName << " line " << _lineNb << " was required several times");
		}
		(*tests)[key] = 0;
	}
	
	static void increaseCounter(std::string _functionName, std::string _fileName, size_t _lineNb) {
		identifier key = identifier(_functionName, _fileName, _lineNb);
		if (tests->count(key) == 0) {
			LOG(fatal, "ie: unknown required test in " << MISC::demangle_cxa(_functionName) << " file " << _fileName << " line " << _lineNb);
		}
		(*tests)[key] += 1;
	}
};
#endif
