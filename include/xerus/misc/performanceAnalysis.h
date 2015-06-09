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

#include "standard.h"
#include <string>

#ifdef PERFORMANCE_ANALYSIS
	#include <map>
	#include <utility>
	#include "timeMeasure.h"
	#define PA_START size_t pa_startTime = misc::uTime()
	#define PA_END(group, name, parameter) { std::pair<size_t, size_t>& pa_call = xerus::misc::performanceAnalysis::calls[group][name][parameter]; pa_call.first++; pa_call.second += misc::uTime() - pa_startTime; }
#else 
	#define PA_START 
	#define PA_END(group, name, parameter)
#endif

namespace xerus {
	namespace misc {
		/// This namespace contains all functions only used for the performance analysism, as well as the global variables.
		namespace performanceAnalysis {
			#ifdef PERFORMANCE_ANALYSIS
				extern const size_t startupTime;
				extern std::map<std::string, std::map<std::string, std::map<std::string, std::pair<size_t, size_t>>>> calls;
			#endif
			
			/// Returns a detailed performance analysis if PERFORMANCE_ANALYSIS is set, an emtpy string otherwise
			std::string get_analysis();
		}
	}
}

	