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
 * @brief Header file for the PerformanceData class.
 */

#pragma once
#include <string>
#include <fstream>
#include "basic.h"
#include "misc/timeMeasure.h"

namespace xerus {


/// @brief storage class for the performance data collected during an algorithm (typically iteration count, time and residual)
class PerformanceData {
public:
	struct DataPoint {
		size_t iterationCount;
		size_t elapsedTime;
		value_t residual;
		DataPoint(size_t _itrCount, size_t _time, value_t _res) 
			: iterationCount(_itrCount), elapsedTime(_time), residual(_res) {}
	};
	std::string additionalInformation;
	std::vector<DataPoint> data;
	size_t startTime;
	bool isLogging;
	
	PerformanceData(bool logging=true) : startTime(~0ul), isLogging(logging) {}
	
	void start() {
		startTime = misc::uTime();
	}
	
	void reset() {
		data.clear();
		additionalInformation.clear();
		startTime = ~0ul;
	}
	
	void add(size_t _itrCount, value_t _residual);
	
	void add(value_t _residual);
	
	operator bool() const {
		return isLogging;
	}
	
	/// @brief The pipe operator allows to add everything that can be converted to string to the additional information in the header. 
	template<class T>
	PerformanceData& operator<<(const T &_info) noexcept {
		additionalInformation += misc::to_string(_info);
		return *this;
	}
	
	void dumpToFile(const std::string &_fileName) const;
};

extern PerformanceData NoPerfData;

}