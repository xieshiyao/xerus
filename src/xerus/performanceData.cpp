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
 * @brief Implementation of the PerformanceData class.
 */

#include <string>
#include <fstream>
#include <xerus/performanceData.h>

namespace xerus {
	
	void PerformanceData::add(size_t _itrCount, value_t _residual) {
		if (isLogging) {
			if (startTime == ~0ul) {
				start();
			}
			data.emplace_back(_itrCount, (misc::uTime()-startTime), _residual);
		}
	}
	
	void PerformanceData::add(value_t _residual) {
		if (data.empty()) {
			add(0, _residual);
		} else {
			add(data.back().iterationCount+1, _residual);
		}
	}

	void PerformanceData::dumpToFile(const std::string &_fileName) const {
		std::string header;
		header += "# ";
		header += additionalInformation;
		misc::replace(header, "\n", "\n# ");
		header += "\n# \n#itr \ttime[us] \tresidual\n";
		std::ofstream out(_fileName);
		out << header;
		for (const DataPoint &d : data) {
			out << d.iterationCount << " " << d.elapsedTime << " " << d.residual << '\n';
		}
		out.close();
	}

	PerformanceData NoPerfData(false);

}
