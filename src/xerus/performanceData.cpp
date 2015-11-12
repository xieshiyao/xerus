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
#include <xerus/misc/missingFunctions.h>

namespace xerus {
	
	PerformanceData::Histogram::Histogram(const std::vector< xerus::PerformanceData::DataPoint >& _data, const xerus::value_t _base) : base(_base), totalTime(0) {
		for (size_t i = 1; i<_data.size(); ++i) {
			if (_data[i].residual >= _data[i-1].residual) {
				continue;
			}
			
			// assume x_2 = x_1 * 2^(-alpha * delta-t)
			value_t relativeChange = _data[i].residual/_data[i-1].residual;
			value_t exponent = log(relativeChange) / log(2);
			size_t delta_t = _data[i].elapsedTime - _data[i-1].elapsedTime;
			value_t rate = - exponent / value_t(delta_t);
			int logRate = int(log(rate)/log(base));
			buckets[logRate] += delta_t;
			totalTime += delta_t;
		}
	}
	
	PerformanceData::Histogram::Histogram(const xerus::value_t _base) : base(_base), totalTime(0) {}
	
	PerformanceData::Histogram PerformanceData::Histogram::operator+=(const Histogram &_other) {
		REQUIRE(misc::approx_equal(_other.base, base), "only histograms of identical base can be added");
		for (const auto &h : _other.buckets) {
			buckets[h.first] += h.second;
		}
		totalTime += _other.totalTime;
		return *this;
	}
	
	void PerformanceData::Histogram::dump_to_file(const std::string &_fileName) const {
		std::ofstream out(_fileName);
		for (auto &h : buckets) {
			out << pow(base, h.first) << " " << double(h.second)/double(totalTime) << '\n';
		}
		out.close();
	}
	
	void PerformanceData::add(const size_t _itrCount, const xerus::value_t _residual, const std::vector<size_t> _ranks, const size_t _flags) {
		if (active) {
			if (startTime == ~0ul) {
				start();
			}
			data.emplace_back(_itrCount, get_runtime(), _residual, _ranks, _flags);
			
			if(printProgress) {
				LOG(PerformanceData, "Iteration " << std::setw(4) << std::setfill(' ') << _itrCount 
					<< " Time: " <<  std::setw(6) << std::setfill(' ') << std::fixed << double(data.back().elapsedTime)*1e-6 
					<< "s Residual: " <<  std::setw(11) << std::setfill(' ') << std::scientific << _residual << " Flags: " << _flags << " Ranks: " << _ranks);
			}
		}
	}
	
	void PerformanceData::add(const xerus::value_t _residual, const TensorNetwork::RankTuple _ranks, const size_t _flags) {
		if (active) {
			if (data.empty()) {
				add(0, _residual, _ranks, _flags);
			} else {
				add(data.back().iterationCount+1, _residual, _ranks, _flags);
			}
		}
	}

	void PerformanceData::dump_to_file(const std::string &_fileName) const {
		std::string header;
		header += "# ";
		header += additionalInformation;
		misc::replace(header, "\n", "\n# ");
		header += "\n# \n#itr \ttime[us] \tresidual \tflags \tranks...\n";
		std::ofstream out(_fileName);
		out << header;
		for (const DataPoint &d : data) {
			out << d.iterationCount << '\t' << d.elapsedTime << '\t' << d.residual << '\t' << d.flags;
			for (size_t r : d.ranks) {
				out << '\t' << r;
			}
			out << '\n';
		}
		out.close();
	}
	
	PerformanceData::Histogram PerformanceData::get_histogram(const xerus::value_t _base) const {
		return Histogram(data, _base);
	}

	PerformanceData NoPerfData(false);
}
