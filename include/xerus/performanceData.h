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
#include <map>
#include "basic.h"
#include "tensorNetwork.h"
#include "misc/timeMeasure.h"

namespace xerus {


	/// @brief Storage class for the performance data collected during an algorithm (typically iteration count, time and residual)
	class PerformanceData {
	public:
		struct DataPoint {
			size_t iterationCount;
			size_t elapsedTime;
			double residual;
			TensorNetwork::RankTuple ranks;
			size_t flags;
			
			DataPoint(const size_t _itrCount, const size_t _time, const value_t _res, const TensorNetwork::RankTuple _ranks, const size_t _flags) 
				: iterationCount(_itrCount), elapsedTime(_time), residual(_res), ranks(_ranks), flags(_flags) {}
		};
		
		struct Histogram {
			value_t base;
			std::map<int, size_t> buckets;
			size_t totalTime;
			
			///@brief Creates a histogram of convergence rates based on the data given. Assumes that the limit of the residuals is equal to 0!
			explicit Histogram(const std::vector<DataPoint> &_data, const value_t _base);
			
			explicit Histogram(const value_t _base);
			
			Histogram operator+=(const Histogram &_other);
			
			void dump_to_file(const std::string &_fileName) const;
		};
		
		const bool active;
		
		bool printProgress;
		
		size_t startTime;
		size_t stopTime;
		std::vector<DataPoint> data;
		
		std::string additionalInformation;
		
		
		explicit PerformanceData(const bool _printProgress = false, const bool _active = true) : active(_active), printProgress(_printProgress), startTime(~0ul), stopTime(~0ul) {}
		
		void start() {
			if (active) {
				startTime = misc::uTime();
			}
		}
		
		void stop_timer() {
			if (active) {
				stopTime = misc::uTime();
			}
		}
		
		void continue_timer() {
			if (active) {
				size_t currtime = misc::uTime();
				startTime += currtime - stopTime;
				stopTime = ~0ul;
			}
		}
		
		void reset() {
			if (active) {
				data.clear();
				additionalInformation.clear();
				startTime = ~0ul;
				stopTime = ~0ul;
			}
		}
		
		size_t get_elapsed_time() const {
			return misc::uTime() - startTime;
		}
		
		size_t get_runtime() const {
			if (stopTime != ~0ul) {
				return stopTime - startTime;
			} else {
				return misc::uTime() - startTime;
			}
		}
		
		void add(const size_t _itrCount, const value_t _residual, const TensorNetwork::RankTuple _ranks = TensorNetwork::RankTuple(), const size_t _flags = 0);
		
		void add(const value_t _residual, const TensorNetwork::RankTuple _ranks = TensorNetwork::RankTuple(), const size_t _flags = 0);
		
		operator bool() const {
			return active;
		}
		
		/// @brief The pipe operator allows to add everything that can be converted to string to the additional information in the header. 
		template<class T>
		PerformanceData& operator<<(const T &_info) noexcept {
			if (active) {
				std::string string = misc::to_string(_info);
				additionalInformation += string;
				
				if(printProgress) {
					LOG(PerformanceData, string);
				}
			}
			return *this;
		}
		
		void dump_to_file(const std::string &_fileName) const;
		
		Histogram get_histogram(const value_t _base) const;
	};

	extern PerformanceData NoPerfData;

}