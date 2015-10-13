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
 * @brief Header file for the different measurment classes.
 */

#pragma once

#include "basic.h"
#include "misc/missingFunctions.h"
#include "fullTensor.h"

namespace xerus {
	/** 
	* @brief Class used to represent a single point measurments.
	*/
    struct SinglePointMeasurment {
		std::vector<size_t> positions;
		value_t value;
		
		SinglePointMeasurment() = default;
		
		SinglePointMeasurment(const std::vector<size_t>& _positions, const value_t _value);
		
		struct Comparator {
			const size_t split_position;
			Comparator(const size_t _splitPos = ~0ul) : split_position(_splitPos)  {}
			bool operator()(const SinglePointMeasurment &_lhs, const SinglePointMeasurment &_rhs) const;
		};
		
		template<class random_engine, class distribution>
		static std::vector<SinglePointMeasurment> create_set(size_t _degree, size_t _n, size_t _numMeasurements, distribution _dist, random_engine _rnd) {
			REQUIRE(misc::pow(_n, _degree) >= _numMeasurements, "it's impossible to perform as many measurements as requested");
			std::uniform_int_distribution<size_t> indexDist(0,_n-1);
			std::set<std::vector<size_t>> measuredPositions;
			std::vector<SinglePointMeasurment> result;
			while (result.size() < _numMeasurements) {
				std::vector<size_t> pos;
				for (size_t i=0; i<_degree; ++i) {
					pos.emplace_back(indexDist(_rnd));
				}
				if (measuredPositions.count(pos) > 0) continue;
				measuredPositions.insert(pos);
				result.emplace_back(pos, 0.0);
			}
			return result;
		}
	};
	
	void sort(std::vector<SinglePointMeasurment>& _set, const size_t _splitPos = ~0ul);
	
	class RankOneMeasurmentSet {
	public:
		std::vector<std::vector<FullTensor>> positions;
		std::vector<value_t> measuredValues;
		
		RankOneMeasurmentSet() = default;
		RankOneMeasurmentSet(const RankOneMeasurmentSet&  _other) = default;
		RankOneMeasurmentSet(      RankOneMeasurmentSet&& _other) = default;
		
		struct Comparator {
		const size_t split_position;
		Comparator(const size_t _splitPos) : split_position(_splitPos)  {}
		bool operator()(const std::vector<FullTensor>& _lhs, const std::vector<FullTensor>& _rhs) const;
	};
		
		void add_measurment(const std::vector<FullTensor>& _position, const value_t _measuredValue);
	};
	
	void sort(RankOneMeasurmentSet& _set, const size_t _splitPos = ~0ul);
}
