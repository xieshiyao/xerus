// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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

#include <set>
#include <random>

#include "basic.h"
 
#include "misc/math.h"

namespace xerus {
	class Tensor;
	template<bool isOperator> class TTNetwork;
	
	/** 
	* @brief Class used to represent a single point measurments.
	*/
    class SinglePointMeasurmentSet {
	public:
		std::vector<std::vector<size_t>> positions;
		std::vector<value_t> measuredValues;
		
		SinglePointMeasurmentSet() = default;
		SinglePointMeasurmentSet(const SinglePointMeasurmentSet&  _other) = default;
		SinglePointMeasurmentSet(      SinglePointMeasurmentSet&& _other) = default;
		
		template<class random_engine>
		static SinglePointMeasurmentSet random(const size_t _degree, const size_t _n, const size_t _numMeasurements, random_engine _rnd) {
			REQUIRE(misc::pow(_n, _degree) >= _numMeasurements, "It's impossible to perform as many measurements as requested");
			std::uniform_int_distribution<size_t> indexDist(0, _n-1);
			std::set<std::vector<size_t>> measuredPositions;
			SinglePointMeasurmentSet result;
			while (result.size() < _numMeasurements) {
				std::vector<size_t> pos;
				for (size_t i = 0; i < _degree; ++i) {
					pos.emplace_back(indexDist(_rnd));
				}
				if (measuredPositions.count(pos) > 0) continue;
				measuredPositions.insert(pos);
				result.add_measurment(pos, 0.0);
			}
			return result;
		}
		
		void add_measurment(const std::vector<size_t>& _position, const value_t _measuredValue);
		
		size_t size() const;
		
		size_t degree() const;
		
		struct Comparator {
			const size_t split_position;
			Comparator(const size_t _splitPos) : split_position(_splitPos)  {}
			bool operator()(const std::vector<size_t>& _lhs, const std::vector<size_t>& _rhs) const;
		};
	};
	
	void sort(SinglePointMeasurmentSet& _set, const size_t _splitPos = ~0ul);
	
	
	class RankOneMeasurmentSet {
	public:
		std::vector<std::vector<Tensor>> positions;
		std::vector<value_t> measuredValues;
		
		RankOneMeasurmentSet() = default;
		RankOneMeasurmentSet(const RankOneMeasurmentSet&  _other) = default;
		RankOneMeasurmentSet(      RankOneMeasurmentSet&& _other) = default;
		
		RankOneMeasurmentSet(const SinglePointMeasurmentSet&  _other, const std::vector<size_t> _dimensions);
		
		struct Comparator {
			const size_t split_position;
			Comparator(const size_t _splitPos) : split_position(_splitPos)  {}
			bool operator()(const std::vector<Tensor>& _lhs, const std::vector<Tensor>& _rhs) const;
		};
		
		void add_measurment(const std::vector<Tensor>& _position, const value_t _measuredValue);
		
		size_t size() const;
		
		size_t degree() const;
		
		value_t test_solution(const TTNetwork<false>& _solution) const;
	};
	
	void sort(RankOneMeasurmentSet& _set, const size_t _splitPos = ~0ul);
	
}
