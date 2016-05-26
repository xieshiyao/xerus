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
#include "misc/containerSupport.h"

namespace xerus {
	class Tensor;
	class TensorNetwork;
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
		
		SinglePointMeasurmentSet& operator=(const SinglePointMeasurmentSet&  _other) = default;
		SinglePointMeasurmentSet& operator=(      SinglePointMeasurmentSet&& _other) = default;
		
		template<class random_engine>
		static SinglePointMeasurmentSet random(const std::vector<size_t> &_dim, const size_t _numMeasurements, random_engine _rnd) {
			REQUIRE(misc::product(_dim) >= _numMeasurements, "It's impossible to perform as many measurements as requested. " << _numMeasurements << " > " << _dim);
			std::vector<std::uniform_int_distribution<size_t>> indexDist;
			for (size_t i=0; i<_dim.size(); ++i) {
				indexDist.emplace_back(0, _dim[i]-1);
			}
			std::set<std::vector<size_t>> measuredPositions;
			SinglePointMeasurmentSet result;
			while (result.size() < _numMeasurements) {
				std::vector<size_t> pos;
				for (size_t i = 0; i < _dim.size(); ++i) {
					pos.emplace_back(indexDist[i](_rnd));
				}
				if (measuredPositions.count(pos) > 0) continue;
				measuredPositions.insert(pos);
				result.add(pos, 0.0);
			}
			return result;
		}
		
		void add(const std::vector<size_t>& _position, const value_t _measuredValue);
		
		size_t size() const;
		
		size_t degree() const;
		
		value_t test_solution(const TensorNetwork& _solution) const;
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

		RankOneMeasurmentSet& operator=(const RankOneMeasurmentSet&  _other) = default;
		RankOneMeasurmentSet& operator=(      RankOneMeasurmentSet&& _other) = default;
		
		void add(const std::vector<Tensor>& _position, const value_t _measuredValue);
		
		size_t size() const;
		
		size_t degree() const;
		
		value_t test_solution(const TTNetwork<false>& _solution) const;
	};
	
	void sort(RankOneMeasurmentSet& _set, const size_t _splitPos = ~0ul);
	
}
