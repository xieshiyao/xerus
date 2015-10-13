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
 * @brief Header file for some example measure sets for tensor completion.
 */

#pragma once
#include "../ttNetwork.h"

namespace xerus { namespace examples { namespace completion {
	
	/**
	 * @brief creates @a _numMeasurements measurements whose values are equal to @$ 1/\sqrt{\sum i_\mu^2} @$
	 * @details the corresponding tensor is only approximately low-rank, cf. Hackbusch, 2012 "Tensor spaces and numerical tensor calculus" and http://www.mis.mpg.de/scicomp/EXP_SUM
	 * @note that index (0,0,0,0...) is singular for @a _additiveConst = 0
	 */
	template<class random_engine, class distribution>
	std::vector<SinglePointMeasurment> inverse_index_norm(size_t _degree, size_t _n, size_t _numMeasurements, distribution _dist, random_engine _rnd, value_t _additiveConst=1.0) {
		REQUIRE(misc::pow(_n, _degree) >= _numMeasurements, "it's impossible to perform as many measurements as requested");
		std::uniform_int_distribution<size_t> indexDist(0,_n-1);
		std::set<std::vector<size_t>> measuredPositions;
		std::vector<SinglePointMeasurment> result;
		while (result.size() < _numMeasurements) {
			std::vector<size_t> position;
			value_t normSqr=0;
			for (size_t i=0; i<_n; ++i) {
				size_t id = indexDist(_rnd);
				position.emplace_back(id);
				normSqr += misc::sqr(double(id)+_additiveConst);
			}
			if (std::abs(normSqr) < 1 || measuredPositions.count(position) > 0) continue;
			
			measuredPositions.insert(position);
			result.emplace_back(position, 1/std::sqrt(normSqr));
		}
		return result;
	}
	
	
	
	/**
	 * @brief creates @a _numMeasurements measurements whose values are equal to @$ (\alpha + \sum i_\mu/(i_{\mu+1} + \alpha - 1))^{-1} @$
	 * @details the corresponding tensor is only approximately low-rank, cf. Grasedyck, Kluge, Kraemer 2015
	 */
	template<class random_engine, class distribution>
	std::vector<SinglePointMeasurment> inverse_index_ratios(size_t _degree, size_t _n, size_t _numMeasurements, distribution _dist, random_engine _rnd, value_t _additiveConst=1.0) {
		REQUIRE(misc::pow(_n, _degree) >= _numMeasurements, "it's impossible to perform as many measurements as requested");
		std::uniform_int_distribution<size_t> indexDist(0,_n-1);
		std::set<std::vector<size_t>> measuredPositions;
		std::vector<SinglePointMeasurment> result;
		while (result.size() < _numMeasurements) {
			std::vector<size_t> position;
			for (size_t i=0; i<_n; ++i) {
				size_t id = indexDist(_rnd);
				position.emplace_back(id);
			}
			if (measuredPositions.count(position) > 0) continue;
			
			value_t sum=0;
			for (size_t i=0; i<_n-1; ++i) {
				sum += (position[i] + 1) / (position[i+1] + _additiveConst);
			}
			
			measuredPositions.insert(position);
			result.emplace_back(position, 1/(_additiveConst + sum));
		}
		return result;
	}
	
}}}

