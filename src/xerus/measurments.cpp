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
 * @brief Implementation of the measurment classes class.
 */

#include <xerus/misc/check.h>
#include <xerus/measurments.h>
 
#include <xerus/misc/sort.h>

#include <xerus/index.h>
#include <xerus/tensor.h> 
#include <xerus/ttNetwork.h>
#include <xerus/indexedTensor.h>
 

namespace xerus {
	// --------------------- SinglePointMeasurmentSet -----------------
	
	size_t SinglePointMeasurmentSet::size() const {
		REQUIRE(positions.size() == measuredValues.size(), "I.E.");
		return positions.size();
	}
	
	size_t SinglePointMeasurmentSet::degree() const {
		IF_CHECK(
			for(size_t i = 0; i+1 < positions.size(); ++i) {
				REQUIRE(positions[i].size() == positions[i+1].size(), "Inconsitent degrees in measurment set.");
			}
		)
		return positions.empty() ? 0 : positions[0].size();
	}
	
	void SinglePointMeasurmentSet::add_measurment(const std::vector<size_t>& _position, const value_t _measuredValue) {
		positions.emplace_back(_position);
		measuredValues.emplace_back(_measuredValue);
	}
	
	
	void sort(SinglePointMeasurmentSet& _set, const size_t _splitPos) {
		misc::simultaneous_sort(_set.positions, _set.measuredValues, [_splitPos](const std::vector<size_t>& _lhs, const std::vector<size_t>& _rhs) {
		for (size_t i = 0; i < _splitPos && i < _lhs.size(); ++i) {
			if (_lhs[i] < _rhs[i]) return true;
			if (_lhs[i] > _rhs[i]) return false;
		}
		for (size_t i = _lhs.size(); i > _splitPos; --i) {
			if (_lhs[i-1] < _rhs[i-1]) return true;
			if (_lhs[i-1] > _rhs[i-1]) return false;
		}
		LOG(fatal, "Measurments must not appear twice. ");
		return false; // equality
	});
	}
	
	
	// --------------------- RankOneMeasurmentSet -----------------
	
	size_t RankOneMeasurmentSet::size() const {
		return measuredValues.size();
	}
	
	size_t RankOneMeasurmentSet::degree() const {
		return positions[0].size();
	}
	
	RankOneMeasurmentSet::RankOneMeasurmentSet(const SinglePointMeasurmentSet&  _other, const std::vector<size_t> _dimensions) {
		std::vector<Tensor> zeroPosition;
		for(size_t j = 0; j < +_other.degree(); ++j) {
			zeroPosition.emplace_back(Tensor({_dimensions[j]}));
		}
			
		for(size_t i = 0; i < _other.size(); ++i) {
			add_measurment(zeroPosition, _other.measuredValues[i]);
			for(size_t j = 0; j < _other.degree(); ++j) {
				positions.back()[j][_other.positions[i][j]] = 1.0;
			}
		}
	}
	
	void RankOneMeasurmentSet::add_measurment(const std::vector<Tensor>& _position, const value_t _measuredValue) {
		IF_CHECK(
			REQUIRE(positions.size() == measuredValues.size(), "Internal Error.");
			if(size() > 0) {
				for(size_t i = 0; i < degree(); ++i) {
					REQUIRE(positions[0][i].dimensions == _position[i].dimensions, "Inconsitend dimensions obtained");
				}
			}
			for (const Tensor & f : _position) {
				REQUIRE(f.degree() == 1, "illegal measurement");
			}
		);
		positions.emplace_back(_position);
		measuredValues.emplace_back(_measuredValue);
	}
	
	
	value_t RankOneMeasurmentSet::test_solution(const TTTensor& _solution) const {
		value_t residualNorm = 0.0;
		value_t measurementNorm = 0.0;
		
		const Index r1, r2, i1;
		
		std::vector<Tensor> reshuffledComponents(degree());
		
		for(size_t d = 0; d < degree(); ++d) {
			reshuffledComponents[d](i1, r1, r2) = _solution.get_component(d)(r1, i1, r2);
		}
		
		// TODO beautify
		#pragma omp parallel reduction(+: residualNorm, measurementNorm)
		{
			std::unique_ptr<Tensor[]> stackMem(new Tensor[degree()+1]);
			Tensor* const stack = stackMem.get()+1;
			stackMem[0] = Tensor::ones({1});
			
			for(size_t d = 0; d+1 < degree(); ++d) {
				stack[d].reset({_solution.rank(d)}, Tensor::Initialisation::None);
			}
			
			stack[degree()-1].reset({1}, Tensor::Initialisation::None);
			
			std::vector<Tensor> intermediates(degree());
			
			for(size_t d = 0; d < degree(); ++d) {
				intermediates[d].reset({_solution.get_component(d).dimensions[0], _solution.get_component(d).dimensions[2]}, Tensor::Initialisation::None);
			}
			
			REQUIRE(stack[degree()-1].size == 1 , "IE");
			
			#pragma omp for schedule(static)
			for(size_t i = 0; i < size(); ++i) {
				for(size_t d = 0; d < degree(); ++d) {
					contract(intermediates[d], positions[i][d], false, reshuffledComponents[d], false, 1);
					contract(stack[d], stack[d-1], false, intermediates[d], false, 1);
				}
				
				
				residualNorm += misc::sqr(measuredValues[i] - stack[degree()-1][0]);
				measurementNorm += misc::sqr(measuredValues[i]);
			}
		}
		
		return std::sqrt(residualNorm)/std::sqrt(measurementNorm);
	}
	
	
	void sort(RankOneMeasurmentSet& _set, const size_t _splitPos) {
		misc::simultaneous_sort(_set.positions, _set.measuredValues, [_splitPos](const std::vector<Tensor>& _lhs, const std::vector<Tensor>& _rhs) {
		for (size_t i = 0; i < _splitPos && i < _lhs.size(); ++i) {
			REQUIRE(_lhs[i].size == _rhs[i].size && _lhs[i].degree() == 1 && _rhs[i].degree() == 1, "");
			for(size_t j = 0; j < _lhs[i].size; ++j) {
				if (_lhs[i][j] < _rhs[i][j]) return true;
				if (_lhs[i][j] > _rhs[i][j]) return false;
			}
		}
		for (size_t i = _lhs.size(); i > _splitPos; --i) {
			REQUIRE(_lhs[i].size == _rhs[i].size && _lhs[i].degree() == 1 && _rhs[i].degree() == 1, "");
			for(size_t j = 0; j < _lhs[i].size; ++j) {
				if (_lhs[i-1][j] < _rhs[i-1][j]) return true;
				if (_lhs[i-1][j] > _rhs[i-1][j]) return false;
			}
		}
		LOG(fatal, "Measurments must not appear twice. ");
		return false; // equality
	});
	}
}
