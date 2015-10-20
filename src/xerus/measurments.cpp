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
 * @brief Implementation of the measurment classes class.
 */

#include <xerus/misc/check.h>
#include <xerus/measurments.h>
#include <xerus/misc/missingFunctions.h>

namespace xerus {
	// --------------------- SinglePointMeasurment -----------------
	
	SinglePointMeasurment::SinglePointMeasurment(const std::vector<size_t>& _positions, const value_t _value) : positions(_positions), value(_value) {}
	
	bool SinglePointMeasurment::Comparator::operator()(const SinglePointMeasurment &_lhs, const SinglePointMeasurment &_rhs) const {
		REQUIRE(_lhs.positions.size() == _rhs.positions.size(), "");
		for (size_t i = 0; i < split_position && i < _lhs.positions.size(); ++i) {
			if (_lhs.positions[i] < _rhs.positions[i]) return true;
			if (_lhs.positions[i] > _rhs.positions[i]) return false;
		}
		for (size_t i = _lhs.positions.size(); i > split_position; --i) {
			if (_lhs.positions[i-1] < _rhs.positions[i-1]) return true;
			if (_lhs.positions[i-1] > _rhs.positions[i-1]) return false;
		}
		return false; // equality
	}
	
	void sort(std::vector<SinglePointMeasurment>& _set, const size_t _splitPos) {
		SinglePointMeasurment::Comparator comp(_splitPos);
// 		auto perm = misc::create_sort_permutation(_set, comp);
// 		misc::apply_permutation(_set, perm);
		std::sort(_set.begin(), _set.end(), comp);
	}
	
	
	
	// --------------------- SinglePointMeasurmentSet -----------------
	
	SinglePointMeasurmentSet::SinglePointMeasurmentSet(const std::vector<SinglePointMeasurment>& _measurments) {
		for(const SinglePointMeasurment& meas : _measurments) {
			positions.emplace_back(meas.positions);
			measuredValues.emplace_back(meas.value);
		}
	}
	
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
	
	bool SinglePointMeasurmentSet::Comparator::operator()(const std::vector<size_t>& _lhs, const std::vector<size_t>& _rhs) const {
		REQUIRE(_lhs.size() == _rhs.size(), "");
		for (size_t i = 0; i < split_position && i < _lhs.size(); ++i) {
			if (_lhs[i] < _rhs[i]) return true;
			if (_lhs[i] > _rhs[i]) return false;
		}
		for (size_t i = _lhs.size(); i > split_position; --i) {
			if (_lhs[i-1] < _rhs[i-1]) return true;
			if (_lhs[i-1] > _rhs[i-1]) return false;
		}
		return false; // equality
	}
	
	void sort(SinglePointMeasurmentSet& _set, const size_t _splitPos) {
		SinglePointMeasurmentSet::Comparator comp(_splitPos);
		misc::simultaneous_sort(_set.positions, _set.measuredValues, comp);
	}
	
	
	// --------------------- RankOneMeasurmentSet -----------------
	
	void RankOneMeasurmentSet::add_measurment(const std::vector<FullTensor>& _position, const value_t _measuredValue) {
		positions.emplace_back(_position);
		measuredValues.emplace_back(_measuredValue);
	}
	
	bool RankOneMeasurmentSet::Comparator::operator()(const std::vector<FullTensor>& _lhs, const std::vector<FullTensor>& _rhs) const {
		REQUIRE(_lhs.size() == _rhs.size(), "I.E.");
		for (size_t i = 0; i < split_position && i < _lhs.size(); ++i) {
			REQUIRE(_lhs[i].size == _rhs[i].size && _lhs[i].degree() == 1 && _rhs[i].degree() == 1, "");
			for(size_t j = 0; j < _lhs[i].size; ++j) {
				if (_lhs[i][j] < _rhs[i][j]) return true;
				if (_lhs[i][j] > _rhs[i][j]) return false;
			}
		}
		for (size_t i = _lhs.size(); i > split_position; --i) {
			REQUIRE(_lhs[i].size == _rhs[i].size && _lhs[i].degree() == 1 && _rhs[i].degree() == 1, "");
			for(size_t j = 0; j < _lhs[i].size; ++j) {
				if (_lhs[i-1][j] < _rhs[i-1][j]) return true;
				if (_lhs[i-1][j] > _rhs[i-1][j]) return false;
			}
		}
		return false; // equality
	}
	
	void sort(RankOneMeasurmentSet& _set, const size_t _splitPos) {
		RankOneMeasurmentSet::Comparator comp(_splitPos);
		misc::simultaneous_sort(_set.positions, _set.measuredValues, comp);
	}
}
