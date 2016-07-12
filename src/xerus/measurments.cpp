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
#include <xerus/misc/random.h>

#include <xerus/index.h>
#include <xerus/tensor.h> 
#include <xerus/tensorNetwork.h>
#include <xerus/ttNetwork.h>
#include <xerus/indexedTensor.h>
#include <xerus/misc/internal.h>
 

namespace xerus {
	// --------------------- SinglePointMeasurementSet -----------------
	
	SinglePointMeasurementSet SinglePointMeasurementSet::random(const size_t _numMeasurements, const std::vector<size_t>& _dimensions) {
		SinglePointMeasurementSet result;
		result.create_random_positions(_numMeasurements, _dimensions);
		result.measuredValues.resize(_numMeasurements, 0);
		return result;
	}
	
	
	SinglePointMeasurementSet SinglePointMeasurementSet::random(const size_t _numMeasurements, const Tensor& _solution) {
		SinglePointMeasurementSet result;
		result.create_random_positions(_numMeasurements, _solution.dimensions);
		result.measure(_solution);
		return result;
	}
		
	SinglePointMeasurementSet SinglePointMeasurementSet::random(const size_t _numMeasurements, const TensorNetwork& _solution) {
		SinglePointMeasurementSet result;
		result.create_random_positions(_numMeasurements, _solution.dimensions);
		result.measure(_solution);
		return result;
	}
	
	SinglePointMeasurementSet SinglePointMeasurementSet::random(const size_t _numMeasurements, const std::vector<size_t>& _dimensions, std::function<value_t(const std::vector<size_t>&)> _callback) {
		SinglePointMeasurementSet result;
		result.create_random_positions(_numMeasurements, _dimensions);
		result.measure(_callback);
		return result;
	}
	
	
	size_t SinglePointMeasurementSet::size() const {
		REQUIRE(positions.size() == measuredValues.size(), "Inconsitend SinglePointMeasurementSet encountered.");
		return positions.size();
	}
	
	
	size_t SinglePointMeasurementSet::degree() const {
		return positions.empty() ? 0 : positions[0].size();
	}
	
	
	void SinglePointMeasurementSet::add(std::vector<size_t> _position, const value_t _measuredValue) {
		REQUIRE(positions.empty() || _position.size() == positions.back().size(), "Given _position has incorrect degree " << _position.size() << ". Expected " << positions.back().size() << ".");
		positions.emplace_back(_position);
		measuredValues.emplace_back(_measuredValue);
	}
	
	
	value_t SinglePointMeasurementSet::frob_norm() const {
		const auto cSize = size();
		double norm = 0.0;
		for(size_t i = 0; i < cSize; ++i) {
			norm += misc::sqr(measuredValues[i]);
		}
		return std::sqrt(norm);
	}
	
	
	void SinglePointMeasurementSet::measure(const Tensor& _solution) {
		const auto cSize = size();
		for(size_t i = 0; i < cSize; ++i) {
			measuredValues[i] = _solution[positions[i]];
		}
	}
	
	// TODO not tested! --> Move sort to creation
	void SinglePointMeasurementSet::measure(const TensorNetwork& _solution) {
		REQUIRE(_solution.degree() == degree(), "Degrees of solution and measurements must match!");
		std::vector<TensorNetwork> stack(degree()+1);
		stack[0] = _solution;
		stack[0].reduce_representation();
		
		
		// Handle first measurment
		for(size_t i = 0; i < degree(); ++i) {
			stack[i+1] = stack[i];
			stack[i+1].fix_mode(0, positions[0][i]);
			stack[i+1].reduce_representation();
		}
		measuredValues[0] = stack.back()[0];
		
		const auto cSize = size();
		for(size_t j = 1; j < cSize; ++j) {
			REQUIRE(positions[j-1] != positions[j], "There were two identical measurements?");
			
			// Find the maximal recyclable stack position
			size_t rebuildIndex = 0;
			for(; rebuildIndex < degree(); ++rebuildIndex) {
				if(positions[j-1][rebuildIndex] != positions[j][rebuildIndex]) {
					break;
				}
			}
			
			// Rebuild stack
			for(size_t i = rebuildIndex; i < degree(); ++i) {
				stack[i+1] = stack[i];
				stack[i+1].fix_mode(0, positions[j][i]);
				stack[i+1].reduce_representation();
			}
			
			measuredValues[j] = stack.back()[0];
		}
	}
	
	
	void SinglePointMeasurementSet::measure(std::function<value_t(const std::vector<size_t>&)> _callback) {
		const auto cSize = size();
		for(size_t i = 0; i < cSize; ++i) {
			measuredValues[i] = _callback(positions[i]);
		}
	}
	
	
	double SinglePointMeasurementSet::test(const Tensor& _solution) const {
		const auto cSize = size();
		double error = 0.0, norm = 0.0;
		for(size_t i = 0; i < cSize; ++i) {
			error += misc::sqr(measuredValues[i] - _solution[positions[i]]);
			norm += misc::sqr(measuredValues[i]);
		}
		return std::sqrt(error/norm);
	}
	
	
	double SinglePointMeasurementSet::test(const TensorNetwork& _solution) const {
		const auto cSize = size();
		double error = 0.0, norm = 0.0;
		for(size_t i = 0; i < cSize; ++i) {
			error += misc::sqr(measuredValues[i] - _solution[positions[i]]);
			norm += misc::sqr(measuredValues[i]);
		}
		return std::sqrt(error/norm);
	}
	
	
	double SinglePointMeasurementSet::test(std::function<value_t(const std::vector<size_t>&)> _callback) const {
		const auto cSize = size();
		double error = 0.0, norm = 0.0;
		for(size_t i = 0; i < cSize; ++i) {
			error += misc::sqr(measuredValues[i] - _callback(positions[i]));
			norm += misc::sqr(measuredValues[i]);
		}
		return std::sqrt(error/norm);
	}
	
	
	
	void sort(SinglePointMeasurementSet& _set, const size_t _splitPos) {
		misc::simultaneous_sort(_set.positions, _set.measuredValues, [_splitPos](const std::vector<size_t>& _lhs, const std::vector<size_t>& _rhs) {
			for (size_t i = 0; i < _splitPos && i < _lhs.size(); ++i) {
				if (_lhs[i] < _rhs[i]) { return true; }
				if (_lhs[i] > _rhs[i]) { return false; }
			}
			
			for (size_t i = _lhs.size(); i > _splitPos; --i) {
				if (_lhs[i-1] < _rhs[i-1]) { return true; }
				if (_lhs[i-1] > _rhs[i-1]) { return false; }
			}
			return false; // equality
		});
	}
	
	
	void SinglePointMeasurementSet::create_random_positions(const size_t _numMeasurements, const std::vector<size_t>& _dimensions) {
		using ::xerus::misc::operator<<;
		XERUS_REQUIRE(misc::product(_dimensions) >= _numMeasurements, "It's impossible to perform as many measurements as requested. " << _numMeasurements << " > " << _dimensions);
		
		// Create distributions
		std::vector<std::uniform_int_distribution<size_t>> indexDist;
		for (size_t i = 0; i < _dimensions.size(); ++i) {
			indexDist.emplace_back(0, _dimensions[i]-1);
		}
		
		std::set<size_t> measuredPositions;
		std::vector<size_t> multIdx(_dimensions.size());
		while (positions.size() < _numMeasurements) {
			size_t pos = 0;
			for (size_t i = 0; i < _dimensions.size(); ++i) {
				multIdx[i] = indexDist[i](misc::randomEngine);
				pos *= _dimensions[i]; pos += multIdx[i];
			}
			if (!misc::contains(measuredPositions, pos)) {
				measuredPositions.insert(pos);
				positions.push_back(multIdx);
			}
		}
	}
	
	
	
	
	// --------------------- RankOneMeasurementSet -----------------
	
	size_t RankOneMeasurementSet::size() const {
		return measuredValues.size();
	}
	
	size_t RankOneMeasurementSet::degree() const {
		return positions[0].size();
	}
	
	RankOneMeasurementSet::RankOneMeasurementSet(const SinglePointMeasurementSet&  _other, const std::vector<size_t> &_dimensions) {
		std::vector<Tensor> zeroPosition;
		for(size_t j = 0; j < +_other.degree(); ++j) {
			zeroPosition.emplace_back(Tensor({_dimensions[j]}));
		}
			
		for(size_t i = 0; i < _other.size(); ++i) {
			add(zeroPosition, _other.measuredValues[i]);
			for(size_t j = 0; j < _other.degree(); ++j) {
				positions.back()[j][_other.positions[i][j]] = 1.0;
			}
		}
	}
	
	void RankOneMeasurementSet::add(const std::vector<Tensor>& _position, const value_t _measuredValue) {
		IF_CHECK(
			INTERNAL_CHECK(positions.size() == measuredValues.size(), "Internal Error.");
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
	
	
	value_t RankOneMeasurementSet::test_solution(const TTTensor& _solution) const {
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
			
			INTERNAL_CHECK(stack[degree()-1].size == 1 , "IE");
			
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
	
	
	void sort(RankOneMeasurementSet& _set, const size_t _splitPos) {
		misc::simultaneous_sort(_set.positions, _set.measuredValues, [_splitPos](const std::vector<Tensor>& _lhs, const std::vector<Tensor>& _rhs) {
			for (size_t i = 0; i < _splitPos && i < _lhs.size(); ++i) {
				REQUIRE(_lhs[i].size == _rhs[i].size && _lhs[i].degree() == 1 && _rhs[i].degree() == 1, "");
				for(size_t j = 0; j < _lhs[i].size; ++j) {
					if (_lhs[i][j] < _rhs[i][j]) { return true; }
					if (_lhs[i][j] > _rhs[i][j]) { return false; }
				}
			}
			for (size_t i = _lhs.size(); i > _splitPos; --i) {
				REQUIRE(_lhs[i].size == _rhs[i].size && _lhs[i].degree() == 1 && _rhs[i].degree() == 1, "");
				for(size_t j = 0; j < _lhs[i].size; ++j) {
					if (_lhs[i-1][j] < _rhs[i-1][j]) { return true; }
					if (_lhs[i-1][j] > _rhs[i-1][j]) { return false; }
				}
			}
			return false; // equality
		});
	}
} // namespace xerus
