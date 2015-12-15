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
 * @brief Implementation of the ALS variants.
 */

#include <xerus/algorithms/als.h>
#include <xerus/basic.h>
#include <xerus/tensorNetwork.h>
 
#include <xerus/indexedTensor_tensor_factorisations.h>

namespace xerus {
	///@brief Creates a set of @a _num unique random index tuples, using the @a _dimensions starting at @a _position +1.
	std::vector<std::vector<size_t>> create_random_tuples(const std::vector<size_t>& _dimensions, const size_t _position, const size_t _num) {
		const size_t N = _dimensions.size()-_position;
		
		std::random_device rd;
		std::mt19937_64 rnd(rd());
		std::vector<std::uniform_int_distribution<size_t>> dists;
		for(size_t d = 0; d < N; ++d) {
			dists.emplace_back(0, _dimensions[_position+d]);
		}

		std::vector<std::vector<size_t>> tuples(_num, std::vector<size_t>(N, 0));
		
		for(size_t i = 0; i < _num; ++i) {
			// Create a random tuple.
			for(size_t d = 0; d < N; ++d) {
				tuples[i][d] = dists[d](rnd);
			}
			
			// Redo this one if it is non unique.
			for(size_t j = 0; j < i; ++j) {
				if(tuples[i] == tuples[j]) { --i; break; }
			}
		}
		
		return tuples;
	}
	
	TTTensor cross_approximation(const Tensor& _input, const std::vector<size_t> _ranks) {
		std::random_device rd;
		std::mt19937_64 rnd(rd());
		
		TTTensor reconstruction(_input.dimensions.size());
		
		
		std::vector<std::vector<size_t>> leftTuples;
		
		reconstruction.set_component(0, Tensor({1, _input.dimensions[0], _ranks[0]}, Tensor::Representation::Dense, Tensor::Initialisation::None));
		std::vector<std::vector<size_t>> tuples = create_random_tuples(_input.dimensions, 0, _ranks[0]);
		
		for(size_t i = 0; i < _input.dimensions[0]; ++i) {
			for(size_t j = 0; j < tuples.size(); ++j) {
				std::vector<size_t> index(1, i);
				index.insert(index.end(), tuples[j].begin(), tuples[j].end());
				reconstruction.component(0)[{0,i,j}] = _input[index];
			}
		}
		
		std::uniform_int_distribution<size_t> dist(0, _input.dimensions[0]);
		
		for(size_t i = 0; i < _ranks[0]; ++i) {
			leftTuples.emplace_back(dist(rnd));
		}
		
		for(size_t position = 0; position < _input.degree()-1; ++position) {
// 			std::vector<std::vector<size_t>> tuples = create_random_tuples(_input.dimensions, position, _ranks[position]);
			
// 			reconstruction.set_component(position, Tensor({_ranks[position-1], _input.dimensions[position], _ranks[position]}, Tensor::Initialisation::None));
		}
		
		return reconstruction;
	}
	
}
