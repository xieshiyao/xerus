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
 * @brief Implementation of some basic greedy contraction heuristics.
 */

#include <xerus/contractionHeuristic.h>

namespace xerus {
    namespace internal {
		
		template<double (*scoreFct)(double, double, double, double, double)>
		void greedy_heuristic(double &_bestCost, std::vector<std::pair<size_t,size_t>> &_contractions, TensorNetwork _network) {
			double bestScore, ourCost=0;
			double ourFinalCost=0;
			std::vector<std::pair<size_t,size_t>> ourContractions;
			size_t bestId1, bestId2;
			do {
				bestScore = std::numeric_limits<double>::max();
				for (size_t i=0; i<_network.nodes.size(); ++i) {
					if (_network.nodes[i].erased) continue;
					TensorNetwork::TensorNode &ni = _network.nodes[i];
					for (size_t j=i+1; j<_network.nodes.size(); ++j) {
						if (_network.nodes[j].erased) continue;
						TensorNetwork::TensorNode &nj = _network.nodes[j];
						/* possible candidate (i.e. link to a later node) */
						/* calculate n,m,r */
						double m=1,n=1,r=1;
						for (size_t d=0; d<ni.degree(); ++d) {
							if (ni.neighbors[d].other == j) {
								r *= (double)ni.neighbors[d].dimension;
							} else {
								m *= (double)ni.neighbors[d].dimension;
							}
						}
						for (size_t d=0; d<nj.degree(); ++d) {
							if (nj.neighbors[d].other != i) {
								n *= (double)nj.neighbors[d].dimension;
							}
						}
						double tmpscore = scoreFct(m,n,r,0.0,0.0);
						if (tmpscore < bestScore) {
							bestScore = tmpscore;
							ourCost = contraction_cost(m,n,r,0.0,0.0);
							bestId1 = i;
							bestId2 = j;
						}
					}
				}
				if (bestScore < std::numeric_limits<double>::max()) {
					ourFinalCost += ourCost;
					if (ourFinalCost > _bestCost) {
						return;
					}
					ourContractions.emplace_back(bestId1,bestId2);
					_network.contract(bestId1,bestId2);
				}
			} while (bestScore < std::numeric_limits<double>::max());
			if (ourFinalCost < _bestCost) {
				_bestCost = ourFinalCost;
				_contractions = std::move(ourContractions);
			}
		}
		
		
		
		
		double contraction_cost(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			return _m*_n*_r; // TODO sparse
		}
		
		
		
		
		double score_size(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			return _n*_m-(_n+_m)*_r;
		}
		double score_mn(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			return _m*_n;
		}
		double score_speed(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			return (_n*_m-(_n+_m)*_r)/(_n*_m*_r);
		}
		double score_r(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			return -_r;
		}
		double score_big_tensor(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			if (_n*_m<(_n+_m)*_r) {
				return -1e10 + _n*_m*_r;
			} else {
				return _n*_m-(_n+_m)*_r;
			}
		}
		double score_littlestep(double _m, double _n, double _r, double _sparsity1, double _sparsity2) {
			if (_n*_m<(_n+_m)*_r) {
				return -std::max(_n,_m)*_r;
			} else {
				return _n*_m-(_n+_m)*_r;
			}
		}
    }

}
