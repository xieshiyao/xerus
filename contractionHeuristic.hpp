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

#pragma once

#include "xerus.h"

namespace xerus {
    std::vector<ContractionHeuristic>* ContractionHeuristic::list = nullptr;
	static __attribute__((destructor)) void destroy_list() {
		delete ContractionHeuristic::list;
		ContractionHeuristic::list = nullptr;
	}

	//TODO non-quadratic
#define GREEDY(name, alg) \
	void name(float &_score, std::vector<std::pair<size_t,size_t>> &_contractions, TensorNetwork &_tn) { \
	float best = 1e32f; \
	size_t bestId1, bestId2; \
	do { \
		best = 1e32f; \
		for (size_t i=0; i<_tn.nodes.size(); ++i) { \
			if (_tn.nodes[i].erased) continue; \
			TensorNode &ni = _tn.nodes[i]; \
			for (size_t j=i+1; j<_tn.nodes.size(); ++j) { \
				if (_tn.nodes[j].erased) continue; \
				TensorNode &nj = _tn.nodes[j]; \
				/* possible candidate (i.e. link to a later node) */\
				/* calculate n,m,r */ \
				float m=1,n=1,r=1; \
				for (size_t d=0; d<ni.degree(); ++d) { \
					if (ni.neighbors[d].other == j) { \
						r *= (float)ni.neighbors[d].dimension; \
					} else { \
						m *= (float)ni.neighbors[d].dimension; \
					} \
				} \
				for (size_t d=0; d<nj.degree(); ++d) { \
					if (nj.neighbors[d].other != i) { \
						n *= (float)nj.neighbors[d].dimension; \
					} \
				} \
				float tmpscore = alg; \
				if (tmpscore < best) { \
					best = tmpscore; \
					bestId1 = i; \
					bestId2 = j; \
				} \
			} \
		} \
		if (best < 1e32f) { \
			_score += (float)_tn.contraction_cost(bestId1,bestId2); \
			_contractions.emplace_back(bestId1,bestId2); \
			_tn.contract(bestId1,bestId2); \
		} \
	} while (best < 1e32f); \
}

GREEDY(greedy_size, n*m-(n+m)*r)
GREEDY(greedy_speed, (n*m-(n+m)*r)/(n*m*r))
GREEDY(greedy_r, -r)
// GREEDY(greedy_bs, r)

namespace {
ContractionHeuristic::AddToVector g("greedy_size", greedy_size);
ContractionHeuristic::AddToVector s("greedy_speed", greedy_speed);
ContractionHeuristic::AddToVector r("greedy_r", greedy_r);
// ContractionHeuristic::AddToVector bs("greedy_bullshit", greedy_bs);
}

}
