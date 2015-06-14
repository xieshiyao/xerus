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

#include <xerus/contractionHeuristic.h>

namespace xerus {
    ContractionHeuristic::ContractionHeuristic(std::string _name, std::function<void(double &, std::vector<std::pair<size_t,size_t>> &, TensorNetwork &)> _scoreFct) 
        : name(_name), scoreFct(_scoreFct) {}
    
    double ContractionHeuristic::rescore(TensorNetwork _tn) { // NOTE take as value to get a deep copy instead of reference!
        score=0;
        contractions.clear();
        scoreFct(score, contractions, _tn);
        return score;
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - External Stuff - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    std::vector<ContractionHeuristic>* ContractionHeuristic::list = nullptr;
    
	static __attribute__((destructor)) void destroy_list() {
		delete ContractionHeuristic::list;
		ContractionHeuristic::list = nullptr;
	}

	//TODO non-quadratic
    #define GREEDY(name, alg) \
        void name(double &_score, std::vector<std::pair<size_t,size_t>> &_contractions, TensorNetwork &_tn) { \
        double best = std::numeric_limits<double>::max(); \
        size_t bestId1, bestId2; \
        do { \
            best = std::numeric_limits<double>::max(); \
            for (size_t i=0; i<_tn.nodes.size(); ++i) { \
                if (_tn.nodes[i].erased) continue; \
                TensorNetwork::TensorNode &ni = _tn.nodes[i]; \
                for (size_t j=i+1; j<_tn.nodes.size(); ++j) { \
                    if (_tn.nodes[j].erased) continue; \
                    TensorNetwork::TensorNode &nj = _tn.nodes[j]; \
                    /* possible candidate (i.e. link to a later node) */\
                    /* calculate n,m,r */ \
                    double m=1,n=1,r=1; \
                    for (size_t d=0; d<ni.degree(); ++d) { \
                        if (ni.neighbors[d].other == j) { \
                            r *= (double)ni.neighbors[d].dimension; \
                        } else { \
                            m *= (double)ni.neighbors[d].dimension; \
                        } \
                    } \
                    for (size_t d=0; d<nj.degree(); ++d) { \
                        if (nj.neighbors[d].other != i) { \
                            n *= (double)nj.neighbors[d].dimension; \
                        } \
                    } \
                    double tmpscore = alg; \
                    if (tmpscore < best) { \
                        best = tmpscore; \
                        bestId1 = i; \
                        bestId2 = j; \
                    } \
                } \
            } \
            if (best < std::numeric_limits<double>::max()) { \
                _score += _tn.contraction_cost(bestId1,bestId2); \
                _contractions.emplace_back(bestId1,bestId2); \
                _tn.contract(bestId1,bestId2); \
            } \
        } while (best < std::numeric_limits<double>::max()); \
    }

    

    namespace internal {
        GREEDY(greedy_size, n*m-(n+m)*r)
        GREEDY(greedy_speed, (n*m-(n+m)*r)/(n*m*r))
        GREEDY(greedy_r, -r)
        
        static ContractionHeuristic::AddToVector greedy_size_heuristic("greedy_size", greedy_size);
        static ContractionHeuristic::AddToVector greedy_speed_heuristic("greedy_speed", greedy_speed);
        static ContractionHeuristic::AddToVector greedy_rank_heuristic("greedy_r", greedy_r);
    }

}
