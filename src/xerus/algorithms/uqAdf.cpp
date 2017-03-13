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
 * @brief Implementation of the ADF variants. 
 */

#include <xerus/algorithms/uqAdf.h>
 
 #include <xerus/indexedTensorMoveable.h>
 #include <xerus/misc/basicArraySupport.h>
 #include <xerus/misc/simpleNumerics.h>
#include <xerus/misc/internal.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace xerus {
    
    Tensor randVar_to_position(const double _v, const size_t _polyDegree) {
        const std::vector<xerus::misc::Polynomial> stochasticBasis = xerus::misc::Polynomial::build_orthogonal_base(_polyDegree, [](const double){return 1.0;}, -1., 1.);
        
        Tensor p({stochasticBasis.size()});
        for (size_t i = 0; i < stochasticBasis.size(); ++i) {
            p[i] = stochasticBasis[i](_v);
        }
        
        return p;
    }
    
    
    
    void uq_adf(TTTensor& _x, const std::vector<std::vector<double>>& _randomVariables, const std::vector<Tensor>& _solutions, const size_t _polyDegree) {
        REQUIRE(_randomVariables.size() == _solutions.size(), "ERROR");
        const size_t N = _randomVariables.size();
        const size_t d = _x.degree();
		
		Tensor one = Tensor::ones({1});
        
        std::vector<std::vector<Tensor>> leftStack(d, std::vector<Tensor>(N)); // From corePosition 1 to d-2
        std::vector<std::vector<Tensor>> rightStack(d, std::vector<Tensor>(N));  // From corePosition 1 to d-1
        
        std::vector<Tensor> Ax;
        
        // Rebuild right Stack
		LOG(bug, "Rebuild right stack");
        for(size_t corePosition = _x.degree()-1; corePosition > 1; --corePosition) {
			LOG(bug, "corePosition: " << corePosition);
            Tensor tmp;
			
			const Tensor shuffledX = reshuffle(_x.get_component(corePosition), {1, 0, 2});
			
            for(size_t j = 0; j < N; ++j) {
// 				LOG(bug, randVar_to_position(_randomVariables[j][corePosition] ,_polyDegree).dimensions << " vs " << shuffledX.dimensions);
				contract(tmp, randVar_to_position(_randomVariables[j][corePosition] ,_polyDegree), shuffledX, 1);
				
				if(corePosition < d-1) {
// 					LOG(bug, rightStack[corePosition+1][j].dimensions << " vs " << tmp.dimensions);
					contract(rightStack[corePosition][j], tmp, rightStack[corePosition+1][j], 1);
				} else {
// 					LOG(bug, one.dimensions << " vs " << tmp.dimensions);
					contract(rightStack[corePosition][j], tmp, one, 1);
				}
            }
        }
        
        
        for(size_t corePosition = 1; corePosition < _x.degree(); ++corePosition) {
			LOG(bug, "Rebuild left stack @ " << corePosition);
            Tensor tmp;
			const Tensor shuffledX = reshuffle(_x.get_component(corePosition), {1, 0, 2});
			
            for(size_t j = 0; j < N; ++j) {
				contract(tmp, randVar_to_position(_randomVariables[j][corePosition], _polyDegree), shuffledX, 1);
				if(corePosition > 1) {
					contract(leftStack[corePosition][j], leftStack[corePosition-1][j], tmp, 1);
				} else if(corePosition==1) {
					leftStack[corePosition][j] = tmp;
				}
            }
        }
        
        LOG(bla, "Hallo Welt");
    }
	
} // namespace xerus
