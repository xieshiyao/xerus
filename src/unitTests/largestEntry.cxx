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


#include<xerus.h>

#include "../../include/xerus/misc/test.h"
using namespace xerus;


UNIT_TEST(Algorithm, largestEntry,
    //Random numbers
    std::mt19937_64 rnd;
    rnd.seed(73);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	std::uniform_int_distribution<size_t> dimDist(1,4);
	std::uniform_int_distribution<size_t> rankDist(1,6);
    
	const size_t D = 16;
	
	for(size_t k = 0; k < 2; ++k) {
		std::vector<size_t> stateDims;
		stateDims.push_back(dimDist(rnd));
		
		std::vector<size_t> ranks;
		
		for(size_t d = 2; d <= D; ++d) {
			stateDims.push_back(dimDist(rnd));
			ranks.push_back(rankDist(rnd));
			REQUIRE(d == stateDims.size() && d == ranks.size()+1, "IE");
			
			TTTensor X = TTTensor::construct_random(stateDims, ranks, rnd, dist);
			X /= X.frob_norm();
			
			FullTensor fullX(X);
			
			size_t posA = 0, posB = 0;
			for(size_t i = 1; i < fullX.size; ++i) {
				if(std::abs(fullX[i]) >= std::abs(fullX[posA])) {
					posB = posA;
					posA = i;
				}
				else if(std::abs(fullX[i]) >= std::abs(fullX[posB])) {
					posB = i;
				}
			}
			LOG(largestEntry, "Largest entries are: " << fullX[posA] << " and " << fullX[posB] << " at " << posA << " and " << posB);
			
			double alpha = std::abs(fullX[posB]/fullX[posA]);
			double Xn = std::abs(fullX[posA]);
			double tau = (1-alpha)*alpha*Xn*Xn/(2.0*double(d-1));
			LOG(largestEntry, alpha << " >> " << tau);
			
			while(misc::sum(X.ranks()) >= d) {
				X = TTTensor::entrywise_product(X, X);
				LOG(largestEntry, "Before ST: " << X.ranks() << " --- " << X.frob_norm());
				X.soft_threshold(tau, true);
				LOG(largestEntry, "After ST: " << X.ranks() << " --- " << X.frob_norm());
				
				Xn = std::max(0.0, (1-(1-alpha)*alpha/2.0))*Xn*Xn;
				double fNorm = X.frob_norm();
				Xn /= fNorm;
				X /= fNorm;
				tau = (1-alpha)*alpha*Xn*Xn/(2.0*double(d-1));
			}
			
			size_t position = 0;
			size_t factor = misc::product(stateDims);
			for(size_t c = 0; c < d; ++c) {
				factor /= stateDims[c];
				size_t maxPos = 0;
				for(size_t i = 1; i < stateDims[c]; ++i) {
					if(std::abs(X.get_component(c)[i]) > std::abs(X.get_component(c)[maxPos])) {
						maxPos = i;
					}
				}
				position += maxPos*factor;
			}
			
			LOG(largestEntry, "Result: " << fullX[position] << " vs " << fullX[posA] << " at positions " << position << " and " << posA);
			TEST(position == posA);
			if(position != posA) {
				LOG(omg, fullX.to_string());
			}
		}
	}
)
