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
    
    Index k,l,m,n,o,p;
	const size_t D = 10;
    
	const std::vector<size_t> stateDims(D, 2);
	
    TTTensor X = TTTensor::construct_random(stateDims, 4, rnd, dist);
	X /= X.frob_norm();
	FullTensor fullX(X);
	
	size_t posA = 0, posB = 0;
	
	for(size_t i = 1; i < fullX.size; ++i) {
		if(std::abs(fullX[i]) > std::abs(fullX[posA])) {
			posB = posA;
			posA = i;
		}
	}
	LOG(largestTwo, fullX[posA] << " and " << fullX[posB]);
	
	const double st = (std::abs(fullX[posA])-std::abs(fullX[posB]))/20;
// 	LOG(bla, "SoftT: " << st);
	
	while(true) {
		X = TTTensor::entrywise_product(X, X);
		X.soft_threshold(st);
		LOG(before, X.ranks() << " --- " << X.frob_norm());
		X /= X.frob_norm();
		LOG(after, X.ranks() << " --- " << X.frob_norm());
		
		if(misc::sum(X.ranks()) == D-1) {
			LOG(finished, "Wuhu");
			break;
		}
	}
)
