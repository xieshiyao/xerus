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
    
	const size_t D = 26;
	const size_t N = 2;
    
	const std::vector<size_t> stateDims(D, N);
	
    TTTensor X = TTTensor::construct_random(stateDims, 6, rnd, dist);
	X /= X.frob_norm();
	
	FullTensor fullX(X);
	
	size_t posA = 0, posB = 0;
	LOG(creation, "finished");
	for(size_t i = 1; i < fullX.size; ++i) {
		if(std::abs(fullX[i]) > std::abs(fullX[posA])) {
			posB = posA;
			posA = i;
		}
	}
	LOG(largestTwo, fullX[posA] << " and " << fullX[posB]);
	
	double alpha = std::abs(fullX[posA]/fullX[posB]);
	double beta = std::sqrt(alpha);
	double Xb = std::abs(fullX[posB]);
	double st = (beta - 1)*beta*beta*beta*Xb*Xb/(2*(D-1));
	LOG(bsdgfla, alpha << " >> " << beta << " >> " << st);
	
	while(misc::sum(X.ranks()) >= D) {
		X = TTTensor::entrywise_product(X, X);
		LOG(beforeST, X.ranks() << " --- " << X.frob_norm());
		X.soft_threshold(st, true);
		LOG(afterST, X.ranks() << " --- " << X.frob_norm());
		alpha = std::abs(X[posA]/X[posB]);
		beta = std::sqrt(alpha);
		Xb = std::abs(X[posB]);
		st = (beta - 1)*beta*beta*beta*Xb*Xb/(2*(D-1));
	
		LOG(largestTwo, X[posA] << " and " << X[posB] << " next ST " << (D-1)*st << " from approx Xb: " << Xb);
	}
	
	size_t position = 0;
	size_t factor = misc::pow(N, D-1);
	for(size_t c = 0; c < D; ++c) {
		size_t maxPos = 0;
		for(size_t i = 1; i < N; ++i) {
			if(std::abs(X.get_component(c)[i]) > std::abs(X.get_component(c)[maxPos])) {
				maxPos = i;
			}
		}
		position += maxPos*factor;
		factor /= N;
	}
	
	LOG(ResultEntry, fullX[position]);
	TEST(position == posA);
)
