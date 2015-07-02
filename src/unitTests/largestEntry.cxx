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

UNIT_TEST(Test, bla,
	const size_t N = 2000;
	Index i,j,k,l;
	FullTensor U({N,N}, [](const std::vector<size_t>& _pos) {
		if(_pos[0] == 0) {return _pos[1] == 0 ? -1.0/std::sqrt(N) : -1/std::sqrt(_pos[1]*(_pos[1]+1)); }
		else if(_pos[1] == 0) { return -1.0/std::sqrt(N); }
		else if(_pos[0] == _pos[1]) { return double(_pos[0])/std::sqrt(_pos[1]*(_pos[1]+1)); }
		else if(_pos[0] < _pos[1]) { return -1.0/std::sqrt(_pos[1]*(_pos[1]+1)); }
		else { return 0.0; }});
	U(i,j) = -1.0*U(j,i);
	FullTensor S({N,N}, [](const std::vector<size_t>& _pos){return _pos[0] == _pos[1] ? 1.0 : 0.0 ;});
	S[0] = 100000.0;
	FullTensor Vt({N,N}, [](const std::vector<size_t>& _pos) {
		if(_pos[0] == 0) {return _pos[1] == 0 ? 1/std::sqrt(N) : 1/std::sqrt(_pos[1]*(_pos[1]+1)); }
		else if(_pos[1] == 0) { return -1.0/std::sqrt(N); }
		else if(_pos[0] == _pos[1]) { return double(_pos[0])/std::sqrt(_pos[1]*(_pos[1]+1)); }
		else if(_pos[0] < _pos[1]) { return -1.0/std::sqrt(_pos[1]*(_pos[1]+1)); }
		else { return 0.0; }});
// 	LOG(test, "U: " << std::endl << U.to_string());
// 	LOG(test, "S: " << std::endl << S.to_string());
// 	LOG(test, "Vt: " << std::endl << Vt.to_string());
	
	
/*	
	FullTensor UU;
	UU(i,k) = U(j,i) * U(j,k);
	LOG(test, "UU: " << std::endl << UU.to_string());
	
	UU(i,k) = U(i,j) * U(k,j);
	LOG(test, "UU: " << std::endl << UU.to_string());*/
	
	FullTensor X;
	
	X(i,l) = U(i,j)*S(j,k)*Vt(k,l);
	LOG(test, "X: " << std::endl << X[0] << " <=> " << X[X.size-1]);
// 	LOG(test, "X: " << std::endl << X.to_string());
	
	FullTensor M, P, L;
	(M(i,j), L(j,k), P(k,l)) = SVD(X(i,l), 1.0); 
	
// 	LOG(test, "M: " << std::endl << M.to_string());
// 	LOG(test, "L: " << std::endl << L.to_string());
// 	LOG(test, "P: " << std::endl << P.to_string());
	
	
	X(i,l) = M(i,j)* L(j,k)* P(k,l);
	
	LOG(test, "X: " << std::endl << X[0] << " <=> " << X[X.size-1]);
)

UNIT_TEST(Algorithm, largestEntry,
    //Random numbers
    std::mt19937_64 rnd;
    rnd.seed(73);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	std::uniform_int_distribution<size_t> dimDist(1,4);
	std::uniform_int_distribution<size_t> rankDist(1,6);
    
	const size_t D = 15;
	
	for(size_t k = 0; k < 250; ++k) {
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
