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
#include "../../include/xerus/misc/missingFunctions.h"
using namespace xerus;


UNIT_TEST(Algorithm, adf_completion,
	const size_t D = 4;
	const size_t N = 12;
	const size_t R = 3; 
	const size_t CS = 10;
// 	std::random_device rd;
// 	std::mt19937_64 rnd(rd());
	std::mt19937_64 rnd;
	std::uniform_int_distribution<size_t> dist(0, N-1);
	std::uniform_real_distribution<value_t> distF(-0.5 ,0.5);
	TTTensor trueSolution = TTTensor::random(std::vector<size_t>(D, N), std::vector<size_t>(D-1, R), rnd, distF);
	std::vector<SinglePointMeasurment> measurements;
	std::set<SinglePointMeasurment, SinglePointMeasurment::Comparator> measSet;
	
	for(size_t d = 0; d < D; ++d) {
		for(size_t n = 0; n < N; ++n) {
			while(measSet.size() < d*N*CS*R*R + (n+1)*CS*R*R) {
				std::vector<size_t> pos;
				for (size_t i = 0; i < d; ++i) {
					pos.emplace_back(dist(rnd));
				}
				pos.emplace_back(n);
				for (size_t i = d+1; i < D; ++i) {
					pos.emplace_back(dist(rnd));
				}
				measSet.emplace(pos, 0.0);
			}
			
		}
	}

	measurements.insert(measurements.end(), measSet.begin(), measSet.end());
// 	LOG(bla, "Set size " << measurements.size() << " should be " << CS*D*N*R*R);
	
	trueSolution.measure(measurements);
	bool test = true;
	for (const SinglePointMeasurment &m : measurements) {
		test = test && misc::approx_equal(m.value, trueSolution[m.positions], 1e-10);
	}
	TEST(test);
	
	TTTensor X = TTTensor::random(trueSolution.dimensions, std::vector<size_t>(D-1, 1), rnd, distF);
	X /= X.frob_norm();
	
	const ADFVariant ADF20(10, EPSILON, true);  
	
	for(size_t r = 1; r < R; ++r) {
		ADF20(X, measurements);
		TTTensor rankInc = TTTensor::random(trueSolution.dimensions, std::vector<size_t>(D-1, 1), rnd, distF);
		X = X+rankInc;
	}
	
	ADF(X, measurements);
	
	MTEST(frob_norm(X - trueSolution)/frob_norm(trueSolution) < 1e-13, frob_norm(X - trueSolution)/frob_norm(trueSolution));
)
