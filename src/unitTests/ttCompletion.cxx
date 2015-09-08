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


UNIT_TEST(Algorithm, adf_completion,
	const size_t D = 10;
	std::random_device rd;
	std::mt19937_64 rnd(rd());
	std::uniform_int_distribution<size_t> dist(0,3);
	std::normal_distribution<value_t> distF(0,1);
// 	TTTensor trueSolution(examples::peaking_diagonals(D,4,1.0));
	TTTensor trueSolution = TTTensor::random(std::vector<size_t>(D, 4), std::vector<size_t>(D-1, 5), rnd, distF);
	std::vector<SinglePointMeasurment> measurements;
	std::set<SinglePointMeasurment, SinglePointMeasurment::Comparator> measSet;
	
	for (size_t i=0; i<100*D*5*5*4; ++i) {
		std::vector<size_t> pos;
		for (size_t n=0; n<D; ++n) {
			pos.emplace_back(dist(rnd));
		}
		measSet.emplace(pos, 0.0);
	}
	
	measurements.insert(measurements.begin(), measSet.begin(), measSet.end());
	
	trueSolution.measure(measurements);
	
	bool test = true;
	for (const SinglePointMeasurment &m : measurements) {
		MTEST(misc::approx_equal(m.value, trueSolution[m.positions], 1e-12), m.value - trueSolution[m.positions]); 
	}
	TEST(test);
	
	TTTensor X = TTTensor::random(trueSolution.dimensions, trueSolution.ranks(), rnd, distF);
	X /= X.frob_norm()*1000.0;
	
	ADF(X, measurements);
	
	TEST(misc::approx_equal(frob_norm(X - trueSolution),0.0,1e-13));
)
