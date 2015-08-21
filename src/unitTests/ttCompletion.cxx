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
// 	std::random_device rd;
	std::mt19937_64 rnd(0xDEAFBEEF);
	std::uniform_int_distribution<size_t> dist(0,3);
	std::normal_distribution<value_t> distF(0,1);
	TTTensor trueSolution(examples::peaking_diagonals(10,4,1.0));
	std::vector<SinglePointMeasurment> measurements;
	std::set<SinglePointMeasurment, SinglePointMeasurment::Comparator> measSet;
	
	for (size_t i=0; i<10*10*5*5*4; ++i) {
		std::vector<size_t> pos;
		for (size_t n=0; n<10; ++n) {
			pos.emplace_back(dist(rnd));
		}
		measSet.emplace(pos, 0.0);
	}
	
	measurements.insert(measurements.begin(), measSet.begin(), measSet.end());
	
	trueSolution.measure(measurements);
	
	for (const SinglePointMeasurment &m : measurements) {
		TEST(misc::approx_equal(m.value, trueSolution[m.positions], 1e-13));
	}
	
	TTTensor X = TTTensor::construct_random(trueSolution.dimensions, trueSolution.ranks(), rnd, distF);
	
	ADF(X, measurements);
	
	TEST(misc::approx_equal(frob_norm(X - trueSolution),0.0,1e-13));
)
