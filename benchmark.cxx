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

#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>

#include <vector>
#include <set>

#include "include/xerus.h"

std::random_device rd;
std::mt19937_64 rnd(rd());
std::normal_distribution<double> normalDist(0,1);

using namespace xerus;

struct leastSquaresProblem {
	std::string name;
	std::vector<size_t> dimensions;
	std::vector<size_t> x_ranks;
	std::vector<size_t> b_ranks;
	leastSquaresProblem(const std::string &_name)
		: name(_name) {};
	
	virtual TTOperator get_a() const {
		std::vector<size_t> dim(dimensions);
		dim.insert(dim.end(), dimensions.begin(), dimensions.end());
		return TTOperator::identity(dim);
	}
	virtual TTTensor get_x() const {
		return TTTensor::random(dimensions, x_ranks, rnd, normalDist);
	};
	virtual TTTensor get_b() const {
		return TTTensor::random(dimensions, b_ranks, rnd, normalDist);
	};
};

namespace ls {
	struct approximation : public leastSquaresProblem {
		approximation(size_t _n, size_t _d, size_t _rankB, size_t _rankX)
			: leastSquaresProblem("approximation")
		{
			dimensions = std::vector<size_t>(_d, _n);
			x_ranks = std::vector<size_t>(_d-1, _rankX);
			b_ranks = std::vector<size_t>(_d-1, _rankB);
		};
	};
	
	struct random : public leastSquaresProblem {
		std::vector<size_t> a_ranks;
		
		random(size_t _n, size_t _d, size_t _rankA, size_t _rankB, size_t _rankX)
			: leastSquaresProblem("random")
		{
			dimensions = std::vector<size_t>(_d, _n);
			a_ranks = std::vector<size_t>(_d-1, _rankA);
			x_ranks = std::vector<size_t>(_d-1, _rankX);
			b_ranks = std::vector<size_t>(_d-1, _rankB);
		};
		
		TTOperator get_a() const override {
			std::vector<size_t> dim(dimensions);
			dim.insert(dim.end(), dimensions.begin(), dimensions.end());
			return TTOperator::random(dim, a_ranks, rnd, normalDist);
		}
	};
	
	
}



std::vector<leastSquaresProblem> leastSquaresProblems{
	ls::approximation(2, 10, 4, 2),
	ls::random(2, 10, 3, 3, 3)
};

int main() {
	std::string profileName;
#ifdef TEST_COVERAGE_
	static_assert(false, "test coverage checking nonsensical with benchmark run");
#endif
#ifdef LOW_OPTIMIZATION
	profileName += "lowOpt";
#elif defined(HIGH_OPTIMIZATION)
	profileName += "highOpt";
#elif defined(DANGEROUS_OPTIMIZATION)
	profileName += "dangerousOpt";
#elif defined(RIDICULOUS_OPTIMIZATION)
	profileName += "ridiculousOpt";
#else
	profileName += "noOpt";
#endif
	
#ifdef USE_LTO
	profileName += "_lto";
#endif
#ifdef DISABLE_RUNTIME_CHECKS_
	profileName += "_noChecks";
#endif
#ifdef REPLACE_ALLOCATOR
	profileName += "_replaceAlloc";
#endif
#ifdef PERFORMANCE_ANALYSIS
	profileName += "_perfAnalysis";
#endif
	
	LOG(benchmark, "running profile " << profileName);
	
	
	
	return 0;
}
