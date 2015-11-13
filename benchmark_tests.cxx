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

#include <xerus.h>
#include <timeMeasure.h>

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
	
	virtual TTOperator get_a() const = 0;
	virtual TTTensor get_x() const {
		return TTTensor::random(dimensions, x_ranks, normalDist, rnd);
	};
	virtual TTTensor get_b() const {
		return TTTensor::random(dimensions, b_ranks, normalDist, rnd);
	};
};

namespace ls {
	struct ls_approximation : public leastSquaresProblem {
		ls_approximation(size_t _n, size_t _d, size_t _rankB, size_t _rankX)
			: leastSquaresProblem("approximation")
		{
			dimensions = std::vector<size_t>(_d, _n);
			x_ranks = std::vector<size_t>(_d-1, _rankX);
			b_ranks = std::vector<size_t>(_d-1, _rankB);
		};
		
		TTOperator get_a() const override {
			std::vector<size_t> dim(dimensions);
			dim.insert(dim.end(), dimensions.begin(), dimensions.end());
			return TTOperator::identity(dim);
		}
	};
}



std::vector<leastSquaresProblem> leastSquaresProblems{
	leastSquaresProblem("random_2tothe10_rank3", ),
};

int main() {
    using namespace xerus;
    LOG(benchmark, "Initializing...");
    
    TimeMeasure timer;
    
    LOG(benchmark, "Finished Initializing. Starting first Round...");
    
    for(size_t i = 1; i <= 1000; ++i) {
        for(size_t n = 1; n < 1000; ++n) {
            LOG(warning, "Starting " << n);
            
            std::uniform_int_distribution<size_t> idDist(0, 2*n);
            
            // Create test Indices
            std::vector<Index> testVector;
            while(testVector.size() < 1000) {
                testVector.emplace_back(idDist(rnd), 1);
            }
            
            
            // Create Vector
            std::vector<Index> idxVec;
            while(idxVec.size() < n) {
                size_t id = idDist(rnd);
                if(!contains(idxVec, Index(id, 1))) {
                    idxVec.emplace_back(id, 1);
                }
            }           
            
            // Test Vector
            timer.step();
            size_t vecHits = 0;
            for(const Index& idx : testVector) {
                if(contains(idxVec, idx)) { vecHits++; }
            }
            add_call("Vector", {n}, timer.get());
            
            // Create Set
            std::set<Index> idxSet;
            for(const Index& idx : idxVec) {
                idxSet.insert(idx);
            }
            REQUIRE(idxSet.size() == n, " Fail " << idxVec.size() << " vs " << idxSet.size() << " not " << n);
            
            timer.step();
            size_t setHits = 0;
            for(const Index& idx : testVector) {
                if(contains(idxSet, idx)) { setHits++; }
            }
            add_call("Set", {n}, timer.get());
            
            REQUIRE(vecHits == setHits, "oO");
        }
        
        LOG(warning, "Finished one Round, writing to files...");
        for(const std::pair<std::string, std::map<std::vector<size_t>, size_t>>& call : results) {
            std::fstream file("BenchmarkTests/"+call.first, std::fstream::out);
            for(const std::pair<std::vector<size_t>, size_t>& callEntry : call.second) {
                for(size_t par : callEntry.first) {
                    file << par << " ";
                }
                file << 10*callEntry.second/i << std::endl;
            }
            file.close();
        }
        LOG(warning, "Finished writing to files. Starting next round...");
    }
    
    return 0;
}
