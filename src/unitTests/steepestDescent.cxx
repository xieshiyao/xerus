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
#include<fstream>

#include "../../include/xerus/misc/test.h"
using namespace xerus;


UNIT_TEST(SteepestDescent, random_operator,
	std::random_device rd;
	std::mt19937_64 rnd(rd());
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	const size_t d = 7;
	const std::vector<size_t> stateDims(d, 5);
	const std::vector<size_t> operatorDims(2*d, 5);
	
	TTOperator A = TTOperator::construct_random(operatorDims, std::vector<size_t>(d-1, 4), rnd, dist);
	TTTensor b = TTTensor::construct_random(stateDims, std::vector<size_t>(d-1,4), rnd, dist);
	Index i,j,k;
	b(i&0) = A(i/2,j/2)*b(j&0);
	TTTensor initX = TTTensor::construct_random(stateDims, std::vector<size_t>(d-1,4), rnd, dist);
	TTTensor x(initX);
	PerformanceData perfdata;
	perfdata << "testcase with random operator\n";
	
	SteepestDescentVariant::HOSVDRetraction svdRet(2ul);
	SteepestDescentVariant hosvdGrad(svdRet);
	hosvdGrad.printProgress = true;
	hosvdGrad(A,x,b,1000ul,perfdata);
	std::cout << "HOSVD: " << perfdata.data.size() << " " << perfdata.data.back().elapsedTime << std::endl;
	perfdata.dumpToFile("hosvd.dat");
	
	
	SteepestDescentVariant alsGrad(SteepestDescentVariant::ALSRetraction);
	alsGrad.printProgress = true;
	perfdata.reset();
	x = initX;
	alsGrad(A,x,b,1000ul, perfdata);
	std::cout << "ALS: " << perfdata.data.size()  << " " << perfdata.data.back().elapsedTime << std::endl;
	perfdata.dumpToFile("als.dat");
	
	SteepestDescentVariant submanGrad(SteepestDescentVariant::SubmanifoldRetraction);
	submanGrad.printProgress = true;
	perfdata.reset();
	x = initX;
	submanGrad(A,x,b,1000ul, perfdata);
	std::cout << "Submanifold: " << perfdata.data.size()  << " " << perfdata.data.back().elapsedTime << std::endl;
	perfdata.dumpToFile("submani.dat");
)
