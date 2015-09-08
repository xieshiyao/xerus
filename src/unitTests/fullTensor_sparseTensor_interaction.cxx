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

UNIT_TEST(FullTensor_SparseTensor_Interaction, Assignment, 
    SparseTensor A({2,2,3,1,2});
    FullTensor B;
    FullTensor resF;
    SparseTensor resS;
    
    Index i,j,k,l,m;
    
    A[{0,0,0,0,0}]=1;
    A[{0,0,0,0,1}]=2;
    A[{0,0,1,0,0}]=3;
    A[{0,0,1,0,1}]=4;
    A[{0,0,2,0,0}]=5;
    A[{0,0,2,0,1}]=6;
    A[{0,1,0,0,0}]=7;
    A[{0,1,0,0,1}]=8;
    A[{0,1,1,0,0}]=9;
    A[{0,1,1,0,1}]=10;
    A[{0,1,2,0,0}]=11;
    A[{0,1,2,0,1}]=12;
    A[{1,0,0,0,0}]=13;
    A[{1,0,0,0,1}]=14;
    A[{1,0,1,0,0}]=15;
    A[{1,0,1,0,1}]=16;
    A[{1,0,2,0,0}]=17;
    A[{1,0,2,0,1}]=18;
    A[{1,1,0,0,0}]=19;
    A[{1,1,0,0,1}]=20;
    A[{1,1,1,0,0}]=21;
    A[{1,1,1,0,1}]=22;
    A[{1,1,2,0,0}]=23;
    A[{1,1,2,0,1}]=24;
    
    
    B(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(B.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    
    resF = B + A;
    TEST(resF.compare_to_data({1+1,2+2,3+3,4+4,5+5,6+6,13+7,14+8,15+9,16+10,17+11,18+12,7+13,8+14,9+15,10+16,11+17,12+18,19+19,20+20,21+21,22+22,23+23,24+24}));
    
    resF = B - A;
    TEST(resF.compare_to_data({1-1,2-2,3-3,4-4,5-5,6-6,13-7,14-8,15-9,16-10,17-11,18-12,7-13,8-14,9-15,10-16,11-17,12-18,19-19,20-20,21-21,22-22,23-23,24-24}));
    
)

UNIT_TEST(FullTensor_SparseTensor_Interaction, Product,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Index i,j,k,l,m,n,o,p,q;
    
    SparseTensor AS = SparseTensor::random({2,3,4,3,5}, 23, rnd, dist);
    SparseTensor BS = SparseTensor::random({6,3,4,2,3}, 23, rnd, dist);
    
    FullTensor AF(AS);
    FullTensor BF(BS);

    FullTensor resSF;
    FullTensor resFS;
    FullTensor check;
    
    check(i,j,k,m,n,o,p,q) = AF(i,j,k,l,m)*BF(n,l,o,p,q);
    resSF(i,j,k,m,n,o,p,q) = AS(i,j,k,l,m)*BF(n,l,o,p,q);
    TEST(approx_equal(check, resSF, 1e-14));
    resFS(i,j,k,m,n,o,p,q) = AF(i,j,k,l,m)*BS(n,l,o,p,q);
    TEST(approx_equal(check, resFS, 1e-14));
    
    check(i,j,m,n,p,q) = AF(i,j,k,l,m)*BF(n,l,k,p,q);
    resSF(i,j,m,n,p,q) = AS(i,j,k,l,m)*BF(n,l,k,p,q);
    TEST(approx_equal(check, resSF, 1e-14));
    resFS(i,j,m,n,p,q) = AF(i,j,k,l,m)*BS(n,l,k,p,q);
    TEST(approx_equal(check, resFS, 1e-14));
    
    check(i,j,m,n,p,q) = AF(q,j,k,l,n)*BF(m,l,k,i,p);
    resSF(i,j,m,n,p,q) = AS(q,j,k,l,n)*BF(m,l,k,i,p);
    TEST(approx_equal(check, resSF, 1e-14));
    resFS(i,j,m,n,p,q) = AF(q,j,k,l,n)*BS(m,l,k,i,p);
    TEST(approx_equal(check, resFS, 1e-14));
    
    check(i,m,n,q) = AF(q,j,k,l,n)*BF(m,l,k,i,j);
    resSF(i,m,n,q) = AS(q,j,k,l,n)*BF(m,l,k,i,j);
    TEST(approx_equal(check, resSF, 1e-14));
    resFS(i,m,n,q) = AF(q,j,k,l,n)*BS(m,l,k,i,j);
    TEST(approx_equal(check, resFS, 1e-14));
)

UNIT_TEST(FullTensor_SparseTensor_Interaction, Random,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);
    std::uniform_int_distribution<size_t> intDist (1, 10);

    Index i,j,k,l,m,n,o,p,q;
	
	std::vector<size_t> dimensions;
	for(size_t d = 0; d < 6; ++d) {
		std::vector<size_t> opDim(dimensions);
		opDim.insert(opDim.end(), dimensions.begin(), dimensions.end());
		std::uniform_int_distribution<size_t> numDist (0, misc::product(dimensions));
		
		SparseTensor AS = SparseTensor::random(dimensions, numDist(rnd), rnd, dist);
		SparseTensor BS = SparseTensor::random(dimensions, numDist(rnd), rnd, dist);
		SparseTensor CS = SparseTensor::random(dimensions, numDist(rnd), rnd, dist);
		SparseTensor I = Tensor::identity(opDim);
		
		FullTensor AF(AS);
		FullTensor BF(BS);
		FullTensor CF(CS);
		
		TEST(approx_equal(AS, AF));
		TEST(approx_equal(AF, AS));
		TEST(approx_equal(BS, BF));
		TEST(approx_equal(BF, BS));
		TEST(approx_equal(CS, CF));
		TEST(approx_equal(CF, CS));

		FullTensor resFF;
		FullTensor resFS;
		FullTensor resSF;
		SparseTensor resSS;
		
		// Simple Addition
		resFF = AF + BF;
		resFS = AF + BS;
		resSF = AS + BF;
		resSS = AS + BS;
		
		TEST(approx_equal(resFF, resFS));
		TEST(approx_equal(resFF, resSF));
		TEST(approx_equal(resFF, resSS));
		
		// Simple Subtraction
		resFF = AF - BF;
		resFS = AF - BS;
		resSF = AS - BF;
		resSS = AS - BS;
		
		TEST(approx_equal(resFF, resFS));
		TEST(approx_equal(resFF, resSF));
		TEST(approx_equal(resFF, resSS));
		
		// Addition + Subtraction
		resFF = AF - BF + CF;
		resFS = AF - BS + CF;
		resSF = AS - BF + CS;
		resSS = AS - BS + CS;
		
		TEST(approx_equal(resFF, resFS));
		TEST(approx_equal(resFF, resSF));
		TEST(approx_equal(resFF, resSS));
		
		// Addition + Subtraction + scale A
		resFF = 7.3*AF - 3.7*(BF + (-9.5)*CF);
		resFS = 7.3*AF - 3.7*(BS + (-9.5)*CF);
		resSF = 7.3*AS - 3.7*(BF + (-9.5)*CS);
		resSS = 7.3*AS - 3.7*(BS + (-9.5)*CS);
		
		TEST(approx_equal(resFF, resFS));
		TEST(approx_equal(resFF, resSF));
		TEST(approx_equal(resFF, resSS));
		
		// Addition + Subtraction + scale B
		resFF = 7300.9*AF - 3.7*(0*BF + 1*CF);
		resFS = 7300.9*AF - 3.7*(0*BS + 1*CF);
		resSF = 7300.9*AS - 3.7*(0*BF + 1*CS);
		resSS = 7300.9*AS - 3.7*(0*BS + 1*CS);
		
		TEST(approx_equal(resFF, resFS));
		TEST(approx_equal(resFF, resSF));
		TEST(approx_equal(resFF, resSS));
		
		
		// Addition + Subtraction + scale C
		resFF = 0.9*AF - 3.7*(7*BF + 2.1*CF) + 2037*(0*AF + 0*BF);
		resFS = 0.9*AF - 3.7*(7*BS + 2.1*CF) + 2037*(0*AS + 0*BF);
		resSF = 0.9*AS - 3.7*(7*BF + 2.1*CS) + 2037*(0*AF + 0*BS);
		resSS = 0.9*AS - 3.7*(7*BS + 2.1*CS) + 2037*(0*AS + 0*BS);
		
		TEST(approx_equal(resFF, resFS));
		TEST(approx_equal(resFF, resSF));
		TEST(approx_equal(resFF, resSS));
		
		dimensions.push_back(intDist(rnd));
	}
)
