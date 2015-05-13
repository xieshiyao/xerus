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

#include "../../include/xerus.h"
#include "../../include/misc/test.h"
using namespace xerus;

UNIT_TEST(TT, TTTensor_Creation,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Index i,j,k;
    
    FullTensor A1 = FullTensor::construct_random({2}, rnd, dist);
    TTTensor TTA1(A1, 1e-14);
    FullTensor B1(1);
    B1(i) = TTA1(i);
    TEST(approx_equal(B1,A1, 1e-12));
    
    FullTensor A2 = FullTensor::construct_random({2,2}, rnd, dist);
    TTTensor TTA2(A2, 1e-14);
    FullTensor B2(2);
    B2(j,i) = TTA2(j,i);
    TEST(approx_equal(B2,A2, 1e-12));
    
    FullTensor A3 = FullTensor::construct_random({2,7}, rnd, dist);
    TTTensor TTA3(A3, 1e-14);
    FullTensor B3(2);
    B3(j,i) = TTA3(j,i);
    TEST(approx_equal(B3,A3, 1e-12));
    
    FullTensor A4 = FullTensor::construct_random({2,2,2,2,2,2,2,2}, rnd, dist);
    TTTensor TTA4(A4, 1e-14);
    FullTensor B4(8);
    B4(j,i^7) = TTA4(j,i^7);
    TEST(approx_equal(B4,A4, 1e-12));
    
    FullTensor A5 = FullTensor::construct_random({7,5,3,1,4,2,8,1}, rnd, dist);
    TTTensor TTA5(A5, 1e-14);
    FullTensor B5(8);
    B5(j,i^7) = TTA5(j,i&1);
    TEST(approx_equal(B5,A5, 1e-12));
)


UNIT_TEST(TT, TTOperator_Creation,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Index i,j,k;
    
    FullTensor A1 = FullTensor::construct_random({2,2}, rnd, dist);
    TTOperator TTA1(A1, 1e-14);
    FullTensor B1(2);
    B1(i^2) = TTA1(i^2);
    TEST(approx_equal(B1,A1, 1e-12));
    
    FullTensor A2 = FullTensor::construct_random({2,7}, rnd, dist);
    TTOperator TTA2(A2, 1e-14);
    FullTensor B2(2);
    B2(j,i) = TTA2(j,i);
    TEST(approx_equal(B2,A2, 1e-12));
    
    FullTensor A3 = FullTensor::construct_random({2,7,3,1}, rnd, dist);
    TTOperator TTA3(A3, 1e-14);
    FullTensor B3(4);
    B3(j,i^3) = TTA3(j,i^3);
    TEST(approx_equal(B3,A3, 1e-12));
    
    FullTensor A4 = FullTensor::construct_random({2,2,2,2,2,2,2,2}, rnd, dist);
    TTOperator TTA4(A4, 1e-14);
    FullTensor B4(8);
    B4(j,i^7) = TTA4(j,i&1);
    TEST(approx_equal(B4,A4, 1e-12));
    
    FullTensor A5 = FullTensor::construct_random({7,5,6,3,1,4,2,1}, rnd, dist); 
    TTOperator TTA5(A5, 1e-14);
    FullTensor B5(8);
    B5(j,i^7) = TTA5(j,i^7);
    TEST(approx_equal(B5,A5, 1e-12));
)


UNIT_TEST(TT, creation_with_epsilon,
	const value_t EPS = 0.01;
	//Random numbers
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A = FullTensor::construct_random({5,5,5,5}, rnd, dist);
	TTTensor ttA(A, EPS); 
	TTTensor ttB(A, 0);
	ttB.round(EPS);
	
	Index i;
	
	TEST(frob_norm(A(i&0)-FullTensor(ttA)(i&0)) < 0.1);
	TEST(ttA.ranks()[1] < 25);
	TEST(frob_norm(A(i&0)-FullTensor(ttB)(i&0)) < 0.1);
	TEST(ttB.ranks()[1] < 25);
)

UNIT_TEST(TT, creation_from_fullTensor_5x5x5x5,
	//Random numbers
	std::mt19937_64 rnd;
	rnd.seed(73);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A = FullTensor::construct_random({5,5,5,5}, rnd, dist);
	TTTensor ttA(A);
	FullTensor B(ttA);
	
	Index i;
	
	TEST(approx_equal(A,B, 1e-12));
	TEST(frob_norm(A(i&0)-B(i&0)) < 1e-13*5*5*5*5);
)


