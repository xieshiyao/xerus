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

UNIT_TEST(TT, sum,
	//Random numbers
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
    std::normal_distribution<value_t> dist (0.0, 1.0);
    std::uniform_int_distribution<size_t> intDist (1, 10);
    
    Index i;
    
    std::vector<size_t> dimensions;
    size_t d = 4;
        dimensions.push_back(intDist(rnd));
		dimensions.push_back(intDist(rnd));
		dimensions.push_back(intDist(rnd));
		dimensions.push_back(intDist(rnd));
        FullTensor A = FullTensor::construct_random(dimensions, rnd, dist);
        FullTensor B = FullTensor::construct_random(dimensions, rnd, dist);
        FullTensor C(d);
        TTTensor ttA(A); 
        TTTensor ttB(B); 
        TTTensor ttC(d);
		TTOperator toA(A); 
        TTOperator toB(B); 
        TTOperator toC(d);
        
        C(i&0) = A(i&0) + B(i&0);
        ttC(i&0) = ttA(i&0) + ttB(i&0);
		toC(i&0) = toA(i&0) + toB(i&0);
        TEST(frob_norm(FullTensor(ttC)(i&0) - C(i&0)) < 3.1*1e-13);
		TEST(frob_norm(FullTensor(toC)(i&0) - C(i&0)) < 3.1*1e-13);
	
)

UNIT_TEST(TT, difference,
	//Random numbers
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor B = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor C(4);
	TTTensor ttA(A); 
	TTTensor ttB(B); 
	TTTensor ttC(4); 
	
	Index i;
	C(i&0) = A(i&0) - B(i&0);
	ttC(i&0) = ttA(i&0) - ttB(i&0);
	
	double fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 6*1e-13);
)

UNIT_TEST(TT, real_difference,
    //Random numbers
    std::mt19937_64 rnd;
    rnd.seed(0X5EED);
    std::normal_distribution<value_t> dist (0.0, 1.0);
    
    TTTensor ttA = TTTensor::construct_random({10,10,10,10,10}, {4,4,4,4}, rnd, dist);
    TTTensor ttB = TTTensor::construct_random({10,10,10,10,10}, {4,4,4,4}, rnd, dist); 
    TTTensor ttC(5); 
    
    Index i;
    ttC(i&0) = ttA(i&0) - ttA(i&0);
    LOG(unit_tests, "Frob norm 1 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 1e-11);
    
    ttC(i&0) = ttB(i&0) - ttB(i&0);
    LOG(unit_tests, "Frob norm 2 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 1e-11);
    
    ttC(i&0) = (ttA(i&0) + ttB(i&0)) - (ttA(i&0) + ttB(i&0));
    LOG(unit_tests, "Frob norm 3 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 1e-11);
    
    ttC(i&0) = (ttA(i&0) + ttB(i&0)) - (ttB(i&0) + ttA(i&0));
    LOG(unit_tests, "Frob norm 4 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 1e-11);
    
    ttC(i&0) = (73*ttA(i&0) + ttB(i&0)) - (ttB(i&0) + 73*ttA(i&0));
    LOG(unit_tests, "Frob norm 5 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 5e-10);
)

UNIT_TEST(TT, difference_of_TTStacks,
    //Random numbers
    std::mt19937_64 rnd;
    rnd.seed(0X5EED);
    std::normal_distribution<value_t> dist (0.0, 1.0);
    
    TTOperator ttO = TTOperator::construct_random({10,10,10,10,10,10,10,10,10,10}, {4,4,4,4}, rnd, dist);
    TTTensor ttA = TTTensor::construct_random({10,10,10,10,10}, {4,4,4,4}, rnd, dist);
    TTTensor ttB = TTTensor::construct_random({10,10,10,10,10}, {4,4,4,4}, rnd, dist); 
    TTTensor ttC;
    
    Index i,j,k;
    ttC(i&0) = ttO(i/2, j/2)*ttA(j&0) - ttO(i/2, j/2)*ttA(j&0);
    LOG(unit_tests, "Frob norm 1 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 1e-7);
    
    ttC(i&0) = ttO(i/2, j/2)*ttB(j&0) - ttO(i/2, j/2)*ttB(j&0);
    LOG(unit_tests, "Frob norm 2 " << frob_norm(ttC(i&0)));
    TEST(frob_norm(ttC(i&0)) < 1e-7);
)

UNIT_TEST(TT, special_sum_diff,
	//Random numbers
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A({10,10,10,10}); // NOTE that this is the 0 tensor
	FullTensor B = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor C(4);
	TTTensor ttA(A); 
	TTTensor ttB(B); 
	TTTensor ttC(4); 
	
	
	Index i;
	
	C(i&0) = A(i&0) + B(i&0);
	ttC(i&0) = ttA(i&0) + ttB(i&0);
	TEST(frob_norm(FullTensor(ttC)(i&0) - C(i&0)) < 5*1e-13);
	TEST(frob_norm(FullTensor(ttC)(i&0) - FullTensor(ttB)(i&0)) < 3.1*1e-13);
	
	C(i&0) = B(i&0) + A(i&0);
	ttC(i&0) = ttB(i&0) + ttA(i&0);
	TEST(frob_norm(FullTensor(ttC)(i&0) - C(i&0)) < 5*1e-13);
	TEST(frob_norm(FullTensor(ttC)(i&0) - FullTensor(ttB)(i&0)) < 3.1*1e-13);
	
	C(i&0) = A(i&0) - B(i&0);
	ttC(i&0) = ttA(i&0) - ttB(i&0);
	TEST(frob_norm(FullTensor(ttC)(i&0) - C(i&0)) < 5*1e-13);
	TEST(frob_norm(FullTensor(ttC)(i&0) + FullTensor(ttB)(i&0)) < 3.1*1e-13);
	
	C(i&0) = B(i&0) - A(i&0);
	ttC(i&0) = ttB(i&0) - ttA(i&0);
	TEST(frob_norm(FullTensor(ttC)(i&0) - C(i&0)) < 5*1e-13);
	TEST(frob_norm(FullTensor(ttC)(i&0) - FullTensor(ttB)(i&0)) < 3.1*1e-13);
	
	FullTensor X({10});
	FullTensor Y = FullTensor::construct_random({10}, rnd, dist);
	FullTensor Z(1);
	TTTensor ttX(X);
	TTTensor ttY(Y);
	TTTensor ttZ(1);
	
	Z(i&0) = X(i&0) + Y(i&0);
	ttZ(i&0) = ttX(i&0) + ttY(i&0);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - Z(i&0)) < 3.1*1e-13);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - FullTensor(ttY)(i&0)) < 3.1*1e-13);
	
	Z(i&0) = Y(i&0) + X(i&0);
	ttZ(i&0) = ttY(i&0) + ttX(i&0);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - Z(i&0)) < 3.1*1e-13);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - FullTensor(ttY)(i&0)) < 3.1*1e-13);
	
	Z(i&0) = X(i&0) - Y(i&0);
	ttZ(i&0) = ttX(i&0) - ttY(i&0);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - Z(i&0)) < 3.1*1e-13);
	TEST(frob_norm(FullTensor(ttZ)(i&0) + FullTensor(ttY)(i&0)) < 3.1*1e-13);
	
	Z(i&0) = Y(i&0) - X(i&0);
	ttZ(i&0) = ttY(i&0) - ttX(i&0);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - Z(i&0)) < 3.1*1e-13);
	TEST(frob_norm(FullTensor(ttZ)(i&0) - FullTensor(ttY)(i&0)) < 3.1*1e-13);
)

UNIT_TEST(TT, product,
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	Index i,j,k,l;
	
	TTOperator ttA = TTOperator::construct_random({10,10,10,10}, {1}, rnd, dist);
	TTOperator ttB = TTOperator::construct_random({10,10,10,10}, {1}, rnd, dist);
	TTTensor ttD = TTTensor::construct_random({10,10}, {2}, rnd, dist);
	FullTensor A(ttA);
	FullTensor B(ttB);
	FullTensor D(ttD);
	
	FullTensor C(4);
	TTOperator ttC(4);
	
	D(i^2) = A(i^2,j^2) * D(j^2);
	ttD(i^2) = ttA(i^2,j^2) * ttD(j^2);
	double fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 10*10*10*10*1e-15);
	
	C(i^2,k^2) = A(i^2,j^2) * B(j^2,k^2);
	ttC(i^2,k^2) = ttA(i^2,j^2) * ttB(j^2,k^2);
	TEST(ttC.nodes.size() == 4);
	fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 10*10*10*10*1e-15);
	
	C(i/2,k/2) = A(j/2,i/2) * B(j/2,k/2);
	ttC(i^2,k/2) = ttA(j^2,i/2) * ttB(j^2,k/2);
	fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 10*10*10*10*1e-15);
	
	ttC(i^2,k/2) = ttB(j/2,k/2) * ttA(j^2,i^2);
	fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 10*10*10*10*1e-15);
	
	C(i^2,k^2) = A(i^2,j^2) * B(k^2,j^2);
	ttC(i^2,k^2) = ttA(i^2,j^2) * ttB(k^2,j^2);
	fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 10*10*10*10*1e-15);
	
	C(i^2,k^2) = A(j^2,i^2) * B(k^2,j^2);
	ttC(i^2,k^2) = ttA(j^2,i^2) * ttB(k^2,j^2);
	fnorm = frob_norm(FullTensor(ttC)(i&0) - C(i&0));
	LOG(unit_tests, "frob_norm " << fnorm);
	TEST(fnorm < 10*10*10*10*1e-15);
	
	FAILTEST(ttC(i^2,k^2) = ttA(j^2,i) * ttB(k^2,j^2));
	FAILTEST(ttC(i^2,k^2) = ttA(j^2,i^2) * ttB(k^2,k^2));
)

UNIT_TEST(TT, identities,
	FullTensor I = FullTensor({2,2,2,2},[](const std::vector<size_t>& _idx)->value_t{
		if ((_idx[0] == _idx[2] && _idx[1] == _idx[3]) || (_idx[0] == 0 && _idx[1]==1 &&_idx[2] ==0 && _idx[3]==0)) {
			return 1;
		} else {
			return 0;
		}
	}
	);
	TTOperator ttI(I);
	TTOperator ttC(4);
	FullTensor C(4);
	Index i,j,k;
	
	ttC(i^2,k^2) = ttI(i^2,j^2) * ttI(j^2,k^2);
	C(i^2,k^2) = I(i^2,j^2) * I(j^2,k^2);
	LOG(unit_test, frob_norm(C(i&0) - FullTensor(ttC)(i&0)));
	TEST(approx_equal(C, FullTensor(ttC), 1e-15));
	
	ttC(k^2,i^2) = ttI(i^2,j^2) * ttI(j^2,k^2);
	C(k^2,i^2) = I(i^2,j^2) * I(j^2,k^2);
	LOG(unit_test, frob_norm(C(i&0) - FullTensor(ttC)(i&0)));
	TEST(approx_equal(C, FullTensor(ttC), 1e-15));
	
	ttC(i^2,k^2) = ttI(i^2,j^2) * ttI(k^2,j^2);
	C(i^2,k^2) = I(i^2,j^2) * I(k^2,j^2);
	LOG(unit_test, frob_norm(C(i&0) - FullTensor(ttC)(i&0)));
	TEST(approx_equal(C, FullTensor(ttC), 1e-15));
	
	ttC(k^2,i^2) = ttI(i^2,j^2) * ttI(k^2,j^2);
	C(k^2,i^2) = I(i^2,j^2) * I(k^2,j^2);
	LOG(unit_test, frob_norm(C(i&0) - FullTensor(ttC)(i&0)));
	TEST(approx_equal(C, FullTensor(ttC), 1e-15));
)

UNIT_TEST(TT, transpose,
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	FullTensor A = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor B(4);
	TTOperator ttA(A); 
	Index i,j;
	
	B(i^2,j^2) = A(j^2,i^2);
	ttA.transpose();
	LOG(unit_test, frob_norm(B(i&0) - FullTensor(ttA)(i&0)));
	TEST(approx_equal(B, FullTensor(ttA), 1e-14));
)

UNIT_TEST(TT, ax_b,
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	TTTensor X = TTTensor::construct_random({10,10,10}, {2,2}, rnd, dist);
	TTTensor B = TTTensor::construct_random({10,10,10}, {2,2}, rnd, dist);
    
	FullTensor I({10,10,10,10,10,10}, [](const std::vector<size_t> &_idx) {
		if (_idx[0]==_idx[3] && _idx[1] == _idx[4] && _idx[2] == _idx[5]) {
			return 1.0;
		} else {
			return 0.0;
		}
	});
    
	TTOperator A(I);
	TTTensor T(3);
	TTTensor S(3);
    
	Index i,j,k;
	
	T(i^3) = A(i^3, j^3) * X(j^3);
	T(i^3) = T(i^3) - B(i^3);
	S(i^3) = A(i^3, j^3) * X(j^3) - B(i^3);
	LOG(unit_test, frob_norm(T(i^3)-S(i^3)));
	TEST(frob_norm(T(i^3)-S(i^3)) < 1e-7);
	
	FullTensor fA(A);
	FullTensor fX(X);
	FullTensor fT(3);
	fT(i^3) = fA(i^3, j^3) * fX(j^3);
	TEST(frob_norm(fT(i^3) - fX(i^3))<1e-7);
	
	T(i^3) = A(i^3, j^3) * X(j^3);
    TEST(frob_norm(FullTensor(T) - fT) < 1e-7);
    TEST(frob_norm(FullTensor(T)(i^3) - fT(i^3)) < 1e-7);
	
	T(i^3) = A(i^3, j^3) * X(j^3);
    TEST(frob_norm(A(i^3, j^3) * X(j^3) - T(i^3)) < 1e-7);
    
	TEST(frob_norm(T(i^3) - X(i^3))<1e-7);
	T(i^3) = T(i^3) - X(i^3);
	S(i^3) = A(i^3, j^3) * X(j^3) - X(i^3);
	LOG(unit_test, frob_norm(T(i^3)-S(i^3)));
	TEST(frob_norm(T(i^3)-S(i^3)) < 1e-7);
	TEST(frob_norm(S(i^3)) < 1e-7);
	
	T(i^3) = A(j^3, i^3) * X(j^3);
	TEST(frob_norm(T(i^3) - X(i^3))<1e-7);
	T(i^3) = T(i^3) - X(i^3);
	S(i^3) = A(j^3, i^3) * X(j^3) - X(i^3);
	LOG(unit_test, frob_norm(T(i^3)-S(i^3)));
	TEST(frob_norm(T(i^3)-S(i^3)) < 1e-7);
	TEST(frob_norm(S(i^3)) < 1e-7);
	
	T(i^3) = A(j^3, i^3) * B(j^3);
	T(i^3) = T(i^3) - B(i^3);
	S(i^3) = A(j^3, i^3) * B(j^3) - B(i^3);
	LOG(unit_test, frob_norm(T(i^3)-S(i^3)));
	TEST(frob_norm(T(i^3)-S(i^3)) < 1e-7);
	TEST(frob_norm(S(i^3)) < 1e-7);
)

UNIT_TEST(TT, operator_times_tensor,
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor B = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor C = FullTensor::construct_random({10,10}, rnd, dist);
	FullTensor D(2);
	FullTensor Do(4);
	TTOperator ttA(A); 
	ttA.round(2); A = FullTensor(ttA);
	TTOperator ttB(B); 
	ttB.round(2); B = FullTensor(ttB);
	TTTensor ttC(C); 
	TTTensor ttD(C);
	TTOperator ttDo(A);
	
	Index i,j,k,l,m;
	ttD(i^2) = ttA(i^2,j^2) * ttB(j^2,k^2) * ttC(k^2);
	D(i^2) = A(i^2,j^2) * B(j^2,k^2) * C(k^2);
	LOG(unit_test, frob_norm(D(i&0) - FullTensor(ttD)(i&0)));
	TEST(approx_equal(D, FullTensor(ttD), 1e-15));
	
	ttD(i^2) = ttA(i^2,j^2) * ttB(k^2,j^2) * ttC(k^2);
	D(i^2) = A(i^2,j^2) * B(k^2,j^2) * C(k^2);
	LOG(unit_test, frob_norm(D(i&0) - FullTensor(ttD)(i&0)));
	TEST(approx_equal(D, FullTensor(ttD), 1e-15));
	
	ttDo(i^2,k^2) = ttA(i^2,j^2) * ttB(j^2,k^2);
	Do(i^2,k^2) = A(i^2,j^2) * B(j^2,k^2);
	LOG(unit_test, frob_norm(Do(i&0) - FullTensor(ttDo)(i&0)));
	TEST(approx_equal(Do, FullTensor(ttDo), 1e-15));
	
	ttDo(i^2,k^2) = ttA(i^2,j^2) * ttA(j^2,k^2);
	Do(i^2,k^2) = A(i^2,j^2) * A(j^2,k^2);
	LOG(unit_test, frob_norm(Do(i&0) - FullTensor(ttDo)(i&0)));
	TEST(approx_equal(Do, FullTensor(ttDo), 1e-15));
	
	ttDo(i^2,l^2) = ttA(i^2,j^2) * ttB(j^2,k^2) * ttA(l^2,k^2);
	Do(i^2,l^2) = A(i^2,j^2) * B(j^2,k^2) * A(l^2,k^2);
	LOG(unit_test, frob_norm(Do(i&0) - FullTensor(ttDo)(i&0)));
	TEST(approx_equal(Do, FullTensor(ttDo), 1e-15));
	
	ttDo(i^2,l^2) = ttA(i^2,j^2) * ttB(j^2,k^2) * ttB(l^2,k^2);
	Do(i^2,l^2) = A(i^2,j^2) * B(j^2,k^2) * B(l^2,k^2);
	LOG(unit_test, frob_norm(Do(i&0) - FullTensor(ttDo)(i&0)));
	TEST(approx_equal(Do, FullTensor(ttDo), 1e-15));
	
	ttDo(i^2,m^2) = ttA(i^2,j^2) * ttB(j^2,k^2) * ttB(l^2,k^2) * ttA(l^2,m^2);
	Do(i^2,m^2) = A(i^2,j^2) * B(j^2,k^2) * B(l^2,k^2) * A(l^2,m^2);
	LOG(unit_test, frob_norm(Do(i&0) - FullTensor(ttDo)(i&0)));
	TEST(approx_equal(Do, FullTensor(ttDo), 1e-15));
)

UNIT_TEST(TT, full_contraction,
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	FullTensor B = FullTensor::construct_random({10,10,10,10}, rnd, dist);
	TTTensor ttA(A); 
	TTTensor ttB(B); 
	TTOperator toA(A); 
	TTOperator toB(B); 
	
	Index i;
	TEST(misc::approx_equal(frob_norm(A(i&0)), frob_norm(ttA(i&0)), 2e-13));
	TEST(misc::approx_equal(frob_norm(A(i&0)), frob_norm(toA(i&0)), 1.6e-13));
	TEST(misc::approx_equal(frob_norm(B(i&0)), frob_norm(ttB(i&0)), 1e-13));
	TEST(misc::approx_equal(frob_norm(B(i&0)), frob_norm(toB(i&0)), 1e-13));
	TEST(misc::approx_equal(frob_norm(A(i&0)-B(i&0)), frob_norm(ttA(i&0)-ttB(i&0)), 1e-12));
	TEST(misc::approx_equal(frob_norm(A(i&0)-B(i&0)), frob_norm(toA(i&0)-toB(i&0)), 1e-12));
	FullTensor C(0);
	C() = A(i/1)*B(i&0);
	TTTensor ttC(0);
	ttC() = ttA(i&0)*ttB(i&0);
	TTOperator toC(0);
	toC() = ttA(i&0)*ttB(i&0);
	TEST(misc::approx_equal(C[{}], ttC[{}], 1e-12));
	TEST(misc::approx_equal(C[{}], toC[{}], 1e-12));
)

UNIT_TEST(TT, disjoint_product,
	std::mt19937_64 rnd;
	rnd.seed(0X5EED);
	std::normal_distribution<value_t> dist (0.0, 1.0);
	
	FullTensor A = FullTensor::construct_random({10,10}, rnd, dist);
	FullTensor B = FullTensor::construct_random({10,10}, rnd, dist);
	FullTensor C(4);
	TTTensor ttA(A); 
	TTTensor ttB(B);
	TTTensor ttC(4);
	
	Index i,j;
	
	ttC = TTTensor::dyadic_product(ttA,ttB);
	C(i^2,j^2) = A(i&0)*B(j&0);
	
	LOG(unit_test, frob_norm(C(i&0) - FullTensor(ttC)(i&0)));
	TEST(approx_equal(C, FullTensor(ttC), 1e-13));
)


