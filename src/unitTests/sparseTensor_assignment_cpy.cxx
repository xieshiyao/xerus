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

UNIT_TEST(SparseTensor, Assignment_Trivia, 
    SparseTensor A({2,2,3,1,2});
    SparseTensor res({2,2,3,1,2});
    SparseTensor res2({2,3,2,1,2});
    SparseTensor res3({2,3,1,2,2});

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
    
    res(i,j,k,l,m) = A(i,j,k,l,m);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res(i&0) = A(i&0);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res(i&0) = A(i^5);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res(i^5) = A(i&0);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res(i^3, j&3) = A(i&2,j^2);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    
    res(j,i,k,l,m) = A(j,i,k,l,m);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res(i,j,k,l,m) = A(j,i,k,l,m);
    TEST(res.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    res(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(res.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    
    res2(i,k,j,l,m) = A(i,j,k,l,m);
    TEST(res2.compare_to_data({1,2,7,8,3,4,9,10,5,6,11,12,13,14,19,20,15,16,21,22,17,18,23,24}));
    res2(i,j,k,l,m) = A(i,k,j,l,m);
    TEST(res2.compare_to_data({1,2,7,8,3,4,9,10,5,6,11,12,13,14,19,20,15,16,21,22,17,18,23,24}));
    
    res(i,j,k,l,m) = A(i,j,l,k,m);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res(i,j,l,k,m) = A(i,j,k,l,m);
    TEST(res.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
)


UNIT_TEST(SparseTensor, Assignment_Const, 
    SparseTensor A({2,2,3,1,2});
    SparseTensor res({2,2,3,1,2});

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
    
    res(i,j,k,l,m) = A(i,j,k,l,m);
    const SparseTensor resC1(res);
    TEST(resC1.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    TEST(approx_equal(resC1[{0,0,0,0,0}], 1.0));
    TEST(approx_equal(resC1[{0,1,2,0,0}], 11.0));
    TEST(approx_equal(resC1[{1,1,2,0,0}], 23.0));
    TEST(approx_equal(resC1.at({1,0,0,0,1}), 14.0));
    TEST(approx_equal(resC1.at({1,0,1,0,1}), 16.0));
    
    res(j,i,k,l,m) = A(i,j,k,l,m);
    const SparseTensor resC2(res);
    TEST(resC2.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    TEST(approx_equal(resC2[{0,0,0,0,0}], 1.0));
    TEST(approx_equal(resC2[{0,1,2,0,0}], 17.0));
    TEST(approx_equal(resC2[{1,1,2,0,0}], 23.0));
    TEST(approx_equal(resC2.at({1,0,0,0,1}), 8.0));
    TEST(approx_equal(resC2.at({1,0,1,0,1}), 10.0));
)

UNIT_TEST(SparseTensor, Assignment_Overwriting_Dimensions,
    SparseTensor A({2,2,3,1,2});
    SparseTensor res1({2,2,3,1,2});
    SparseTensor res2({1,3,5,1,7});
    SparseTensor res3({2,3,2,2,3});
    SparseTensor res4({13,9,2,5,3});
    SparseTensor res5({1,1,1,1,1});
    SparseTensor res6;
    SparseTensor res7({1,2,3,4,5});
    SparseTensor res8({5,4,3,2,1});
    SparseTensor res9({3,2,3,2,3});
    SparseTensor res10({4,5,3,1,3});
    SparseTensor res11;
    SparseTensor res12({1,5,3,1,3});

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
    
    res1(i,j,k,l,m) = A(i,j,k,l,m);
    TEST(res1.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res2(i&0) = A(i&0);
    TEST(res2.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res3(i&0) = A(i^5);
    TEST(res3.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res4(i^5) = A(i&0);
    TEST(res4.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res5(i^3, j&3) = A(i&2,j^2);
    TEST(res5.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    
    res6(j,i,k,l,m) = A(j,i,k,l,m);
    TEST(res6.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res7(i,j,k,l,m) = A(j,i,k,l,m);
    TEST(res7.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    res8(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(res8.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    
    res9(i,k,j,l,m) = A(i,j,k,l,m);
    TEST(res9.compare_to_data({1,2,7,8,3,4,9,10,5,6,11,12,13,14,19,20,15,16,21,22,17,18,23,24}));
    res10(i,j,k,l,m) = A(i,k,j,l,m);
    TEST(res10.compare_to_data({1,2,7,8,3,4,9,10,5,6,11,12,13,14,19,20,15,16,21,22,17,18,23,24}));
    
    res11(i,j,k,l,m) = A(i,j,l,k,m);
    TEST(res11.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
    res12(i,j,l,k,m) = A(i,j,k,l,m);
    TEST(res12.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}));
)

UNIT_TEST(SparseTensor, Assignment_LHS_Equals_RHS, 
    SparseTensor B({2,2});
    SparseTensor C({2,2});

    Index i, J, K;
    
    B[{0,0}]=1;
    B[{0,1}]=2;
    B[{1,0}]=3;
    B[{1,1}]=4;
    
    C[{0,0}]=5;
    C[{0,1}]=6;
    C[{1,0}]=7;
    C[{1,1}]=8;
    
    B(i,J) = B(i,J);
    TEST(B.compare_to_data({1,2,3,4}));
    B(i,J) = B(J,i);
    TEST(B.compare_to_data({1,3,2,4}));
)

UNIT_TEST(SparseTensor, Assignment_Fixed_Indices,
    SparseTensor A({2,2,3,1,2});
    SparseTensor res1({2,3,1,2});
    SparseTensor res2;
    SparseTensor res3;
    
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
    
    res1(i,j,k,l) = A(0,i,j,k,l);
    TEST(res1.compare_to_data({1,2,3,4,5,6,7,8,9,10,11,12}));
    res1(i,j,k,l) = A(1,i,j,k,l);
    TEST(res1.compare_to_data({13,14,15,16,17,18,19,20,21,22,23,24}));
    
    res1(i,j,k,l) = A(i,0,j,k,l);
    TEST(res1.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18}));
    res1(i,j,k,l) = A(i,1,j,k,l);
    TEST(res1.compare_to_data({7,8,9,10,11,12,19,20,21,22,23,24}));
    
    res1(i,j,k,l) = A(i,j,k,l,0);
    TEST(res1.compare_to_data({1,3,5,7,9,11,13,15,17,19,21,23}));
    res1(i,j,k,l) = A(i,j,k,l,1);
    TEST(res1.compare_to_data({2,4,6,8,10,12,14,16,18,20,22,24}));
    
    res2(i,j,k) = A(i,j,1,k,1);
    TEST(res2.compare_to_data({4,10,16,22}));
    res2(i,j,k) = A(i,j,2,k,1);
    TEST(res2.compare_to_data({6,12,18,24}));
    res2(i,j,k) = A(j,i,2,k,1);
    TEST(res2.compare_to_data({6,18,12,24}));
    
    res2(i,k,j) = A(j,i,2,k,1);
    TEST(res2.compare_to_data({6,18,12,24}));
    
    res3(i,j) = A(j,i,2,0,1);
    TEST(res3.compare_to_data({6,18,12,24}));
    res3(i,j) = A(1,i,2,j,1);
    TEST(res3.compare_to_data({18,24}));
)

UNIT_TEST(SparseTensor, Assignment_Negatives,
    SparseTensor A({2,2,2,2});
    SparseTensor A2({2,2,2,2});
    SparseTensor B({2,2,2});
    SparseTensor C;
    SparseTensor D({2,2});
    SparseTensor E({});
    
    Index i,j,k,l;
        
    FAILTEST(A(i) = A(i));
    FAILTEST(A(i&1) = A(i&1));
    FAILTEST(A(i^3) = A(i^3));
    FAILTEST(A(i,j^3) = A(i^3,j));
    FAILTEST(A(i^5) = A(i&0));
    FAILTEST(A(i^4,j&5) = A(i^4,j&5));
    FAILTEST(A2(i) = A(i));
    FAILTEST(A2(i&1) = A(i&1));
    FAILTEST(A2(i^3) = A(i^3));
    FAILTEST(A2(i,j^3) = A(i^3,j));
    FAILTEST(A2(i^5) = A(i&0));
    FAILTEST(A2(i^4,j&5) = A(i^4,j&5));
    FAILTEST(A(i,j,k^2) = B(i,k^2));
    FAILTEST(A(i,j,k,l) = B(i,j,k));
    FAILTEST(A(i,j,k,l) = B(j,k,i));
    FAILTEST(C(i,j) = D(i^2));
    FAILTEST(C(i,j) = D(i^2));
	static_assert(!std::is_assignable<decltype(D(i,j) * D(j,k)), decltype(D(i,k))>::value,"");
)


    