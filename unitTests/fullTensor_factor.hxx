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

#pragma once
#include "../xerus.h"

UNIT_TEST(FullTensor, Factors,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::construct_random({2,7,5,5,2,7}, rnd, dist);
    FullTensor B = FullTensor::construct_random({2,7,5,5,2,7}, rnd, dist);
    FullTensor A3 = 3*A;
    FullTensor B7 = 7*B;
	TEST(!A.has_factor());
	TEST(!B.has_factor());
	TEST(A3.has_factor());
	TEST(B7.has_factor());
	
    FullTensor Q(4);
    FullTensor R(4);
    FullTensor res1(4);
    FullTensor res2(2);
    FullTensor res3(4);
    FullTensor res4(6);
    
	
    
    Index i, j, k, l, m, n, o, p, r, s;
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(3*A(i,j,k,l,m,n));
    res4(i,j,k,l,m,n) = 3.7*res1(i,j,k,o)*(res2(o,p)/3.7)*res3(p,l,m,n);
    TEST(approx_equal(res4, A3, 1e-12));
    
    (Q(i,j,k,l), R(l,m,n,r)) = QR(B7(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = (Q(i,j,k,o)/12.5)*(12.5*R(o,m,n,r)/7);
    TEST(approx_equal(res4, B, 1e-12));
)

UNIT_TEST(FullTensor, value_t_Product,
    FullTensor A({4,2,2,7});
    FullTensor B(4);
    FullTensor C(4);
    FullTensor D(4);
    for(size_t i = 0; i < A.size; ++i) {
        A[i] = 73;
    }
    
    Index j;
    
    B(j&0) = A(j&0)*2.0;
    C(j&0) = 3*A(j&0);
    D(j&0) = A(j&0)/73.0;
    A(j&0) = A(j&0)/2;
    
    TEST(B.compare_data(std::vector<value_t>(A.size, 146)));
    TEST(C.compare_data(std::vector<value_t>(A.size, 219)));
    TEST(D.compare_data(std::vector<value_t>(A.size, 1)));
    TEST(A.compare_data(std::vector<value_t>(A.size, 36.5)));
)