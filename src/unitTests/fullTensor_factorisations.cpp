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

UNIT_TEST(FullTensor, SVD_Identity,
    FullTensor A({2,2,2,2});
    FullTensor res1({2,2,4});
    FullTensor res2({4,4});
    FullTensor res3({4,2,2});
    
    A[{0,0,0,0}] = 1;
    A[{0,1,0,1}] = 1;
    A[{1,0,1,0}] = 1;
    A[{1,1,1,1}] = 1;    
    
    Index i, j, k, l, m, n;
    
    (res1(i,j,m), res2(m,n), res3(n,k,l)) = SVD(A(i,j,k,l));
    TEST(res2.compare_to_data(A.data.get()));
)

UNIT_TEST(FullTensor, SVD_Random_512x512,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::construct_random({8,8,8,8,8,8}, rnd, dist);
    FullTensor res1(4);
    FullTensor res2(2);
    FullTensor res3(4);
    FullTensor res4(6);
     
    
    Index i, j, k, l, m, n, o, p, r, s;
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(i,j,k,l,m,n));
    res4(i,j,k,l,m,n) = res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-12));
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(l,k,m,i,j,n));
    res4(l,k,m,i,j,n) =  res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-12));
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(l,i,m,k,j,n));
    res4(k,i,m,l,j,n) =  res1(i,j,l,o)*res2(o,p)*res3(p,k,m,n);
    TEST(approx_equal(res4, A, 1e-12));
)

UNIT_TEST(FullTensor, SVD_Random_Order_Six,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::construct_random({9,7,5,5,9,7}, rnd, dist);
    FullTensor res1(4);
    FullTensor res2(2);
    FullTensor res3(4);
    FullTensor res4(6);
     
    
    Index i, j, k, l, m, n, o, p, r, s;
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(i,j,k,l,m,n));
    res4(i,j,k,l,m,n) = res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-12));
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(m,j,k,l,i,n));
    res4(m,j,k,l,i,n) = res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-12));
)

UNIT_TEST(FullTensor, QR_AND_RQ_Random_Order_Six,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::construct_random({7,5,9,7,5,9}, rnd, dist);
    FullTensor Q(4);
    FullTensor R(4);
    FullTensor Q2(3);
    FullTensor R2(5);
    FullTensor Q3(5);
    FullTensor R3(3);
    FullTensor Q4(4);
    FullTensor res4(6);
    
    Index i, j, k, l, m, n, o, p, q, r;

    
    (Q(i,j,k,l), R(l,m,n,r)) = QR(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q(i,j,k,o)*R(o,m,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (Q(i,j,k,l), R(l,m,n,r)) = QR(A(i,n,k,m,j,r));
    res4(i,n,k,m,j,r) = Q(i,j,k,o)*R(o,m,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (Q2(i,k,l), R2(l,m,j,n,r)) = QR(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q2(i,k,o)*R2(o,m,j,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (Q3(i,m,j,k,l), R3(l,n,r)) = QR(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q3(i,m,j,k,o)*R3(o,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    
    (R(i,j,k,l), Q(l,m,n,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R(i,j,k,o)*Q(o,m,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (R(i,j,k,l), Q(l,m,n,r)) = RQ(A(i,n,k,m,j,r));
    res4(i,n,k,m,j,r) = R(i,j,k,o)*Q(o,m,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (R2(i,m,j,k,l), Q2(l,n,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R2(i,m,j,k,l)*Q2(l,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (R3(i,k,l), Q3(l,m,j,n,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R3(i,k,o)*Q3(o,m,j,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    
    (R3(i,k,l), Q4(l,j,n,r)) = RQ(A(i,j,k,3,n,r));
    res4(i,j,k,n,r) = R3(i,k,o)*Q4(o,j,n,r);
    TEST(frob_norm(A(i,j,k,3,n,r) - res4(i,j,k,n,r)) < 1e-12);
)



// UNIT_TEST(FullTensor, Factorisations_Plus_Arithmetic,
//     std::mt19937_64 rnd;
//     std::normal_distribution<value_t> dist (0.0, 10.0);
// 
//     FullTensor A = FullTensor::construct_random({7,5,9,7,5,9}, rnd, dist);
//     FullTensor B = FullTensor::construct_random({7,5,9,7,5,9}, rnd, dist);
//     FullTensor C = FullTensor::construct_random({7,5,9,7,5,9}, rnd, dist);
//     FullTensor Q1(4);
//     FullTensor Q2(4);
//     FullTensor Q3(4);
//     FullTensor R1(4);
//     FullTensor R2(4);
//     FullTensor R3(4);
//     FullTensor res1(6);
//     FullTensor res2(6);
//      
// )
