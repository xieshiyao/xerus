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
	
    (res1(m,i,j), res2(n,m), res3(k,n,l)) = SVD(A(i,j,k,l));
    TEST(res2.compare_to_data(A.data.get()));
)

UNIT_TEST(FullTensor, SVD_Random_512x512,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::random({8,8,8,8,8,8}, rnd, dist);
    FullTensor res1(4);
    FullTensor res2(2);
    FullTensor res3(4);
    FullTensor res4(6);
     
    
    Index i, j, k, l, m, n, o, p, r, s;
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(i,j,k,l,m,n));
    res4(i,j,k,l,m,n) = res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-14));
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(l,k,m,i,j,n));
    res4(l,k,m,i,j,n) =  res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-14));
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(l,i,m,k,j,n));
    res4(k,i,m,l,j,n) =  res1(i,j,l,o)*res2(o,p)*res3(p,k,m,n);
    TEST(approx_equal(res4, A, 1e-14));
	
	(res1(i,o,k,j), res2(p,o), res3(l,n,m,p)) = SVD(A(l,i,m,k,j,n));
    res4(l,k,m,i,j,n) =  res1(k,o,i,j)*res2(p,o)*res3(l,n,m,p);
    TEST(approx_equal(res4, A, 1e-14));
)

UNIT_TEST(FullTensor, SVD_soft_thresholding,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::random({3,5,2,7,3,12}, rnd, dist);
    FullTensor Ax, U, V, Us, Vs;
    SparseTensor S, Ss;
    
    Index i, j, k, l, m, n, o, p, r, s;
    
    (U(i,j,k,o), S(o,p), V(p,l,m,n)) = SVD(A(i,j,k,l,m,n));
    (Us(i,j,k,o), Ss(o,p), Vs(p,l,m,n)) = SVD(A(i,j,k,l,m,n), 7.3);
	
	U.resize_dimension(U.degree()-1, Ss.dimensions[0]);
	V.resize_dimension(0, Ss.dimensions[0]);
	
    TEST(approx_equal(U, Us, 1e-12));
    TEST(approx_equal(V, Vs, 1e-12));
	
	for(size_t x = 0; x < S.dimensions[0]; ++x) {
		if(x < Ss.dimensions[0]) {
			TEST(misc::approx_equal(Ss[{x,x}], std::max(0.0, S[{x, x}]-7.3), 3e-13));
		} else {
			TEST(S[{x,x}] <= 7.3);
		}
	}
	
	Ax(i,j,k,l,m,n) = U(i,j,k,o)* S(o,p)* V(p,l,m,n);
	TEST(approx_equal(A, Ax, 1e-12));
)

UNIT_TEST(FullTensor, SVD_Random_Order_Six,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::random({9,7,5,5,9,7}, rnd, dist);
    FullTensor res1(4);
    FullTensor res2(2);
    FullTensor res3(4);
    FullTensor res4(6);
     
    
    Index i, j, k, l, m, n, o, p, r, s;
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(i,j,k,l,m,n));
    res4(i,j,k,l,m,n) = res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-14));
    
    (res1(i,j,k,o), res2(o,p), res3(p,l,m,n)) = SVD(A(m,j,k,l,i,n));
    res4(m,j,k,l,i,n) = res1(i,j,k,o)*res2(o,p)*res3(p,l,m,n);
    TEST(approx_equal(res4, A, 1e-14));
)

UNIT_TEST(FullTensor, QR_AND_RQ_Random_Order_Six,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    FullTensor A = FullTensor::random({7,5,9,7,5,9}, rnd, dist);
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
    MTEST(approx_equal(res4, A, 2e-15), "1 " << frob_norm(res4-A) << " / " << frob_norm(A));
	
	res4(l,m) = Q(i,j,k,l) * Q(i,j,k,m);
	res4.modify_diag_elements([](value_t &entry){entry -= 1;});
	TEST(frob_norm(res4) < 1e-12);
    
    (Q(i,j,k,l), R(l,m,n,r)) = QR(A(i,n,k,m,j,r));
    res4(i,n,k,m,j,r) = Q(i,j,k,o)*R(o,m,n,r);
    MTEST(approx_equal(res4, A, 2e-15), "2 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    (Q2(i,k,l), R2(l,m,j,n,r)) = QR(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q2(i,k,o)*R2(o,m,j,n,r);
    MTEST(approx_equal(res4, A, 2e-15), "3 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    (Q3(i,m,j,k,l), R3(l,n,r)) = QR(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q3(i,m,j,k,o)*R3(o,n,r);
    MTEST(approx_equal(res4, A, 2e-15), "4 " << frob_norm(res4-A) << " / " << frob_norm(A));
	
	(Q(i,l,j,k,m), R(l,n,r)) = QR(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q(i,o,j,k,m)*R(o,n,r);
    MTEST(approx_equal(res4, A, 2e-15), "5 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    
    (R(i,j,k,l), Q(l,m,n,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R(i,j,k,o)*Q(o,m,n,r);
    MTEST(approx_equal(res4, A, 5e-15), "6 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    (R(i,j,k,l), Q(l,m,n,r)) = RQ(A(i,n,k,m,j,r));
    res4(i,n,k,m,j,r) = R(i,j,k,o)*Q(o,m,n,r);
    MTEST(approx_equal(res4, A, 5e-15), "7 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    (R2(i,m,j,k,l), Q2(l,n,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R2(i,m,j,k,l)*Q2(l,n,r);
    MTEST(approx_equal(res4, A, 2e-15), "8 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    (R3(i,k,l), Q3(l,m,j,n,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R3(i,k,o)*Q3(o,m,j,n,r);
    MTEST(approx_equal(res4, A, 2e-15), "9 " << frob_norm(res4-A) << " / " << frob_norm(A));
	
	(R(l,i,k), Q(n,m,j,l,r)) = RQ(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = R(o,i,k)*Q(n,m,j,o,r);
    MTEST(approx_equal(res4, A, 2e-15), "10 " << frob_norm(res4-A) << " / " << frob_norm(A));
    
    
    (R3(i,k,l), Q4(l,j,n,r)) = RQ(A(i,j,k,3,n,r));
    res4(i,j,k,n,r) = R3(i,k,o)*Q4(o,j,n,r);
    TEST(frob_norm(A(i,j,k,3,n,r) - res4(i,j,k,n,r)) < 1e-12);
)


UNIT_TEST(FullTensor, QC,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 1.0);

    FullTensor A = FullTensor::random({2,2,2,2,2,2}, rnd, dist);
	FullTensor B({2,3}, [](size_t i){return double(i);});
    FullTensor Q(4);
    FullTensor R(4);
    FullTensor Q2(3);
    FullTensor R2(5);
    FullTensor Q3(5);
    FullTensor R3(3);
    FullTensor Q4(4);
    FullTensor res4(6);
    
    Index i, j, k, l, m, n, o, p, q, r;
	
	(Q(i,j), R(j,k)) = QC(B(i,k));
    
    (Q(i,j,k,l), R(l,m,n,r)) = QC(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q(i,j,k,o)*R(o,m,n,r);
    TEST(approx_equal(res4, A, 1e-15));
    
    (Q(i,j,k,l), R(l,m,n,r)) = QC(A(i,n,k,m,j,r));
    res4(i,n,k,m,j,r) = Q(i,j,k,o)*R(o,m,n,r);
    TEST(approx_equal(res4, A, 1e-15));
    
    (Q2(i,k,l), R2(l,m,j,n,r)) = QC(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q2(i,k,o)*R2(o,m,j,n,r);
    TEST(approx_equal(res4, A, 1e-12));
    
    (Q3(i,m,j,k,l), R3(l,n,r)) = QC(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q3(i,m,j,k,o)*R3(o,n,r);
    TEST(approx_equal(res4, A, 1e-15));
	
	(Q(i,l,j,k,m), R(l,n,r)) = QC(A(i,j,k,m,n,r));
    res4(i,j,k,m,n,r) = Q(i,o,j,k,m)*R(o,n,r);
    TEST(approx_equal(res4, A, 1e-15));
)
