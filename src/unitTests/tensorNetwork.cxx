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

UNIT_TEST(TensorNetwork, contractions_of_4_to_degree_0,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Tensor A = Tensor::random({100,1}, rnd, dist);
	Tensor B = Tensor::random({100,1}, rnd, dist);
	Tensor C = Tensor::random({100,1}, rnd, dist);
	Tensor D = Tensor::random({100,1}, rnd, dist);
	Tensor E;
	Index i1,i2,i3,i4;
	
	
	E() = A(i1,i2) * D(i1,i2);
	double a1 = E[{}];
	E() = B(i3,i4) * C(i3,i4);
	double a2 = E[{}];
	E() = A(i1,i2) * B(i3,i2) * C(i3,i4) * D(i1,i4);
	TEST(misc::approx_equal(E[{}], a1 * a2, 1e-20));
	E() = B(i3,i2) * C(i3,i4) * D(i1,i4) * A(i1,i2);
	TEST(misc::approx_equal(E[{}], a1 * a2, 1e-20));
	E() = B(i3,i2) * D(i1,i4) * C(i3,i4) * A(i1,i2);
	TEST(misc::approx_equal(E[{}], a1 * a2, 1e-20));
)

UNIT_TEST(TensorNetwork, contractions_of_3_to_degree_0,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Tensor A = Tensor::random({1,10}, rnd, dist);
	Tensor B = Tensor::random({10,100}, rnd, dist);
	Tensor C = Tensor::random({100,1}, rnd, dist);
	Tensor E;
	Index i1,i2,i3,i4;
	
	
	E() = A(i1,i2) * B(i2,i3) * C(i3,i1);
	double a1 = E[{}];
	E() = B(i2,i3) * C(i3,i1) * A(i1,i2);
	double a2 = E[{}];
	E() = C(i3,i1) * B(i2,i3) * A(i1,i2);
	double a3 = E[{}];
	LOG(unit_test, a1 << " " << a2 << " " << a3 << " " << a1-a2 << " " << a2-a3);
	TEST(misc::approx_equal(a1, a2, 1e-20*a1)); // effectively ==
	TEST(misc::approx_equal(a2, a3, 1e-20*a2));
)

UNIT_TEST(TensorNetwork, contraction_single_node_trace,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Tensor A = Tensor::random({1,10,10}, rnd, dist);
	Tensor B = Tensor::random({1}, rnd, dist);
	Tensor E;
	Index i1,i2,i3;
	
 	E() = A(i1,i2,i2) * B(i1);
	TEST(std::isnormal(E[{}]));
)

UNIT_TEST(TensorNetwork, contraction_single_network_trace,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Tensor A = Tensor::random({2,2,2,2}, rnd, dist);
	TensorNetwork ATN(A);
	Tensor E;
	TensorNetwork ETN(0);
	Index i1,i2,i3,i4;
	
	ETN() = ATN(i1,i1,i2,i2);
	E = Tensor(ETN);
	TEST(std::isnormal(E[{}]));
 	E() = ATN(i1,i1,i2,i2);
	TEST(std::isnormal(E[{}]));
	
)


UNIT_TEST(TensorNetwork, index_reshuffle2,
	TensorNetwork A, B;
	Index n1,n2,n3,n4,r1,r2;
	Tensor At({2,1,3});
	Tensor Bt({1,4,5,1});
	Tensor Bt2({1,6,7,1});
	A = At;
	B = Bt;
	A.draw(std::string("out")+misc::to_string(0));
	A(n1^(1), n2, r2, n3^(1), n4) = A(n1^(1), r1, n3^(1)) * B(r1, n2, n4, r2);
	A.draw(std::string("out")+misc::to_string(1));
	TEST(A.dimensions == std::vector<size_t>({2,4,1,3,5}));
	B = Bt2;
	A(n1^(2), n2, r2, n3^(2), n4) = A(n1^(2), r1, n3^(2)) * B(r1, n2, n4, r2);
	A.draw(std::string("out")+misc::to_string(2));
	MTEST(A.dimensions == std::vector<size_t>({2,4,6,1,3,5,7}), "found " << A.dimensions);
)


UNIT_TEST(TensorNetwork, triple_indices,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

	TensorNetwork A(3);
	TensorNetwork B(2);
	TensorNetwork C(2);
	TensorNetwork D(2);
	TensorNetwork F(2);
	Tensor E0;
	Tensor E1;
	Tensor E2;
	Index i1,i2,i3,i4;
	
	FAILTEST(E0()   = A(i1,i1,i2)*B(i2,i2));
	FAILTEST(E1(i2) = A(i1,i1,i2)*B(i2,i2));
	FAILTEST(E0()   = A(i1,i2,i2)*B(i2,i1));
	FAILTEST(E1(i2) = A(i1,i2,i2)*B(i2,i1));
	FAILTEST(E0()   = A(i2,i2,i2)*B(i1,i1));
	FAILTEST(E1(i2) = A(i2,i2,i2)*B(i1,i1));
	FAILTEST(E0()   = A(i1,i2,i2)*B(i1,i3)*C(i3,i2));
// 	FAILTEST(E1(i2) = A(i1,i2,i2)*B(i1,i3)*C(i3,i2)); //FEATURE
	FAILTEST(E0()      = B(i1,i2)*C(i2,i3)*D(i3,i2));
// 	FAILTEST(E2(i1,i2) = B(i1,i2)*C(i2,i3)*D(i3,i2)); //FEATURE
	FAILTEST(E0()      = B(i1,i2)*C(i2,i3)*D(i1,i2));
// 	FAILTEST(E2(i2,i3) = B(i1,i2)*C(i2,i3)*D(i1,i2)); //FEATURE
	FAILTEST(E0()      = B(i1,i2)*C(i2,i3)*D(i3,i4)*F(i4,i2));
// 	FAILTEST(E2(i1,i2) = B(i1,i2)*C(i2,i3)*D(i3,i4)*F(i4,i2)); //FEATURE
)

UNIT_TEST(TensorNetwork, contraction_multi_node_trace,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Tensor A = Tensor::random({1,10}, rnd, dist);
	Tensor B = Tensor::random({1,10}, rnd, dist);
	Tensor E;
	Index i1,i2,i3,i4;
	
	TensorNetwork tmp(4);
	tmp(i1,i2,i3,i4) = A(i1,i3) * B(i2,i4);
	E() = tmp(i1,i1,i2,i2);
	TEST(std::isnormal(E[{}]));
)

UNIT_TEST(TensorNetwork, index_reshuffle,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Tensor A = Tensor::random({1,10}, rnd, dist);
	Tensor B = Tensor::random({1,10}, rnd, dist);
	Tensor E;
	Index i1,i2,i3,i4;
	
	TensorNetwork tmp(4);
	tmp(i1,i2,i3,i4) = A(i1,i3) * B(i2,i4);
	tmp(i1,i2,i3,i4) = tmp(i3,i4,i1,i2);
	E() = tmp(i1,i1,i2,i2);
	TEST(std::isnormal(E[{}]));
)

UNIT_TEST(TensorNetwork, Save_Network, 
    Tensor A({2,2});
    Tensor B({2,2});
    Tensor C({2,2});
    Tensor D({2,2});
    Tensor E({2,2});
    Tensor F({2,2});
    Tensor G({2,2});
    TensorNetwork res1(2);
    TensorNetwork res1A(6);
    Tensor res1AF;
    TensorNetwork res2(2);
    TensorNetwork res2A(6);
    Tensor res2AF;
    TensorNetwork res2B(4);
    Tensor res3({2,2});

    Index i,j,k,l,m,n,o;
    
    A[{0,0}]=1;
    A[{0,1}]=2;
    A[{1,0}]=3;
    A[{1,1}]=4;
    
    B[{0,0}]=5;
    B[{0,1}]=6;
    B[{1,0}]=7;
    B[{1,1}]=8;
    
    C[{0,0}]=9;
    C[{0,1}]=10;
    C[{1,0}]=11;
    C[{1,1}]=12;
    
    D[{0,0}]=13;
    D[{0,1}]=14;
    D[{1,0}]=15;
    D[{1,1}]=16;
    
    E[{0,0}]=17;
    E[{0,1}]=18;
    E[{1,0}]=19;
    E[{1,1}]=20;
    
    F[{0,0}]=21;
    F[{0,1}]=22;
    F[{1,0}]=23;
    F[{1,1}]=24;
     
    res2(i,l) = A(i,j) * B(j,k) * C(k,l); 
	res1 = std::move(res2); // res2 should still be valid but not change res1 in the following
    res2(l,o) = D(l,m) * E(m,n) * F(n,o);
    res3(i,o) = res1(i,l) * res2(l,o);
    TEST(approx_entrywise_equal(res3, {20596523, 21531582, 46728183, 48849590}));
    
    res1A(i,j,m,n,k,l) = A(i,j) * E(m,n) * C(k,l); 
    res2A(l,m,j,k,n,o) = D(l,m) * B(j,k) * F(n,o);
    res3(i,o) = res1A(i,j,m,n,k,l) * res2A(l,m,j,k,n,o);
    TEST(approx_entrywise_equal(res3, {20596523, 21531582, 46728183, 48849590}));
	
	res1A(i,j,m,n,k,l) = A(i,j) * E(m,n) * C(k,l); 
    res2A(l,m,j,k,n,o) = D(l,m) * B(j,k) * F(n,o);
    res3(i,o) = res1A(i,j,m,n,k,l) * res2A(l,m,j,k,n,o);
    TEST(approx_entrywise_equal(res3, {20596523, 21531582, 46728183, 48849590}));
    
    res1A(i,l,m,n,j,k) = A(i,j) * E(m,n) * C(k,l); 
    res2A(l,o,m,n,j,k) = D(l,m) * B(j,k) * F(n,o);
    res3(i,o) = res1A(i,l,m,n,j,k) * res2A(l,o,m,n,j,k);
    TEST(approx_entrywise_equal(res3, {20596523, 21531582, 46728183, 48849590}));

    res1AF(i,l,m,n,j,k) = A(i,j) * E(m,n) * C(k,l); 
    res2AF(l,o,m,n,j,k) = D(l,m) * B(j,k) * F(n,o);
    res3(i,o) = res1AF(i,l,m,n,j,k) * res2AF(l,o,m,n,j,k);
    TEST(approx_entrywise_equal(res3, {20596523, 21531582, 46728183, 48849590}));
    
    res1A(i,l,m,n,j,k) = A(i,j) * E(m,n) * C(k,l); 
    res2A(l,o,m,n,j,k) = D(l,m) * B(j,k) * F(n,o);
    res3(i,o) = res1A(i,l,m,n,j,k) * res2A(l,o,m,n,j,k);
    TEST(approx_entrywise_equal(res3, {20596523, 21531582, 46728183, 48849590}));
)
