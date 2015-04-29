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

UNIT_TEST(SparseTensor, Arithmetic_Negatives,
    SparseTensor A({2,2,2});
    SparseTensor B({2,2});
    SparseTensor B2({3,3});
    SparseTensor B3({3,2});
    SparseTensor B4({2,3});
    SparseTensor B5({2,2,2});
    SparseTensor C({2});
    
    Index i,j,k;

    FAILTEST(C(i) = B(i,j) * B2(j,j));
    FAILTEST(C(i) = B(i,j) * B3(j,j));
    FAILTEST(C(i) = B(i,j) * B4(j,j));
    FAILTEST(C(i) = B(i,j) * B5(j,j,j));
    FAILTEST(B(i,j) = B(i,j) + B2(j,j));
    FAILTEST(B(i,j) = B(i,j) + B3(j,j));
    FAILTEST(B(i,j) = B(i,j) + B4(j,j));
    FAILTEST(B(i,j) = B(i,j) + B5(j,j,j));
)

UNIT_TEST(SparseTensor, triple_indices,
	std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

	SparseTensor A({1,1,1});
	SparseTensor B({1,1});
	SparseTensor C({1,1});
	SparseTensor D({1,1});
	SparseTensor F({1,1});
	SparseTensor E0;
	SparseTensor E1({1});
	SparseTensor E2({1});
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
    