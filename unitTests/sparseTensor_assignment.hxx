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

UNIT_TEST(SparseTensor, Assignment_Trivia2, 
    SparseTensor A({2,2,3,1,2});
    SparseTensor res({2,2,3,1,2});

    Index i,j,k,l,m;
    
    A[{0,1,0,0,0}] = 73;
    
    res(i,j,k,l,m) = A(i,j,k,l,m);
    TEST(approx_equal(res[{0,1,0,0,0}], 73.0, 1e-14));
    
    res(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(approx_equal(res[{1,0,0,0,0}], 73.0, 1e-14));
)

UNIT_TEST(SparseTensor, Assignment_Sparse_To_Sparse, 
    SparseTensor A({2,2,3,1,2});
    SparseTensor res({2,2,3,1,2});
    SparseTensor res2({2,3,2,1,2});
    SparseTensor res3({2,3,1,2,2});

    Index i,j,k,l,m;
    
    A[{0,0,0,0,0}]=1;
    A[{0,0,1,0,1}]=4;
    A[{0,1,0,0,1}]=8;
    A[{0,1,2,0,0}]=11;
    A[{1,0,0,0,1}]=14;
    A[{1,1,0,0,0}]=19;
    A[{1,1,2,0,0}]=23;
    
    res(i,j,k,l,m) = A(i,j,k,l,m);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    
    res(i&0) = A(i&0);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    
    res(i&0) = A(i^5);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    
    res(i^5) = A(i&0);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    
    res(i^3, j&3) = A(i&2,j^2);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    
    res(j,i,k,l,m) = A(j,i,k,l,m);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    res(i,j,k,l,m) = A(j,i,k,l,m);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,14,0,0,0,0,0,8,0,0,11,0,19,0,0,0,23,0}));
    res(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(FullTensor(res).compare_to_data({1,0,0,4,0,0,0,14,0,0,0,0,0,8,0,0,11,0,19,0,0,0,23,0}));
    
    res2(i,k,j,l,m) = A(i,j,k,l,m);
    TEST(FullTensor(res2).compare_to_data({1,0,0,8,0,4,0,0,0,0,11,0,0,14,19,0,0,0,0,0,0,0,23,0}));
    res2(i,j,k,l,m) = A(i,k,j,l,m);
    TEST(FullTensor(res2).compare_to_data({1,0,0,8,0,4,0,0,0,0,11,0,0,14,19,0,0,0,0,0,0,0,23,0}));
)

UNIT_TEST(SparseTensor, Assignment_Sparse_To_Full, 
    SparseTensor A({2,2,3,1,2});
    FullTensor res(5);
    FullTensor res2(4);
    FullTensor res3(3);

    Index i,j,k,l,m;
    
    A[{0,0,0,0,0}]=1;
    A[{0,0,1,0,1}]=4;
    A[{0,1,0,0,1}]=8;
    A[{0,1,2,0,0}]=11;
    A[{1,0,0,0,1}]=14;
    A[{1,1,0,0,0}]=19;
    A[{1,1,2,0,0}]=23;
    
    res(i,j,k,l,m) = A(i,j,k,l,m);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    res(i&0) = A(i&0);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    res(i&0) = A(i^5);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    res(i^5) = A(i&0);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    res(i^3, j&3) = A(i&2,j^2);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    
    res(j,i,k,l,m) = A(j,i,k,l,m);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,8,0,0,11,0,0,14,0,0,0,0,19,0,0,0,23,0}));
    res(i,j,k,l,m) = A(j,i,k,l,m);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,14,0,0,0,0,0,8,0,0,11,0,19,0,0,0,23,0}));
    res(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(res.compare_to_data({1,0,0,4,0,0,0,14,0,0,0,0,0,8,0,0,11,0,19,0,0,0,23,0}));
    
    res(i,k,j,l,m) = A(i,j,k,l,m);
    TEST(res.compare_to_data({1,0,0,8,0,4,0,0,0,0,11,0,0,14,19,0,0,0,0,0,0,0,23,0}));
    res(i,j,k,l,m) = A(i,k,j,l,m);
    TEST(res.compare_to_data({1,0,0,8,0,4,0,0,0,0,11,0,0,14,19,0,0,0,0,0,0,0,23,0}));
    
    res2(i,k,j,l) = A(i,j,k,0,l);
    TEST(res2.compare_to_data({1,0,0,8,0,4,0,0,0,0,11,0,0,14,19,0,0,0,0,0,0,0,23,0}));
    res2(i,j,k,l) = A(i,k,j,0,l);
    TEST(res2.compare_to_data({1,0,0,8,0,4,0,0,0,0,11,0,0,14,19,0,0,0,0,0,0,0,23,0}));
    
    res3(i,k,j) = A(l,l,i,j,k);
    TEST(res3.compare_to_data({20,0,0,4,23,0}));
)
