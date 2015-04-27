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

UNIT_TEST(FullTensor, Traces,
    FullTensor A({2,2});
    FullTensor B({2,2,2});
    FullTensor C({2,2,2,2});
    FullTensor res1({});
    FullTensor res2({2});
    FullTensor res3({2,2});
    
    Index i, j, k;
    
    A[{0,0}] = 1;
    A[{0,1}] = 2;
    A[{1,0}] = 4;
    A[{1,1}] = 8;
    
    
    B[{0,0,0}] = 1;
    B[{0,0,1}] = 2;
    B[{0,1,0}] = 4;
    B[{0,1,1}] = 8;
    B[{1,0,0}] = 16;
    B[{1,0,1}] = 32;
    B[{1,1,0}] = 64;
    B[{1,1,1}] = 128;
    
    C[{0,0,0,0}] = 1;
    C[{0,0,0,1}] = 2;
    C[{0,0,1,0}] = 4;
    C[{0,0,1,1}] = 8;
    C[{0,1,0,0}] = 16;
    C[{0,1,0,1}] = 32;
    C[{0,1,1,0}] = 64;
    C[{0,1,1,1}] = 128;
    C[{1,0,0,0}] = 256;
    C[{1,0,0,1}] = 512;
    C[{1,0,1,0}] = 1024;
    C[{1,0,1,1}] = 2048;
    C[{1,1,0,0}] = 4096;
    C[{1,1,0,1}] = 8192;
    C[{1,1,1,0}] = 16384;
    C[{1,1,1,1}] = 32768;
    
    res1() = A(i,i);
    TEST(res1.compare_to_data({9}));
    
    res2(j) = B(i,i,j);
    TEST(res2.compare_to_data({65, 130}));
    res2(j) = B(i,j,i);
    TEST(res2.compare_to_data({33, 132}));
    res2(j) = B(j,i,i);
    TEST(res2.compare_to_data({9, 144}));
    
    res3(j,k) = C(i,i,j,k);
    TEST(res3.compare_to_data({4097, 8194, 16388, 32776}));
    res3(j,k) = C(i,j,i,k);
    TEST(res3.compare_to_data({1025, 2050, 16400, 32800}));
    res3(j,k) = C(i,j,k,i);
    TEST(res3.compare_to_data({513, 2052, 8208, 32832}));
    res3(j,k) = C(j,i,i,k);
    TEST(res3.compare_to_data({65, 130, 16640, 33280}));
    res3(j,k) = C(j,i,k,i);
    TEST(res3.compare_to_data({33, 132, 8448, 33792}));
    res3(j,k) = C(j,k,i,i);
    TEST(res3.compare_to_data({9, 144, 2304, 36864}));
    
    res3(k,j) = C(i,i,j,k);
    TEST(res3.compare_to_data({4097, 16388, 8194, 32776}));
    res3(k,j) = C(i,j,i,k);
    TEST(res3.compare_to_data({1025, 16400, 2050, 32800}));
    res3(k,j) = C(i,j,k,i);
    TEST(res3.compare_to_data({513, 8208, 2052, 32832}));
    res3(k,j) = C(j,i,i,k);
    TEST(res3.compare_to_data({65, 16640, 130, 33280}));
    res3(k,j) = C(j,i,k,i);
    TEST(res3.compare_to_data({33, 8448, 132, 33792}));
    res3(k,j) = C(j,k,i,i);
    TEST(res3.compare_to_data({9, 2304, 144, 36864}));
    
    res2(k) = C(0,k,i,i);
    TEST(res2.compare_to_data({9, 144}));
    res2(k) = C(1,k,i,i);
    TEST(res2.compare_to_data({2304,36864}));
    res2(j) = C(j,0,i,i);
    TEST(res2.compare_to_data({9, 2304}));
    res2(j) = C(j,1,i,i);
    TEST(res2.compare_to_data({144, 36864}));
    
    res1() = C(i,i,j,j);
    TEST(res1.compare_to_data({1+8+4096+32768}));
    res1() = C(i,j,i,j);
    TEST(res1.compare_to_data({1+32+1024+32768}));
    res1() = C(i,j,j,i);
    TEST(res1.compare_to_data({1+64+512+32768}));
)
