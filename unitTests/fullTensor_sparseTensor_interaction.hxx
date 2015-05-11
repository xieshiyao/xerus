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

#include <type_traits>

UNIT_TEST(FullTensor_SparseTensor_Interaction, Basics, 
    SparseTensor A({2,2,3,1,2});
    FullTensor B;
    FullTensor resF;
    SparseTensor resS;
    
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
    
    
    B(j,i,k,l,m) = A(i,j,k,l,m);
    TEST(B.compare_to_data({1,2,3,4,5,6,13,14,15,16,17,18,7,8,9,10,11,12,19,20,21,22,23,24}));
    
    resF = B + A;
    TEST(resF.compare_to_data({1+1,2+2,3+3,4+4,5+5,6+6,13+7,14+8,15+9,16+10,17+11,18+12,7+13,8+14,9+15,10+16,11+17,12+18,19+19,20+20,21+21,22+22,23+23,24+24}));
    
    resF = B - A;
    TEST(resF.compare_to_data({1-1,2-2,3-3,4-4,5-5,6-6,13-7,14-8,15-9,16-10,17-11,18-12,7-13,8-14,9-15,10-16,11-17,12-18,19-19,20-20,21-21,22-22,23-23,24-24}));
    
)