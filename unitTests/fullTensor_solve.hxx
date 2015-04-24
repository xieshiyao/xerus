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

UNIT_TEST(FullTensor, solve_Ax_equals_b,
    Index i,j,k;
          
    FullTensor A1({4,2,2});
    A1[{0,0,0}] = 1;
    A1[{1,0,1}] = 1;
    A1[{2,1,0}] = 1;
    A1[{3,1,1}] = 1;
    
    FullTensor b1({4});
    b1[{0}] = 73;
    b1[{1}] = -73;
    b1[{2}] = 128;
    b1[{3}] = 93;
    
    FullTensor x1({2,2});
    
    
    x1(k,i) = (b1(j)/A1(j,k,i));
    
    TEST((b1[{0}] - x1[{0,0}]) < 1e-14);
    TEST((b1[{1}] - x1[{0,1}]) < 1e-14);
    TEST((b1[{2}] - x1[{1,0}]) < 1e-14);
    TEST((b1[{3}] - x1[{1,1}]) < 1e-14);
    
    FullTensor A2({4,2,2});
    A2[{0,0,0}] = 1;
    A2[{1,0,1}] = 1;
    A2[{2,1,0}] = 0;
    A2[{3,1,1}] = 0;
    
    FullTensor b2({4});
    b2[{0}] = 73;
    b2[{1}] = -73;
    b2[{2}] = 0;
    b2[{3}] = 0;
    
    FullTensor x2({2,2});
    
    x2(k,i) = (b2(j)/A2(j,k,i));
    
    TEST((b2[{0}] - x2[{0,0}]) < 1e-14);
    TEST((b2[{1}] - x2[{0,1}]) < 1e-14);
    TEST((x2[{1,0}]) < 1e-14);
    TEST((x2[{1,1}]) < 1e-14);
)
