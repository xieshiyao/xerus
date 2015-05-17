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

UNIT_TEST(TensorNetwork, element_access,
    FullTensor A({1,2});
    FullTensor B({2,3});
    TensorNetwork res;

    Index i,j,k,l;
    
    A[{0,0}] = 1;
    A[{0,1}] = 2;
    
    B[{0,0}] = 3;
    B[{0,1}] = 4;
    B[{0,2}] = 5;
    B[{1,0}] = 6;
    B[{1,1}] = 7;
    B[{1,2}] = 8;
    
    //No Index contracted
    res(i,j,k,l) = A(i,j) * B(k,l);
    TEST(approx_equal(res[0], 3.0));
    TEST(approx_equal(res[1], 4.0));
    
    /*
    TEST(res1.compare_to_data({3,4,5,6,7,8,6,8,10,12,14,16}));
    res1(i,j,k,l) = A(i,k) * B(j,l);
    TEST(res1.compare_to_data({3,4,5,6,8,10,6,7,8,12,14,16}));
    res1(i,k,j,l) = A(i,j) * B(k,l);
    TEST(res1.compare_to_data({3,4,5,6,8,10,6,7,8,12,14,16}));
    
    //One Index contracted
    res2(i,k) = A(i,j) * B(j,k);
    TEST(res2.compare_to_data({15,18,21}));
    res3(k,i) = A(i,j) * B(j,k);
    TEST(res2.compare_to_data({15,18,21}));
    res2(i,k) = B(j,k) * A(i,j);
    TEST(res2.compare_to_data({15,18,21}));
    res3(k,i) = B(j,k) * A(i,j);
    TEST(res2.compare_to_data({15,18,21}));*/
)