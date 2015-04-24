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

UNIT_TEST(TT, TTTensor_Rounding,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);

    Index i,j,k;
    
    FullTensor A1 = FullTensor::construct_random({2}, rnd, dist);
    FullTensor B1(1);
    TTTensor TTA1(A1, 1e-14);
    B1(i) = TTA1(i);
    TEST(approx_equal(B1,A1, 1e-12));
    
    TTA1.round(1e-14);
    B1(i) = TTA1(i);
    TEST(approx_equal(B1,A1, 1e-12));
    
    TTA1.round(1);
    B1(i) = TTA1(i);
    TEST(approx_equal(B1,A1, 1e-12));
    
    
    FullTensor A2 = FullTensor::construct_random({2,2}, rnd, dist);
    FullTensor B2(2);
    TTTensor TTA2(A2, 1e-14);
    B2(j,i) = TTA2(j,i);
    TEST(approx_equal(B2,A2, 1e-12));
    
    TTA2.round(1e-14);
    B2(j,i) = TTA2(j,i);
    TEST(approx_equal(B2,A2, 1e-12));
    
    TTA2.round(2);
    B2(j,i) = TTA2(j,i);
    TEST(approx_equal(B2,A2, 1e-12));
    
    
    FullTensor A3 = FullTensor::construct_random({2,7}, rnd, dist);
    FullTensor B3(2);
    TTTensor TTA3(A3, 1e-14);
    B3(j,i) = TTA3(j,i);
    TEST(approx_equal(B3,A3, 1e-12));
    
    TTA3.round(1e-14);
    B3(j,i) = TTA3(j,i);
    TEST(approx_equal(B3,A3, 1e-12));
    
    TTA3.round(2);
    B3(j,i) = TTA3(j,i);
    TEST(approx_equal(B3,A3, 1e-12));
    
    
    FullTensor A4 = FullTensor::construct_random({2,2,2,2,2,2,2,2}, rnd, dist);
    FullTensor B4(8);
    TTTensor TTA4(A4, 1e-14);
    B4(j,i^7) = TTA4(j,i^7);
    TEST(approx_equal(B4,A4, 1e-12));
    
    TTA4.round(1e-14);
    B4(j,i^7) = TTA4(j,i^7);
    TEST(approx_equal(B4,A4, 1e-12));
    
    TTA4.round(512);
    B4(j,i^7) = TTA4(j,i^7);
    TEST(approx_equal(B4,A4, 1e-12));
    
    
    FullTensor A5 = FullTensor::construct_random({5,6,3,1,4,2,8,1}, rnd, dist);
    FullTensor B5(8);
    TTTensor TTA5(A5, 1e-14);
    B5(j,i^7) = TTA5(j,i^7);
    TEST(approx_equal(B5,A5, 1e-12));
    
    TTA5.round(1e-14);
    B5(j,i^7) = TTA5(j,i^7);
    TEST(approx_equal(B5,A5, 1e-12));
    
    TTA5.round(576);
    B5(j,i^7) = TTA5(j,i^7);
    TEST(approx_equal(B5,A5, 1e-12));
)

