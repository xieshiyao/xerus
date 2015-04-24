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

UNIT_TEST(SparseTensor, Creation,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0); 
    
    FullTensor fullA = FullTensor::construct_random({7,13,2,9,3}, rnd, dist);
    FullTensor fullB = FullTensor::construct_random({7,13,2,9,3}, rnd, dist);
    FullTensor fullX({7,13,2,9,3});
    
    SparseTensor sparseA = SparseTensor(fullA);
    SparseTensor sparseB = SparseTensor(fullB);
    SparseTensor sparseX({7,13,2,9,3});
    
    
    TEST(approx_equal(fullA, sparseA, 6e-14));
    TEST(approx_equal(fullB, sparseB, 6e-14));
    TEST(approx_equal(fullA, FullTensor(sparseA), 6e-14));
    TEST(approx_equal(fullB, FullTensor(sparseB), 6e-14));
    
    fullX += fullA;
    sparseX += sparseA;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX -= fullB;
    sparseX -= sparseB;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX = fullA + fullB;
    sparseX = sparseA + sparseB;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX = fullA - fullB;
    sparseX = sparseA - sparseB;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX = 2.0*fullA;
    sparseX = 2.0*sparseA;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX = 10.0*fullA*2;
    sparseX = 10.0*sparseA*2;
    TEST(approx_equal(fullX, sparseX, 6e-14));
     
    fullX = fullA/10.0;
    sparseX = sparseA/10.0;
    TEST(approx_equal(fullX, sparseX, 6e-14));

    fullX = 0*fullA + fullB;
    sparseX = 0*sparseA + sparseB;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX = 7.3*fullA + fullB*5;
    sparseX = 7.3*sparseA + sparseB*5;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    fullX = 7.9*fullA/13.7 + fullB*5;
    sparseX = 7.9*sparseA/13.7 + sparseB*5;
    TEST(approx_equal(fullX, sparseX, 6e-14));
    
    TEST(approx_equal(fullA, sparseA, 6e-14));
    TEST(approx_equal(fullB, sparseB, 6e-14));
    
    // Two times to be sure
    
    fullX += fullA;
    sparseX += sparseA;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX -= fullB;
    sparseX -= sparseB;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX = fullA + fullB;
    sparseX = sparseA + sparseB;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX = fullA - fullB;
    sparseX = sparseA - sparseB;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX = 2.0*fullA;
    sparseX = 2.0*sparseA;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX = 10.0*fullA*2;
    sparseX = 10.0*sparseA*2;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
     
    fullX = fullA/10.0;
    sparseX = sparseA/10.0;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));

    fullX = 0*fullA + fullB;
    sparseX = 0*sparseA + sparseB;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX = 7.3*fullA + fullB*5;
    sparseX = 7.3*sparseA + sparseB*5;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    fullX = 7.9*fullA/13.7 + fullB*5;
    sparseX = 7.9*sparseA/13.7 + sparseB*5;
    TEST(approx_equal(fullX, FullTensor(sparseX), 6e-14));
    
    TEST(approx_equal(fullA, sparseA, 6e-14));
    TEST(approx_equal(fullB, sparseB, 6e-14));
    
)
