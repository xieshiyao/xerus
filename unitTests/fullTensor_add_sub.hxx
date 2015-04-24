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

UNIT_TEST(TensorSum, matrix_2x2,
    FullTensor res({2,2});
    FullTensor B({2,2});
    FullTensor C({2,2});

    Index i, J;
    
    B[{0,0}]=1;
    B[{0,1}]=2;
    B[{1,0}]=3;
    B[{1,1}]=4;
    
    C[{0,0}]=5;
    C[{0,1}]=6;
    C[{1,0}]=7;
    C[{1,1}]=8;
    
    res(i,J) = B(i,J) + C(i,J);
    TEST(compare_memory_to_vector(res.data.get(), {6,8,10,12}));
    res(i,J) = B(i,J) + C(J,i);
    TEST(compare_memory_to_vector(res.data.get(), {6,9,9,12}));
)
 
UNIT_TEST(TensorSum, lhs_equals_rhs,
    FullTensor B({2,2});
    FullTensor C({2,2});

    Index i, J;
    
    B[{0,0}]=1;
    B[{0,1}]=2;
    B[{1,0}]=3;
    B[{1,1}]=4;
    
    C[{0,0}]=5;
    C[{0,1}]=6;
    C[{1,0}]=7;
    C[{1,1}]=8;
    
    B(i,J) = B(i,J) + C(i,J);
    ASSERT(compare_memory_to_vector(B.data.get(), {6,8,10,12}));
    B(i,J) = B(i,J) + B(J,i);
    ASSERT(compare_memory_to_vector(B.data.get(), {12,18,18,24}));
)

namespace ___I_AM_A_NEW_NAMESPACE____________ {
_pure_  value_t filler1(const std::vector<size_t> &_idx) {
    return double(_idx[0] + _idx[1]);
}
_pure_  value_t filler2(const std::vector<size_t> &_idx) {
    return double(_idx[0] * _idx[1]);
}
_pure_ value_t filler12(const std::vector<size_t> &_idx) {
    return double(_idx[0] + _idx[1] + _idx[0] * _idx[1]);
}

UNIT_TEST(TensorSum, matrix_1000x1000,
    FullTensor res({1024,1024});
    FullTensor A({1024,1024}, filler1);
    FullTensor B({1024,1024}, filler2);
    FullTensor C({1024,1024}, filler12);

    Index i, J;
    
    res(i,J) = A(i,J) + B(i,J);
    TEST(memcmp(res.data.get(), C.data.get(), sizeof(value_t)*1024*1024)==0);
    res(J,i) = A(J,i) + B(i,J);
    TEST(memcmp(res.data.get(), C.data.get(), sizeof(value_t)*1024*1024)==0);
)
}

UNIT_TEST(TensorSum, dyadic,
    FullTensor res({2,2});
    FullTensor B({2});
    FullTensor C({2});

    Index i, J, K;
    
    B[{0}]=1;
    B[{1}]=2;
    
    C[{0}]=5;
    C[{1}]=9;
    
    FAILTEST(res(i,J) = B(i) + C(J));
//     TEST(compare_memory_to_vector(res.data.get(), {6,10,7,11}));
)

UNIT_TEST(TensorSum, threefold_sum,
    FullTensor res({2});
    FullTensor B({2});
    FullTensor C({2});
    FullTensor D({2});

    Index i, J, K;
    
    B[{0}]=1;
    B[{1}]=2;
    
    C[{0}]=5;
    C[{1}]=9;
    
    D[{0}]=7;
    D[{1}]=13;
    
    res(i) = B(i) + C(i) + D(i);
    TEST(compare_memory_to_vector(res.data.get(), {13,24}));
)
