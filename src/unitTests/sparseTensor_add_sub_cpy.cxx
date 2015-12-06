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

UNIT_TEST(SparseTensor, sum_matrix_2x2,
    SparseTensor res({2,2});
    SparseTensor B({2,2});
    SparseTensor C({2,2});

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
    TEST(res.compare_to_data({6,8,10,12}));
    res(i,J) = B(i,J) + C(J,i);
    TEST(res.compare_to_data({6,9,9,12}));
)
 
UNIT_TEST(SparseTensor, sum_lhs_equals_rhs,
    SparseTensor B({2,2});
    SparseTensor C({2,2});

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
	TEST(B.compare_to_data({6,8,10,12}));
    B(i,J) = B(i,J) + B(J,i);
    TEST(B.compare_to_data({12,18,18,24}));
)


UNIT_TEST(SparseTensor, sum_dyadic,
    SparseTensor res({2,2});
    SparseTensor B({2});
    SparseTensor C({2});

    Index i, J, K;
    
    FAILTEST(res(i,J) = B(i) + C(J));
)

UNIT_TEST(SparseTensor, sum_threefold,
    SparseTensor res({2});
    SparseTensor B({2});
    SparseTensor C({2});
    SparseTensor D({2});

    Index i, J, K;
    
    B[0]=1;
    B[1]=2;
    
    C[0]=5;
    C[1]=9;
    
    D[0]=7;
    D[1]=13;
    
    res(i) = B(i) + C(i) + D(i);
    TEST(res.compare_to_data({13,24}));
)
