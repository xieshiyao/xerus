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

UNIT_TEST(FullTensor, Contraction_Order_0,
    FullTensor A({});
    FullTensor B({});
    FullTensor res1({});
    
    A[{}] = 42;
    B[{}] = 73;
    
    contract(res1(), A(), B());
    TEST(compare_memory_to_vector(res1.data.get(), {42*73}));
)

UNIT_TEST(FullTensor, Contraction_Order_1,
    FullTensor A({2});
    FullTensor B({2});
    FullTensor C({3});
    FullTensor res1({2,2});
    FullTensor res2({});
    FullTensor res3({2,3});
    FullTensor res4({3,2});
    
    Index i,j;
    
    A[{0}] = 1;
    A[{1}] = 2;
    
    B[{0}] = 3;
    B[{1}] = 4;
    
    C[{0}] = 5;
    C[{1}] = 6;
    C[{2}] = 7;
    
    // Same dimensions
    // Contraction with no index being contracted
    contract(res1(i,j), A(i), B(j));
    TEST(compare_memory_to_vector(res1.data.get(), {3,4,6,8}));
    contract(res1(j,i), A(j), B(i));
    TEST(compare_memory_to_vector(res1.data.get(), {3,4,6,8}));
    contract(res1(i,j), A(j), B(i));
    TEST(compare_memory_to_vector(res1.data.get(), {3,6,4,8}));
    contract(res1(j,i), A(i), B(j));
    TEST(compare_memory_to_vector(res1.data.get(), {3,6,4,8}));
    
    // Contraction with one index being contracted
    contract(res2(), A(i), B(i));
    TEST(compare_memory_to_vector(res2.data.get(), {11}));
    
    // Different dimensions
    // Contraction with no index being contracted
    contract(res3(i,j), A(i), C(j));
    TEST(compare_memory_to_vector(res3.data.get(), {5,6,7,10,12,14}));
    contract(res3(j,i), A(j), C(i));
    TEST(compare_memory_to_vector(res3.data.get(), {5,6,7,10,12,14}));
    contract(res4(i,j), C(i), A(j));
    TEST(compare_memory_to_vector(res4.data.get(), {5,10,6,12,7,14}));
    contract(res4(j,i), C(j), A(i));
    TEST(compare_memory_to_vector(res4.data.get(), {5,10,6,12,7,14}));
    
    contract(res4(i,j), A(j), C(i));
    TEST(compare_memory_to_vector(res4.data.get(), {5,10,6,12,7,14}));
    contract(res4(j,i), A(i), C(j));
    TEST(compare_memory_to_vector(res4.data.get(), {5,10,6,12,7,14}));
    contract(res3(i,j), C(j), A(i));
    TEST(compare_memory_to_vector(res3.data.get(), {5,6,7,10,12,14}));
    contract(res3(j,i), C(i), A(j));
    TEST(compare_memory_to_vector(res3.data.get(), {5,6,7,10,12,14}));
)

UNIT_TEST(FullTensor, Contraction_Order_2_Same_Dimensions,
    FullTensor A({2,2});
    FullTensor B({2,2});
    FullTensor res1({2,2,2,2});
    FullTensor res2({2,2});
    FullTensor res3({});

    Index i,j,k,l;
    
    A[{0,0}]=1;
    A[{0,1}]=2;
    A[{1,0}]=3;
    A[{1,1}]=4;
    
    B[{0,0}]=5;
    B[{0,1}]=6;
    B[{1,0}]=7;
    B[{1,1}]=8;
    
    //Most possible contractions with no index being contracted
    //Switch pairs (i,j) and (k,l)
    contract(res1(i,j,k,l), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,7,8,10,12,14,16,15,18,21,24,20,24,28,32}));
    contract(res1(k,l,i,j), A(k,l), B(i,j));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,7,8,10,12,14,16,15,18,21,24,20,24,28,32}));
    contract(res1(i,j,k,l), A(k,l), B(i,j));
    TEST(compare_memory_to_vector(res1.data.get(), {5,10,15,20,6,12,18,24,7,14,21,28,8,16,24,32}));
    contract(res1(k,l,i,j), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,10,15,20,6,12,18,24,7,14,21,28,8,16,24,32}));
    
    //Switch i/j or k/l
    contract(res1(i,j,k,l), A(j,i), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,7,8,15,18,21,24,10,12,14,16,20,24,28,32}));
    contract(res1(j,i,k,l), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,7,8,15,18,21,24,10,12,14,16,20,24,28,32}));
    contract(res1(i,j,k,l), A(i,j), B(l,k));
    TEST(compare_memory_to_vector(res1.data.get(), {5,7,6,8,10,14,12,16,15,21,18,24,20,28,24,32}));
    contract(res1(i,j,l,k), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,7,6,8,10,14,12,16,15,21,18,24,20,28,24,32}));
    
    //Switch i/l or j/k
    contract(res1(i,j,k,l), A(l,j), B(k,i));
    TEST(compare_memory_to_vector(res1.data.get(), {5,15,7,21,10,20,14,28,6,18,8,24,12,24,16,32}));
    contract(res1(l,j,k,i), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,15,7,21,10,20,14,28,6,18,8,24,12,24,16,32}));
    contract(res1(i,j,k,l), A(i,k), B(j,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,10,12,7,8,14,16,15,18,20,24,21,24,28,32}));
    contract(res1(i,k,j,l), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,10,12,7,8,14,16,15,18,20,24,21,24,28,32}));
    
    //Switch i/k or j/l
    contract(res1(i,j,k,l), A(k,j), B(i,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,15,18,10,12,20,24,7,8,21,24,14,16,28,32}));
    contract(res1(k,j,i,l), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,15,18,10,12,20,24,7,8,21,24,14,16,28,32}));
    contract(res1(i,j,k,l), A(i,l), B(k,j));
    TEST(compare_memory_to_vector(res1.data.get(), {5,10,7,14,6,12,8,16,15,20,21,28,18,24,24,32}));
    contract(res1(i,l,k,j), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {5,10,7,14,6,12,8,16,15,20,21,28,18,24,24,32}));
    
    //Switch pairs using multi indices
    contract(res1(i^2, j^2), A(i^2), B(j^2));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,7,8,10,12,14,16,15,18,21,24,20,24,28,32}));
    contract(res1(j^2, i^2), A(j^2), B(i^2));
    TEST(compare_memory_to_vector(res1.data.get(), {5,6,7,8,10,12,14,16,15,18,21,24,20,24,28,32}));
    contract(res1(i^2, j^2), A(j^2), B(i^2));
    TEST(compare_memory_to_vector(res1.data.get(), {5,10,15,20,6,12,18,24,7,14,21,28,8,16,24,32}));
    contract(res1(j^2, i^2), A(i^2), B(j^2));
    TEST(compare_memory_to_vector(res1.data.get(), {5,10,15,20,6,12,18,24,7,14,21,28,8,16,24,32}));
    
    
    //All possible contractions with one index being contracted
    contract(res2(i, k), A(i, j), B(j, k));
    TEST(compare_memory_to_vector(res2.data.get(), {19,22,43,50}));
    contract(res2(i, k), A(j, i), B(j, k));
    TEST(compare_memory_to_vector(res2.data.get(), {26,30,38,44}));
    contract(res2(i, k), A(i, j), B(k, j));
    TEST(compare_memory_to_vector(res2.data.get(), {17,23,39,53}));
    contract(res2(i, k), A(j, i), B(k, j));
    TEST(compare_memory_to_vector(res2.data.get(), {23,31,34,46}));
    
    contract(res2(k, i), A(i, j), B(j, k));
    TEST(compare_memory_to_vector(res2.data.get(), {19,43,22,50}));
    contract(res2(k, i), A(j, i), B(j, k));
    TEST(compare_memory_to_vector(res2.data.get(), {26,38,30,44}));
    contract(res2(k, i), A(i, j), B(k, j));
    TEST(compare_memory_to_vector(res2.data.get(), {17,39,23,53}));
    contract(res2(k, i), A(j, i), B(k, j));
    TEST(compare_memory_to_vector(res2.data.get(), {23,34,31,46}));
    
    //All possible contractions with two indices being contracted
    contract(res3(), A(i, j), B(i, j));
    TEST(compare_memory_to_vector(res3.data.get(), {70}));
    contract(res3(), A(i, j), B(j, i));
    TEST(compare_memory_to_vector(res3.data.get(), {69}));
    contract(res3(), A(j, i), B(i, j));
    TEST(compare_memory_to_vector(res3.data.get(), {69}));
    contract(res3(), A(i^2), B(i^2));
    TEST(compare_memory_to_vector(res3.data.get(), {70}));
)

UNIT_TEST(FullTensor, Contraction_Order_2_Different_Dimensions,
    FullTensor A({1,2});
    FullTensor B({2,3});
    FullTensor res1({1,2,2,3});
    FullTensor res2({1,3});
    FullTensor res3({3,1});

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
    contract(res1(i,j,k,l), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {3,4,5,6,7,8,6,8,10,12,14,16}));
    contract(res1(i,j,k,l), A(i,k), B(j,l));
    TEST(compare_memory_to_vector(res1.data.get(), {3,4,5,6,8,10,6,7,8,12,14,16}));
    contract(res1(i,k,j,l), A(i,j), B(k,l));
    TEST(compare_memory_to_vector(res1.data.get(), {3,4,5,6,8,10,6,7,8,12,14,16}));
    
    //One Index contracted
    contract(res2(i,k), A(i,j), B(j,k));
    TEST(compare_memory_to_vector(res2.data.get(), {15,18,21}));
    contract(res3(k,i), A(i,j), B(j,k));
    TEST(compare_memory_to_vector(res2.data.get(), {15,18,21}));
    contract(res2(i,k), B(j,k), A(i,j));
    TEST(compare_memory_to_vector(res2.data.get(), {15,18,21}));
    contract(res3(k,i), B(j,k), A(i,j));
    TEST(compare_memory_to_vector(res2.data.get(), {15,18,21}));
)

UNIT_TEST(FullTensor, Contraction_Order_3_Same_Dimensions,
    FullTensor A({2,2,2});
    FullTensor B({2,2,2});
    FullTensor res2({2,2});
    FullTensor res3({});

    Index i,j,k,l,m;
    
    A[{0,0,0}]=1;
    A[{0,1,0}]=2;
    A[{1,0,0}]=3;
    A[{1,1,0}]=4;
    A[{0,0,1}]=5;
    A[{0,1,1}]=6;
    A[{1,0,1}]=7;
    A[{1,1,1}]=8;
    
    B[{0,0,0}]=5;
    B[{0,0,1}]=9;
    B[{0,1,0}]=6;
    B[{0,1,1}]=10;
    B[{1,0,0}]=7;
    B[{1,0,1}]=11;
    B[{1,1,0}]=8;
    B[{1,1,1}]=12;
    
    contract(res2(i,j), A(l,m,j), B(m,i,l));
    TEST(compare_memory_to_vector(res2.data.get(), {5+2*7+3*9+4*11, 5*5+6*7+7*9+8*11, 6+2*8+3*10+4*12, 5*6+6*8+7*10+8*12}));
    contract(res3(), A(i,j,k), B(i,k,j));
    TEST(approx_equal(res3[0], 5.0+5*6+2*9+6*10+3*7+7*8+4*11+8*12, 1e-13));
)

UNIT_TEST(FullTensor, Contraction_Order_3_Different_Dimensions,
    FullTensor A({1,2,3});
    FullTensor B({2,3,4});
    FullTensor C({2,1,3});
    FullTensor res1({1,2,3,2,3,4});
    FullTensor res2({1,3,4,3});
    FullTensor res3({1,3,3,4});
    FullTensor res4({1,4});
    FullTensor res5({1,1});
    FullTensor res6({1,3,3,1});
    
    FullTensor C1(4);

    Index i,j,k,l,m,n;
    
    A[{0,0,0}]=1;
    A[{0,0,1}]=2;
    A[{0,0,2}]=3;
    A[{0,1,0}]=4;
    A[{0,1,0}]=5; //NOTE
    A[{0,1,2}]=6;

    B[{0,0,0}]=7;
    B[{0,0,1}]=8;
    B[{0,0,2}]=9;
    B[{0,0,3}]=10;
    B[{0,1,0}]=11;
    B[{0,1,1}]=12;
    B[{0,1,2}]=13;
    B[{0,1,3}]=14;
    B[{0,2,0}]=15;
    B[{0,2,1}]=16;
    B[{0,2,2}]=17;
    B[{0,2,3}]=18;
    B[{1,0,0}]=19;
    B[{1,0,1}]=20;
    B[{1,0,2}]=21;
    B[{1,0,3}]=22;
    B[{1,1,0}]=23;
    B[{1,1,1}]=24;
    B[{1,1,2}]=25;
    B[{1,1,3}]=26;
    B[{1,2,0}]=27;
    B[{1,2,1}]=28;
    B[{1,2,2}]=29;
    B[{1,2,3}]=30;
    
    C[{0,0,0}]=31;
    C[{0,0,1}]=32;
    C[{0,0,2}]=33;
    C[{1,0,0}]=34;
    C[{1,0,1}]=35;
    C[{1,0,2}]=36;
    
    //WARNING These results are obtained by an older version of the library and may contain errors themself!
    contract(res1(i,j,k,l,m,n), A(i,j,k), B(l,m,n));
    TEST(compare_memory_to_vector(res1.data.get(), {7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,42,48,54,60,66,72,78,84,90,96,102,108,114,120,126,132,138,144,150,156,162,168,174,180}));
    
    contract(res1(i,l,m,j,k,n), A(i,j,k), B(l,m,n));
    TEST(compare_memory_to_vector(res1.data.get(), {7,8,9,10,14,16,18,20,21,24,27,30,35,40,45,50,0,0,0,0,42,48,54,60,11,12,13,14,22,24,26,28,33,36,39,42,55,60,65,70,0,0,0,0,66,72,78,84,15,16,17,18,30,32,34,36,45,48,51,54,75,80,85,90,0,0,0,0,90,96,102,108,19,20,21,22,38,40,42,44,57,60,63,66,95,100,105,110,0,0,0,0,114,120,126,132,23,24,25,26,46,48,50,52,69,72,75,78,115,120,125,130,0,0,0,0,138,144,150,156,27,28,29,30,54,56,58,60,81,84,87,90,135,140,145,150,0,0,0,0,162,168,174,180}));
    
    contract(res2(i,m,n,k), A(i,j,k), B(j,m,n));
    TEST(compare_memory_to_vector(res2.data.get(),{102,14,135,108,16,144,114,18,153,120,20,162,126,22,171,132,24,180,138,26,189,144,28,198,150,30,207,156,32,216,162,34,225,168,36,234}));
    
    contract(res3(i,k,m,n), A(i,j,k), B(j,m,n));
    TEST(compare_memory_to_vector(res3.data.get(),{102,108,114,120,126,132,138,144,150,156,162,168,14,16,18,20,22,24,26,28,30,32,34,36,135,144,153,162,171,180,189,198,207,216,225,234}));
    
    contract(res4(i,n), A(i,j,k), B(j,k,n));
    TEST(compare_memory_to_vector(res4.data.get(),{331,348,365,382}));
    contract(res4(i,n), A(i,j^2), B(j^2,n));
    TEST(compare_memory_to_vector(res4.data.get(),{331,348,365,382}));
    
    contract(res5(i,n), A(i,j,k), C(j,n,k));
    TEST(compare_memory_to_vector(res5.data.get(),{580}));
    
    contract(res6(i,j,m,n), A(i,k,m), C(k,n,j));
    TEST(compare_memory_to_vector(res6.data.get(),{201,62,297,207,64,306,213,66,315}));
)
