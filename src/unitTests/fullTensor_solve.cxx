// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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

UNIT_TEST(Tensor, solve_Ax_equals_b,
    Index i,j,k;
          
    Tensor A1({4,2,2});
    A1[{0,0,0}] = 1;
    A1[{1,0,1}] = 1;
    A1[{2,1,0}] = 1;
    A1[{3,1,1}] = 1;
    
    Tensor b1({4});
    b1[0] = 73;
    b1[1] = -73;
    b1[2] = 128;
    b1[3] = 93;
    
    Tensor x1({2,2});
    
    
    x1(k,i) = (b1(j)/A1(j,k,i));
    
    TEST((b1[0] - x1[{0,0}]) < 1e-14);
    TEST((b1[1] - x1[{0,1}]) < 1e-14);
    TEST((b1[2] - x1[{1,0}]) < 1e-14);
    TEST((b1[3] - x1[{1,1}]) < 1e-14);
    
    Tensor A2({4,2,2});
    A2[{0,0,0}] = 1;
    A2[{1,0,1}] = 1;
    A2[{2,1,0}] = 0;
    A2[{3,1,1}] = 0;
    
    Tensor b2({4});
    b2[0] = 73;
    b2[1] = -73;
    b2[2] = 0;
    b2[3] = 0;
    
    Tensor x2({2,2});
    
    x2(k,i) = (b2(j)/A2(j,k,i));
    
    TEST((b2[0] - x2[{0,0}]) < 1e-14);
    TEST((b2[1] - x2[{0,1}]) < 1e-14);
    TEST((x2[{1,0}]) < 1e-14);
    TEST((x2[{1,1}]) < 1e-14);
)

UNIT_TEST(Tensor, solve_sparse,
	std::mt19937_64 rnd(0x5EED);
	std::normal_distribution<double> dist(0.0, 1.0);
	const size_t N = 100;
	std::uniform_int_distribution<size_t> eDist(1,N*N-1);
	
	Index i,j,k;
	
	Tensor id = Tensor::identity({N,N});
	Tensor r = Tensor({N}, [](size_t _i)->value_t{return double(_i);});
	Tensor x;
	x(i) = r(j) / id(j,i);
	MTEST(frob_norm(x-r) < 1e-14, "d " << frob_norm(x-r));
	
	r.use_sparse_representation();
	x(i) = r(j) / id(j,i);
	MTEST(frob_norm(x-r) < 1e-14, "d " << frob_norm(x-r));
	r.use_dense_representation();
	
	// consistency with dense solve:
	for (size_t n=0; n<N*3; ++n) {
		id[eDist(rnd)] = dist(rnd);
	}
	id.use_sparse_representation();
	
	// test faithful reconstruction
	internal::CholmodSparse idt(id.get_unsanitized_sparse_data(), N, N, false);
	Tensor id2({N,N});
	id2.get_unsanitized_sparse_data() = idt.to_map();
	MTEST(frob_norm(id-id2) < 1e-15, frob_norm(id-id2)); 
	
	Tensor fid(id);
	fid.use_dense_representation();
	Tensor fx;
// 	id.use_sparse_representation();
	TEST(id.is_sparse());
	
	fx(i) = r(j) / fid(j,i);
	x(i) = r(j) / id(j,i);
	MTEST(frob_norm(fx-x)/frob_norm(x)<3e-14, frob_norm(fx-x)/frob_norm(x));
)
