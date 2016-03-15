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


UNIT_TEST(Tensor, read_write_file,
	std::mt19937_64 rnd(0x77778888);
	std::normal_distribution<double> dist(0.0,1.0);
	Tensor A = Tensor::random({12,13,14}, rnd, dist);
	
	A.save_to_file("test.dat", FileFormat::TSV);
	Tensor Ab = Tensor::load_from_file("test.dat");
	MTEST(approx_equal(A, Ab), "dense tsv " << frob_norm(A-Ab));
	
	A.save_to_file("test.dat", FileFormat::BINARY);
	Ab = Tensor::load_from_file("test.dat");
	MTEST(approx_equal(A, Ab), "dense bin " << frob_norm(A-Ab));
	
	Tensor S({100,234,567}, Tensor::Representation::Sparse);
	S[{5,7,99}] = 1; S[{80,123,5}] = 6; S[{1,2,3}] = 4; S[{99,233,566}] = 5;
	S[{12,12,12}] = 8; S[{15,15,15}] = 7; S[{99,99,99}] = 3; S[{65,65,65}] = 2;
	
	S.save_to_file("test.dat", FileFormat::BINARY);
	Ab = Tensor::load_from_file("test.dat");
	MTEST(frob_norm(S - Ab) < 1e-16, "sparse bin " << frob_norm(S-Ab));
	
	S.save_to_file("test.dat", FileFormat::TSV);
	Ab = Tensor::load_from_file("test.dat");
	TEST(Ab.is_sparse());
	MTEST(frob_norm(S - Ab) < 1e-16, "sparse tsv " << frob_norm(S-Ab));
)

UNIT_TEST(TensorNetwork, read_write_file,
	std::mt19937_64 rnd(0x77778888);
	std::normal_distribution<double> dist(0.0,1.0);
	Tensor A = Tensor::random({12,13,14}, rnd, dist);
	Tensor B = Tensor::random({12,13,14}, rnd, dist);
	Tensor C = Tensor::random({12,13,14}, rnd, dist);
	Index i,j,k,l,m,n;
	TensorNetwork T;
	T(k,l,n) = A(i,j,k) * B(i,l,m) * C(n,j,m);
	T.save_to_file("test.dat", FileFormat::TSV);
	TensorNetwork Tb = TensorNetwork::load_from_file("test.dat");
	MTEST(T.dimensions == Tb.dimensions, T.dimensions << " vs " << Tb.dimensions);
	Tb.require_valid_network();
	MTEST(frob_norm(T(i&0)-Tb(i&0))/frob_norm(T) < 5e-16, frob_norm(T(i&0)-Tb(i&0))/frob_norm(T));
)
