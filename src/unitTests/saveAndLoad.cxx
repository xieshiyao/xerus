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

UNIT_TEST2(SaveAndLoad, TensorTSV) {
	std::uniform_int_distribution<size_t> dimDist(1, 3);
	std::vector<size_t> dims1, dims2, dimsX, dimsA;
	
	for(size_t d = 1; d <= 8; ++d) {
		Tensor sA = Tensor::random(dimsA, misc::product(dimsA)/9+1, rnd, normalDist);
		Tensor sX = Tensor::random(dimsX, misc::product(dimsX)/9+1, rnd, normalDist);
		
		Tensor A = sA.dense_copy();
		Tensor X = sX.dense_copy();
		
		Tensor rA, rX;

		xerus::misc::exec("mkdir -p unitTestFiles");
		
		misc::save_to_file(A, "unitTestFiles/A.dat", misc::FileFormat::TSV);
		rA = misc::load_from_file<Tensor>("unitTestFiles/A.dat");
		MTEST(approx_equal(A, rA), frob_norm(A-rA));
		
		misc::save_to_file(X, "unitTestFiles/X.dat", misc::FileFormat::TSV);
		rX = misc::load_from_file<Tensor>("unitTestFiles/X.dat");
		MTEST(approx_equal(X, rX), frob_norm(X-rX));
		
		misc::save_to_file(sA, "unitTestFiles/sparseA.dat", misc::FileFormat::TSV);
		rA = misc::load_from_file<Tensor>("unitTestFiles/sparseA.dat");
		MTEST(approx_equal(sA, rA), frob_norm(sA-rA));
		
		misc::save_to_file(sX, "unitTestFiles/sparseX.dat", misc::FileFormat::TSV);
		rX = misc::load_from_file<Tensor>("unitTestFiles/sparseX.dat");
		MTEST(approx_equal(sX, rX), frob_norm(sX-rX));
		
		// Add a new dimension
		dims1.push_back(dimDist(rnd));
		dimsX = dims1;
		dimsA = dims1 | dims1;
	}
}});



UNIT_TEST2(SaveAndLoad, TensorBinary) {
	std::uniform_int_distribution<size_t> dimDist(1, 3);
	std::vector<size_t> dims1, dims2, dimsX, dimsA;
	
	for(size_t d = 1; d <= 8; ++d) {
		Tensor sA = Tensor::random(dimsA, misc::product(dimsA)/9+1, rnd, normalDist);
		Tensor sX = Tensor::random(dimsX, misc::product(dimsX)/9+1, rnd, normalDist);
		
		Tensor A = sA.dense_copy();
		Tensor X = sX.dense_copy();
		
		Tensor rA, rX;
		
		xerus::misc::exec("mkdir -p unitTestFiles");
		
		misc::save_to_file(A, "unitTestFiles/A_binary.dat", misc::FileFormat::BINARY);
		rA = misc::load_from_file<Tensor>("unitTestFiles/A_binary.dat");
		MTEST(approx_equal(A, rA), frob_norm(A-rA));
		
		misc::save_to_file(X, "unitTestFiles/X_binary.dat", misc::FileFormat::BINARY);
		rX = misc::load_from_file<Tensor>("unitTestFiles/X_binary.dat");
		MTEST(approx_equal(X, rX), frob_norm(X-rX));
		
		misc::save_to_file(sA, "unitTestFiles/sparseA_binary.dat", misc::FileFormat::BINARY);
		rA = misc::load_from_file<Tensor>("unitTestFiles/sparseA_binary.dat");
		MTEST(approx_equal(sA, rA), frob_norm(sA-rA));
		
		misc::save_to_file(sX, "unitTestFiles/sparseX_binary.dat", misc::FileFormat::BINARY);
		rX = misc::load_from_file<Tensor>("unitTestFiles/sparseX_binary.dat");
		MTEST(approx_equal(sX, rX), frob_norm(sX-rX));
		
		// Add a new dimension
		dims1.push_back(dimDist(rnd));
		dimsX = dims1;
		dimsA = dims1 | dims1;
	}
}});

