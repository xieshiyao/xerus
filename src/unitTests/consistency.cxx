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


UNIT_TEST2(Consistency, sum_and_difference) {
	std::uniform_int_distribution<size_t> dimDist(1, 3);
	
	std::vector<size_t> dims1, dims2, dims3, dimsX, dimsY, dimsA, dimsB;
	
	const Index i, j, k;
	
	for(size_t d = 1; d <= 8; ++d) {
		Tensor A = Tensor::random(dimsA, rnd, normalDist);
		Tensor B = Tensor::random(dimsB, rnd, normalDist);
		Tensor X = Tensor::random(dimsX, rnd, normalDist);
		Tensor Y = Tensor::random(dimsX, rnd, normalDist);
		Tensor C;
		
		TTOperator ttA(A, 0.2);
		TTOperator ttB(B, 0.2);
		TTOperator ttoC;
		TTTensor ttX(X, 0.2);
		TTTensor ttY(Y, 0.2);
		TTTensor tttC;
		
		TensorNetwork tnA(ttA);
		TensorNetwork tnB(ttB);
		TensorNetwork tnX(ttX);
		TensorNetwork tnY(ttY);
		TensorNetwork tnC;
		
		A = ttA;
		B = ttB;
		X = ttX;
		Y = ttY;
		
		Tensor sA = A.sparse_copy();
		Tensor sB = B.sparse_copy();
		Tensor sX = X.sparse_copy();
		Tensor sY = Y.sparse_copy();
		Tensor sC;
		
		TEST(approx_equal(A, sA));
		TEST(approx_equal(A, tnA));
		TEST(approx_equal(A, ttA));
		TEST(approx_equal(sA, tnA));
		TEST(approx_equal(sA, ttA));
		TEST(approx_equal(tnA, ttA));
		
		TEST(approx_equal(B, sB));
		TEST(approx_equal(B, tnB));
		TEST(approx_equal(B, ttB));
		TEST(approx_equal(sB, tnB));
		TEST(approx_equal(sB, ttB));
		TEST(approx_equal(tnB, ttB));
		
		TEST(approx_equal(X, sX));
		TEST(approx_equal(X, tnX));
		TEST(approx_equal(X, ttX));
		TEST(approx_equal(sX, tnX));
		TEST(approx_equal(sX, ttX));
		TEST(approx_equal(tnX, ttX));
		
		TEST(approx_equal(Y, sY));
		TEST(approx_equal(Y, tnY));
		TEST(approx_equal(Y, ttY));
		TEST(approx_equal(sY, tnY));
		TEST(approx_equal(sY, ttY));
		TEST(approx_equal(tnY, ttY));
		
		
		C = X+X;
		sC = sX+sX;
		tttC = ttX+ttX;
// 		tnC = tnX+tnX; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, tttC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, tttC, 1e-14));
// 		TEST(approx_equal(tnC, tttC));
		
		C = X+Y;
		sC = sX+sY;
		tttC = ttX+ttY;
// 		tnC = tnX+tnY; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, tttC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, tttC, 1e-14));
// 		TEST(approx_equal(tnC, tttC));
		
		C = X-Y;
		sC = sX-sY;
		tttC = ttX-ttY;
// 		tnC = tnX-tnY; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, tttC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, tttC, 1e-14));
// 		TEST(approx_equal(tnC, tttC));
		
		C = X+Y+X;
		sC = sX+sY+sX;
		tttC = ttX+ttY+ttX;
// 		tnC = tnX+tnY+tnX; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, tttC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, tttC, 1e-14));
// 		TEST(approx_equal(tnC, tttC));
		
		
		C = 3.7*X+Y+X-3*Y;
		sC = 3.7*sX+sY+sX-3*sY;
		tttC = 3.7*ttX+ttY+ttX-3*ttY;
// 		tnC = 3.7*tnX+tnY+tnX-3*tnY; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, tttC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, tttC, 1e-14));
// 		TEST(approx_equal(tnC, tttC));
		
		
		C = A+A;
		sC = sA+sA;
		ttoC = ttA+ttA;
// 		tnC = tnA+tnA; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttoC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttoC, 1e-14));
// 		TEST(approx_equal(tnC, ttoC));
		
		C = A+B;
		sC = sA+sB;
		ttoC = ttA+ttB;
// 		tnC = tnA+tnB; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttoC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttoC, 1e-14));
// 		TEST(approx_equal(tnC, ttoC));
		
		
		C = A-B;
		sC = sA-sB;
		ttoC = ttA-ttB;
// 		tnC = tnA-tnB; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttoC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttoC, 1e-14));
// 		TEST(approx_equal(tnC, ttoC));
		
		
		C = 3.7*A+B-1.2*A;
		sC = 3.7*sA+sB-1.2*sA;
		ttoC = 3.7*ttA+ttB-1.2*ttA;
// 		tnC = 3.7*tnA+tnB-1.2*tnA; TODO
		
		TEST(approx_equal(C, sC));
// 		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttoC, 1e-14));
// 		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttoC, 1e-14));
// 		TEST(approx_equal(tnC, ttoC));
		
		
		// Add a new dimension
		dims1.push_back(dimDist(rnd));
		dimsX = dims1;
		dimsY = dims1;
		dimsA = dims1 | dims1;
		dimsB = dims1 | dims1;
	}
}});


UNIT_TEST2(Consistency, operator_times_tensor) {
	std::uniform_int_distribution<size_t> dimDist(1, 4);
	
	std::vector<size_t> dims1, dims2, dims3, dimsX, dimsY, dimsA, dimsB;
	
	const Index i, j, k;
	
	for(size_t d = 1; d <= 8; ++d) {
		Tensor A = Tensor::random(dimsA, rnd, normalDist);
		Tensor B = Tensor::random(dimsB, rnd, normalDist);
		Tensor X = Tensor::random(dimsX, rnd, normalDist);
		Tensor Y = Tensor::random(dimsX, rnd, normalDist);
		Tensor C;
		
		TTOperator ttA(A, 0.2); 
		TTOperator ttB(B, 0.2); 
		TTTensor ttX(X, 0.2);
		TTTensor ttY(Y, 0.2);
		TTTensor ttC;
		
		TensorNetwork tnA(ttA);
		TensorNetwork tnB(ttB);
		TensorNetwork tnX(ttX);
		TensorNetwork tnY(ttY);
		TensorNetwork tnC;
		
		A = ttA;
		B = ttB;
		X = ttX;
		Y = ttY;
		
		Tensor sA = A.sparse_copy();
		Tensor sB = B.sparse_copy();
		Tensor sX = X.sparse_copy();
		Tensor sY = Y.sparse_copy();
		Tensor sC;
		
		TEST(approx_equal(A, sA));
		TEST(approx_equal(A, tnA));
		TEST(approx_equal(A, ttA));
		TEST(approx_equal(sA, tnA));
		TEST(approx_equal(sA, ttA));
		TEST(approx_equal(tnA, ttA));
		
		TEST(approx_equal(B, sB));
		TEST(approx_equal(B, tnB));
		TEST(approx_equal(B, ttB));
		TEST(approx_equal(sB, tnB));
		TEST(approx_equal(sB, ttB));
		TEST(approx_equal(tnB, ttB));
		
		TEST(approx_equal(X, sX));
		TEST(approx_equal(X, tnX));
		TEST(approx_equal(X, ttX));
		TEST(approx_equal(sX, tnX));
		TEST(approx_equal(sX, ttX));
		TEST(approx_equal(tnX, ttX));
		
		TEST(approx_equal(Y, sY));
		TEST(approx_equal(Y, tnY));
		TEST(approx_equal(Y, ttY));
		TEST(approx_equal(sY, tnY));
		TEST(approx_equal(sY, ttY));
		TEST(approx_equal(tnY, ttY));
		
		
		C(i&0) = A(i/2,j/2)*X(j&0);
		sC(i&0) = sA(i/2,j/2)*sX(j&0);
		ttC(i&0) = ttA(i/2,j/2)*ttX(j&0);
		tnC(i&0) = tnA(i/2,j/2)*tnX(j&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		C(i&0) = B(i/2,j/2)*X(j&0);
		sC(i&0) = sB(i/2,j/2)*sX(j&0);
		ttC(i&0) = ttB(i/2,j/2)*ttX(j&0);
		tnC(i&0) = tnB(i/2,j/2)*tnX(j&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		C(i&0) = B(i/2,j/2)*Y(i&0);
		sC(i&0) = sB(i/2,j/2)*sY(i&0);
		ttC(i&0) = ttB(i/2,j/2)*ttY(i&0);
		tnC(i&0) = tnB(i/2,j/2)*tnY(i&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		C(i&0) = B(j/2,i/2)*A(i/2,j/2)*X(j&0);
		sC(i&0) = sB(j/2,i/2)*sA(i/2,j/2)*sX(j&0);
		ttC(i&0) = ttB(j/2,i/2)*A(i/2,j/2)*ttX(j&0);
		tnC(i&0) = tnB(j/2,i/2)*A(i/2,j/2)*tnX(j&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		// Add a new dimension
		dims1.push_back(dimDist(rnd));
		dims2.push_back(dimDist(rnd));
		dimsX = dims1;
		dimsY = dims2;
		dimsA = dims1 | dims1;
		dimsB = dims2 | dims1;
	}
}});

UNIT_TEST2(Consistency, fix_slate) {
	std::uniform_int_distribution<size_t> dimDist(1, 4);
	
	std::vector<size_t> dims1, dims2, dims3, dimsX, dimsY, dimsA, dimsB;
	
	const Index i, j, k;
	
	for(size_t d = 1; d <= 8; ++d) {
		// Add a new dimension
		dims1.push_back(dimDist(rnd));
		dims2.push_back(dimDist(rnd));
		dimsX = dims1;
		dimsY = dims2;
		dimsA = dims1 | dims1;
		dimsB = dims2 | dims1;
		
		Tensor A = Tensor::random(dimsA, rnd, normalDist);
		Tensor B = Tensor::random(dimsB, rnd, normalDist);
		Tensor X = Tensor::random(dimsX, rnd, normalDist);
		Tensor Y = Tensor::random(dimsX, rnd, normalDist);
		Tensor C;
		
		TTOperator ttA(A, 0.2); 
		TTOperator ttB(B, 0.2); 
		TTTensor ttX(X, 0.2);
		TTTensor ttY(Y, 0.2);
		TTTensor ttC;
		
		TensorNetwork tnA(ttA);
		TensorNetwork tnB(ttB);
		TensorNetwork tnX(ttX);
		TensorNetwork tnY(ttY);
		TensorNetwork tnC;
		
		A = ttA;
		B = ttB;
		X = ttX;
		Y = ttY;
		
		Tensor sA = A.sparse_copy();
		Tensor sB = B.sparse_copy();
		Tensor sX = X.sparse_copy();
		Tensor sY = Y.sparse_copy();
		Tensor sC;
		
		
		std::uniform_int_distribution<size_t> slateSelect(0, d-1);
		const size_t slate = slateSelect(rnd);
		std::uniform_int_distribution<size_t> posSelect(0, dimsX[slate]-1);
		const size_t position = posSelect(rnd);
		
		A.fix_slate(slate+d, position);
		A.fix_slate(slate, position);
		B.fix_slate(slate+d, position);
		B.fix_slate(slate, position);
		X.fix_slate(slate, position);
		Y.fix_slate(slate, position);
		
		sA.fix_slate(slate+d, position);
		sA.fix_slate(slate, position);
		sB.fix_slate(slate+d, position);
		sB.fix_slate(slate, position);
		sX.fix_slate(slate, position);
		sY.fix_slate(slate, position);
		
		ttA = TTOperator(A); // Can't use fix_slate
		ttB = TTOperator(B); // Can't use fix_slate
		ttX = TTTensor(X);
		ttY = TTTensor(Y);
		
		ttA.require_correct_format();
		ttB.require_correct_format();
		ttX.require_correct_format();
		ttY.require_correct_format();
		
		tnA.fix_slate(slate+d, position);
		tnA.fix_slate(slate, position);
		tnB.fix_slate(slate+d, position);
		tnB.fix_slate(slate, position);
		tnX.fix_slate(slate, position);
		tnY.fix_slate(slate, position);
		
		TEST(approx_equal(A, sA));
		TEST(approx_equal(A, tnA));
		TEST(approx_equal(A, ttA));
		TEST(approx_equal(sA, tnA));
		TEST(approx_equal(sA, ttA));
		TEST(approx_equal(tnA, ttA));
		
		TEST(approx_equal(B, sB));
		TEST(approx_equal(B, tnB));
		TEST(approx_equal(B, ttB));
		TEST(approx_equal(sB, tnB));
		TEST(approx_equal(sB, ttB));
		TEST(approx_equal(tnB, ttB));
		
		TEST(approx_equal(X, sX));
		TEST(approx_equal(X, tnX));
		TEST(approx_equal(X, ttX));
		TEST(approx_equal(sX, tnX));
		TEST(approx_equal(sX, ttX));
		TEST(approx_equal(tnX, ttX));
		
		TEST(approx_equal(Y, sY));
		TEST(approx_equal(Y, tnY));
		TEST(approx_equal(Y, ttY));
		TEST(approx_equal(sY, tnY));
		TEST(approx_equal(sY, ttY));
		TEST(approx_equal(tnY, ttY));
		
		
		C(i&0) = A(i/2,j/2)*X(j&0);
		sC(i&0) = sA(i/2,j/2)*sX(j&0);
		ttC(i&0) = ttA(i/2,j/2)*ttX(j&0);
		tnC(i&0) = tnA(i/2,j/2)*tnX(j&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		C(i&0) = B(i/2,j/2)*X(j&0);
		sC(i&0) = sB(i/2,j/2)*sX(j&0);
		ttC(i&0) = ttB(i/2,j/2)*ttX(j&0);
		tnC(i&0) = tnB(i/2,j/2)*tnX(j&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		C(i&0) = B(i/2,j/2)*Y(i&0);
		sC(i&0) = sB(i/2,j/2)*sY(i&0);
		ttC(i&0) = ttB(i/2,j/2)*ttY(i&0);
		tnC(i&0) = tnB(i/2,j/2)*tnY(i&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
		
		C(i&0) = B(j/2,i/2)*A(i/2,j/2)*X(j&0);
		sC(i&0) = sB(j/2,i/2)*sA(i/2,j/2)*sX(j&0);
		ttC(i&0) = ttB(j/2,i/2)*A(i/2,j/2)*ttX(j&0);
		tnC(i&0) = tnB(j/2,i/2)*A(i/2,j/2)*tnX(j&0);
		
		TEST(approx_equal(C, sC));
		TEST(approx_equal(C, tnC));
		TEST(approx_equal(C, ttC));
		TEST(approx_equal(sC, tnC));
		TEST(approx_equal(sC, ttC));
		TEST(approx_equal(tnC, ttC));
	}
}});
