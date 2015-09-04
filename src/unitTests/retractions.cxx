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


UNIT_TEST(algorithms, retractions,
	std::mt19937_64 rnd(0xC0CAC01A);
	std::normal_distribution<double> dist (0.0, 1.0);
	
	std::vector<size_t> stateDims(8,4);
	std::vector<size_t> stateRank(7,2);
	
	TTTensor X = TTTensor::random(stateDims, stateRank, rnd, dist);
	TTTensor zero = TTTensor::random(stateDims, stateRank, rnd, dist);
	zero *= 1e-16;
	
	TTTangentVector tangentZero(X, zero);
	for (size_t i=0; i<stateDims.size(); ++i) {
// 		LOG(ASDAS, frob_norm(tangentZero.components[i]));
		TEST(frob_norm(tangentZero.components[i]) < 1e-16 * (i==0?1000:1));
	}
	
	{
		TTTensor Y(X);
		ALSRetractionI(Y, tangentZero);
		TEST(frob_norm(X-Y) < 1e-8);
	}
	{
		TTTensor Y(X);
		ALSRetractionII(Y, zero);
		TEST(frob_norm(X-Y) < 1e-8);
	}
	{
		TTTensor Y(X);
		SubmanifoldRetractionI(Y, tangentZero);
		TEST(frob_norm(X-Y) < 1e-8);
	}
	{
		TTTensor Y(X);
		SubmanifoldRetractionII(Y, zero);
		TEST(frob_norm(X-Y) < 1e-8);
	}
	HOSVDRetraction hosvd(2);
	{
		TTTensor Y(X);
		hosvd(Y, tangentZero);
		TEST(frob_norm(X-Y) < 1e-8);
	}
	{
		TTTensor Y(X);
		hosvd(Y, zero);
		TEST(frob_norm(X-Y) < 1e-8);
	}
	
	TTTensor change = TTTensor::random(stateDims, stateRank, rnd, dist);
	const value_t EPS = 1e-4;
	change *= 1e-4 / frob_norm(change);
	
	Index i1,i2,i3,i4;
	TTTangentVector tangentChange(X, change);
	for (size_t i=0; i<stateDims.size(); ++i) {
		TEST(frob_norm(tangentChange.components[i]) > 1e-10);
		if (i!=0) {
			TEST(frob_norm(tangentChange.components[i](i1,i2,i3) * X.component(i)(i4,i2,i3)) < 1e-15);
		}
	}
	
	{
		TTTensor Y(X);
		ALSRetractionI(Y, tangentChange);
		value_t fnorm = frob_norm(X-Y);
		LOG(fnorm, fnorm);
		TEST(fnorm > 1e-8);
		TEST(fnorm < EPS);
	}
	{
		TTTensor Y(X);
		ALSRetractionII(Y, change);
		value_t fnorm = frob_norm(X-Y);
		LOG(fnorm, fnorm);
		TEST(fnorm > 1e-8);
		TEST(fnorm < EPS);
	}
	{
		TTTensor Y(X);
		SubmanifoldRetractionI(Y, tangentChange);
		value_t fnorm = frob_norm(X-Y);
		LOG(fnorm, fnorm);
		TEST(fnorm > 1e-8);
		TEST(fnorm < EPS);
	}
	{
		TTTensor Y(X);
		SubmanifoldRetractionII(Y, change);
		value_t fnorm = frob_norm(X-Y);
		LOG(fnorm, fnorm);
		TEST(fnorm > 1e-8);
		TEST(fnorm < EPS);
	}
	{
		TTTensor Y(X);
		hosvd(Y, tangentChange);
		value_t fnorm = frob_norm(X-Y);
		LOG(fnorm, fnorm);
		TEST(fnorm > 1e-8);
		TEST(fnorm < EPS);
	}
	{
		TTTensor Y(X);
		hosvd(Y, change);
		value_t fnorm = frob_norm(X-Y);
		LOG(fnorm, fnorm);
		TEST(fnorm > 1e-8);
		TEST(fnorm < EPS);
	}
)


UNIT_TEST(algorithms, vectorTransport,
	std::mt19937_64 rnd(0xC0CAC01A);
	std::normal_distribution<double> dist (0.0, 1.0);
	
	std::vector<size_t> stateDims(8,4);
	std::vector<size_t> stateRank(7,2);
	
	TTTensor X = TTTensor::random(stateDims, stateRank, rnd, dist);
	TTTensor change = TTTensor::random(stateDims, stateRank, rnd, dist);
	TTTangentVector tangentChange1(X, change);
	LOG(norms, frob_norm(change) << " " << frob_norm(tangentChange1.change_direction(X)));
	TEST(frob_norm(change) > frob_norm(tangentChange1.change_direction(X)));
	TTTangentVector tangentChange2(tangentChange1);
	
	tangentChange2 += tangentChange1 * (-1);
	for (size_t i=0; i<stateDims.size(); ++i) {
		TEST(frob_norm(tangentChange2.components[i]) < 1e-14);
	}
	
	tangentChange2 = TTTangentVector(X, tangentChange1.change_direction(X));
	LOG(asdasd, frob_norm(tangentChange2.change_direction(X) - tangentChange1.change_direction(X)));
	//ProjectiveVectorTransport(X, X, tangentChange2);
// 	tangentChange2 -= tangentChange1;
	for (size_t i=0; i<stateDims.size(); ++i) {
		LOG(aklsjd, frob_norm(tangentChange2.components[i]) / frob_norm(tangentChange1.components[i]));
		TEST(frob_norm(tangentChange2.components[i]) < 1e-14);
	}
)


