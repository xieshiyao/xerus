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
		MTEST(frob_norm(tangentZero.components[i]) < 1e-14, frob_norm(tangentZero.components[i]));
	}
	
	{
		TTTensor Y(X);
		ALSRetractionI(Y, tangentZero);
		MTEST(frob_norm(X-Y) < 1e-8, "alsI " << frob_norm(X-Y));
	}
	{
		TTTensor Y(X);
		ALSRetractionII(Y, zero);
		MTEST(frob_norm(X-Y) < 1e-8, "alsII " << frob_norm(X-Y));
	}
// 	{
// 		TTTensor Y(X);
// 		SubmanifoldRetractionI(Y, tangentZero);
// 		MTEST(frob_norm(X-Y) < 1e-8, "submaniI " << frob_norm(X-Y));
// 	}
// 	{
// 		TTTensor Y(X);
// 		SubmanifoldRetractionII(Y, zero);
// 		MTEST(frob_norm(X-Y) < 1e-8, "submaniII " << frob_norm(X-Y));
// 	}
	HOSVDRetraction hosvd(2);
	{
		TTTensor Y(X);
		hosvd(Y, tangentZero);
		MTEST(frob_norm(X-Y) < 1e-8, "hosvdI " << frob_norm(X-Y));
	}
	{
		TTTensor Y(X);
		hosvd(Y, zero);
		MTEST(frob_norm(X-Y) < 1e-8, "hosvdII " << frob_norm(X-Y));
	}
	
	TTTensor change = TTTensor::random(stateDims, stateRank, rnd, dist);
	const value_t EPS = 1e-4;
	change *= EPS / frob_norm(change);
	
	Index i1,i2,i3,i4;
	TTTangentVector tangentChange(X, change);
	for (size_t i=0; i<stateDims.size(); ++i) {
		TEST(frob_norm(tangentChange.components[i]) > 1e-10);
		if (i!=0) {
			TEST(frob_norm(tangentChange.components[i](i1,i2,i3) * X.component(i)(i4,i2,i3)) < 1e-15);
		}
	}
	
	TTTensor riemannianChange = TTTensor(tangentChange);
	riemannianChange /= frob_norm(riemannianChange);
	
	{
		TTTensor Y(X);
		ALSRetractionI(Y, tangentChange);
		TTTensor actualChange = Y-X;
		value_t fnorm = frob_norm(actualChange);
		value_t angle = std::acos(value_t(riemannianChange(i1&0) * actualChange(i1&0)) / fnorm);
		MTEST(angle < 100 * EPS * EPS, "alsI " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm > EPS * 1e-4, "alsI " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm < EPS, "alsI " << fnorm << " " << angle << " " << EPS);
	}
	{
		TTTensor Y(X);
		ALSRetractionII(Y, change);
		TTTensor actualChange = Y-X;
		value_t fnorm = frob_norm(actualChange);
		value_t angle = std::acos(value_t(riemannianChange(i1&0) * actualChange(i1&0)) / fnorm);
		MTEST(angle < 100 * EPS * EPS, "alsII " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm > EPS * 1e-4, "alsII " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm < EPS, "alsII " << fnorm << " " << angle << " " << EPS);
	}
// 	{
// 		TTTensor Y(X);
// 		SubmanifoldRetractionI(Y, tangentChange);
// 		TTTensor actualChange = Y-X;
// 		value_t fnorm = frob_norm(actualChange);
// 		value_t angle = std::acos(value_t(riemannianChange(i1&0) * actualChange(i1&0)) / fnorm);
// 		MTEST(angle < 100 * EPS * EPS, "submaniI " << fnorm << " " << angle << " " << EPS);
// 		MTEST(fnorm > EPS * 1e-4, "submaniI " << fnorm << " " << angle << " " << EPS);
// 		MTEST(fnorm < EPS, "submaniI " << fnorm << " " << angle << " " << EPS);
// 	}
// 	{
// 		TTTensor Y(X);
// 		SubmanifoldRetractionII(Y, change);
// 		TTTensor actualChange = Y-X;
// 		value_t fnorm = frob_norm(actualChange);
// 		value_t angle = std::acos(value_t(riemannianChange(i1&0) * actualChange(i1&0)) / fnorm);
// 		MTEST(angle < 100 * EPS * EPS, "submaniII " << fnorm << " " << angle << " " << EPS);
// 		MTEST(fnorm > EPS * 1e-4, "submaniII " << fnorm << " " << angle << " " << EPS);
// 		MTEST(fnorm < EPS, "submaniII " << fnorm << " " << angle << " " << EPS);
// 	}
	{
		TTTensor Y(X);
		hosvd(Y, tangentChange);
		TTTensor actualChange = Y-X;
		value_t fnorm = frob_norm(actualChange);
		value_t angle = std::acos(value_t(riemannianChange(i1&0) * actualChange(i1&0)) / fnorm);
		MTEST(angle < 100 * EPS * EPS, "hosvdI " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm > EPS * 1e-4, "hosvdI " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm < EPS, "hosvdI " << fnorm << " " << angle << " " << EPS);
	}
	{
		TTTensor Y(X);
		hosvd(Y, change);
		TTTensor actualChange = Y-X;
		value_t fnorm = frob_norm(actualChange);
		value_t angle = std::acos(value_t(riemannianChange(i1&0) * actualChange(i1&0)) / fnorm);
		MTEST(angle < 100 * EPS * EPS, "hosvdII " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm > EPS * 1e-4, "hosvdII " << fnorm << " " << angle << " " << EPS);
		MTEST(fnorm < EPS, "hosvdII " << fnorm << " " << angle << " " << EPS);
	}
)


UNIT_TEST(algorithms, vectorTransport,
	std::mt19937_64 rnd(0xC0CAC01A);
	std::normal_distribution<double> dist (0.0, 1.0);
	
	std::vector<size_t> stateDims(8,4);
	std::vector<size_t> stateRank(7,2);
	Index j;
	
	TTTensor X = TTTensor::random(stateDims, stateRank, rnd, dist);
	TTTensor change = TTTensor::random(stateDims, stateRank, rnd, dist);
	TTTangentVector tangentChange1(X, change);
	// projection should decrease norm
	MTEST(frob_norm(change) > frob_norm(TTTensor(tangentChange1)), frob_norm(change) << " " << frob_norm(TTTensor(tangentChange1)));
	TTTangentVector tangentChange2(tangentChange1);
	
	// test whether copy construction, += and * work as intended
	tangentChange2 += tangentChange1 * (-1);
	for (size_t i=0; i<stateDims.size(); ++i) {
		TEST(frob_norm(tangentChange2.components[i]) < 1e-14);
	}
	
	// projection P*P = P
	tangentChange2 = TTTangentVector(X, TTTensor(tangentChange1));
	MTEST(frob_norm(TTTensor(tangentChange2) - TTTensor(tangentChange1)) < 1e-13, frob_norm(TTTensor(tangentChange2) - TTTensor(tangentChange1)));
	
	// vector transport from X to X should not change the vector at all
	ProjectiveVectorTransport(X, X, tangentChange2);
	tangentChange2 -= tangentChange1;
	for (size_t i=0; i<stateDims.size(); ++i) {
		MTEST(frob_norm(tangentChange2.components[i]) < 1e-13, frob_norm(tangentChange2.components[i]));
	}
	
	// tangent space of 10*X should be equal to tangent space of X
	tangentChange2 = TTTangentVector(10 * X, TTTensor(tangentChange1));
	MTEST(frob_norm(TTTensor(tangentChange2) - TTTensor(tangentChange1)) < 1e-13, frob_norm(TTTensor(tangentChange2) - TTTensor(tangentChange1)));
	
	value_t normOld = frob_norm(TTTensor(tangentChange1));
	for (value_t eps = 1e-10; eps <= 1e-4; eps *= 10) {
		TTTensor newX(X);
		SubmanifoldRetractionI(newX, tangentChange1 * eps);
		tangentChange2 = tangentChange1;
		ProjectiveVectorTransport(newX, X, tangentChange2);
		value_t normNew = frob_norm(TTTensor(tangentChange2));
		value_t changeAngle = (tangentChange1.scalar_product(tangentChange2) / normOld / normNew);
		MTEST(1-changeAngle < 1e-6, (1-changeAngle));
		MTEST(normNew <= normOld + 1e-9, normNew << " vs " << normOld);
	}
)


