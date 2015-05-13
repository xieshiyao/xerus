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


UNIT_TEST(ALS, identity,
    //Random numbers
    std::mt19937_64 rnd;
    rnd.seed(73);
	std::normal_distribution<value_t> dist (0.0, 1.0);
    
    Index k,l,m,n,o,p;
    
    FullTensor X({10, 10, 10});
    FullTensor B = FullTensor::construct_random({10, 10, 10}, rnd, dist);
    
    FullTensor I({10,10,10,10,10,10}, [](const std::vector<size_t> &_idx) {
		if (_idx[0]==_idx[3] && _idx[1] == _idx[4] && _idx[2] == _idx[5]) {
			return 1.0;
		} else {
			return 0.0;
		}
	});
    
    X(k^3) = I(k^3,l^3)*B(l^3);
    TEST(frob_norm(X(k^3) - B(k^3)) < 1e-13);
    
    TTTensor ttB(B, 0.001);
    TTTensor ttX(X, 0.001);
    TTOperator ttI(I, 0.001);
	
	ttX(k^3) = ttI(k^3,l^3)*ttB(l^3);
    TEST(frob_norm(ttI(k^3,l^3)*ttB(l^3) - ttB(k^3)) < 1e-13);
	TEST(frob_norm(ttI(k^3,l^3)*ttX(l^3) - ttX(k^3)) < 1e-13);
    
    std::vector<double> perfdata;
    
    TEST(ALS(ttI, ttX, ttB, 0.001, &perfdata) < 0.01);
    TEST(frob_norm(FullTensor(ttX)(k^3) - FullTensor(ttB)(k^3)) < 1e-13 * 1000);
	LOG(unit_test, "perf: " << perfdata);
    perfdata.clear();
	
    ttX = TTTensor::construct_random(ttX.dimensions, ttX.ranks(), rnd, dist);
    TEST(ALS(ttI, ttX, ttB, 0.001, &perfdata) < 0.01);
	LOG(unit_test, "perf: " << perfdata);
	LOG(unit_test, "norm: " << frob_norm(FullTensor(ttX)(k^3) - FullTensor(ttB)(k^3)));
    TEST(frob_norm(FullTensor(ttX)(k^3) - FullTensor(ttB)(k^3)) < 1e-9); // approx 1e-16 * dim * max_entry
    /*
    LOG(fsda, "Perform");
    for(FullTensor& part : ttX.representation) { memset(part.data.get(), 0, part.dataSize*sizeof(double)); }
    performALS(ttI, ttX, ttB, 0.001);
    LOG(fsda, "TEST");
    TEST(frob_norm(FullTensor(ttX)(k^3) - FullTensor(ttB)(k^3)) < 1e-13);
*/

)


// UNIT_TEST(ALS, random,
// 	std::mt19937_64 rnd;
//     rnd.seed(73);
// 	std::normal_distribution<value_t> dist (0.0, 1.0);
//     
//     Index i,j,k,l,m,n,o,p;
//     
// 	FullTensor A = FullTensor::construct_random({256, 256}, rnd, dist);
// 	A(k,l) = A(m,k)*A(m,l);
//     FullTensor B = FullTensor::construct_random({256}, rnd, dist);
//     TTVector X1(B);
// 	std::vector<double> perf;
// 	{
// 		TTMatrix ttA(A);
// 		TTVector ttB(B);
// 		TEST(ALS(ttA, X1, ttB, 1e-8, &perf) < 1e-6);
// 		LOG(unit_test, "perf " << perf);
// 	}
// 	FullTensor resFull(1);
//     resFull(i) = A(i,j)*X1.representation[0](0,j,0)-B(i);
//     LOG(sdf, "True: " << frob_norm(resFull(i)));
//     FullTensor tmpX(X1);
// 	A.reinterpret_dimensions(std::vector<uint>(8,4));
// 	B.reinterpret_dimensions(std::vector<uint>(4,4));
//     tmpX.reinterpret_dimensions(std::vector<uint>(4,4));
//     FullTensor resFuller(4);
//     resFuller(i^4) = A(i^4,j^4)*tmpX(j^4)-B(i^4);
//     LOG(sdf, "Reimpertrete: " << frob_norm(resFuller(i^4)));
//     
// 	TTVector X2(tmpX, -1e-27);//(B);
//     TTMatrix ttA(A, -1e-27);
//     TTVector ttB(B, -1e-27);
//     resFuller(i^4) = FullTensor(ttA)(i^4,j^4)*FullTensor(X2)(j^4)-FullTensor(ttB)(i^4);
//     LOG(sdf, "TT and back: " << frob_norm(resFuller(i^4)));
//     
// 	{	
// 		perf.clear();
//         TTVector res(4);
//         res(i^4) = ttA(i^4, j^4)*X2(j^4) - ttB(i^4);
//         LOG(sdf, "TT: " << frob_norm(res(i^4)));
//         LOG(sdf, "TT (anderes frob norm): " << frob_norm(FullTensor(res)(i^4)));
//         res.round(1e-13);
//         LOG(sdf, "TT (noch andere frob norm): " << frob_norm(res(i^4)));
// 		TEST(performALS(ttA, X2, ttB, 1e-8, &perf));
// 		LOG(unit_test, "perf " << perf);
// 	}
// 	FullTensor fX(X2);
// 	fX.reinterpret_dimensions({256});
// 	double fnorm = frob_norm(FullTensor(X1)(k&0)-fX(k&0));
// 	LOG(unit_test, "norm: " << fnorm);
// 	TEST(fnorm < 1e-4);
// )
