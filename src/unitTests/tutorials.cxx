#include<xerus.h>

#include "../../include/xerus/misc/test.h"
using namespace xerus;



UNIT_TEST(Tutorials, quick_start,
	xerus::FullTensor A({512,512}, [](const std::vector<size_t> &idx){
		if (idx[0] == idx[1]) {
			return 2.0;
		} else if (idx[1] == idx[0]+1 || idx[1]+1 == idx[0]) {
			return -1.0;
		} else {
			return 0.0;
		}
	});
	
	A *= 512*512;
	
	A.reinterpret_dimensions(std::vector<size_t>(18, 2));
	xerus::TTOperator ttA(A);
	
// 	std::cout << "ttA ranks: " << ttA.ranks() << std::endl;
// 	TEST(ttA.ranks()
	
	xerus::FullTensor b({512}, []() {
		return 1.0;
	});
	
	b.reinterpret_dimensions(std::vector<size_t>(9, 2));
	xerus::TTTensor ttb(b);
	
	std::mt19937_64 rnd(0x5EED);
	std::normal_distribution<double> dist (0.0, 1.0);
	xerus::TTTensor ttx = xerus::TTTensor::random(std::vector<size_t>(9, 2), std::vector<size_t>(8, 3), rnd, dist);
	
	xerus::ALS(ttA, ttx, ttb);
	
	xerus::Index i,j,k;
	
	double residual = frob_norm( ttA(i^9,j^9)*ttx(j^9) - ttb(i^9) );
	TEST(residual < 1e-9);
// 	std::cout << "residual: " << residual << std::endl;
	
	xerus::FullTensor x;
	x(j^9) = b(i^9) / A(i^9, j^9);
	TEST(frob_norm(x(i&0) - ttx(i&0)) < 1e-9);
// 	std::cout << FullTensor(ttx).to_string() << std::endl;
// 	std::cout << "error: " << frob_norm(x(i&0) - ttx(i&0)) << std::endl;
)
