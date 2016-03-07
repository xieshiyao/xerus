#include<xerus.h>

#include "../../include/xerus/misc/test.h"
using namespace xerus;

UNIT_TEST(TT, entrywise_product,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 1.0);

    Index i,j,k;
    
	TTTensor A = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	TTTensor B = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	
	Tensor Af(A);
	Tensor Bf(B);
	
	TTOperator Ao(Af);
	TTOperator Bo(Bf);
	
	TTTensor C = entrywise_product(A, B);
	TTOperator Co = entrywise_product(Ao, Bo);
	Tensor Cf = entrywise_product(Af, Bf);
	
	TEST(frob_norm(Cf - Tensor(Co))/frob_norm(Cf) < 1e-14);
	TEST(frob_norm(Cf - Tensor(C))/frob_norm(Cf) < 1e-14);
	
	
	TTTensor D1 = entrywise_product(A, A);
	TTTensor D2(A);
	D2.entrywise_square();
	
	
	TTOperator Do1 = entrywise_product(Ao, Ao);
	TTOperator Do2(Ao);
	Do2.entrywise_square();
	
	Tensor Df = entrywise_product(Af, Af);

	TEST(approx_equal(Df, Tensor(D1), 1e-14));
	TEST(approx_equal(Df, Tensor(D2), 1e-14));
	TEST(approx_equal(Df, Tensor(Do1), 1e-14));
	TEST(approx_equal(Df, Tensor(Do2), 1e-14));
)

UNIT_TEST(TT, soft_thresholding,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 1.0);

    Index i,j,k;
    
	TTTensor A = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	TTTensor B = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	
	
	Tensor Af(A);
	Tensor Bf(B);
	
	TTOperator Ao(Af);
	TTOperator Bo(Bf);
	
	TTTensor C = entrywise_product(A, B);
	TTOperator Co = entrywise_product(Ao, Bo);
	Tensor Cf = entrywise_product(Af, Bf);
	
	TEST(frob_norm(Cf - Tensor(Co))/frob_norm(Cf) < 1e-14);
	TEST(frob_norm(Cf - Tensor(C))/frob_norm(Cf) < 1e-14);
)
