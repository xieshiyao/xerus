#include<xerus.h>

#include "../../include/xerus/misc/test.h"
using namespace xerus;

UNIT_TEST(TT, entrywise_product,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 1.0);

    Index i,j,k;
    
	TTTensor A = TTTensor::construct_random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	TTTensor B = TTTensor::construct_random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	
	FullTensor Af(A);
	FullTensor Bf(B);
	
	TTOperator Ao(Af);
	TTOperator Bo(Bf);
	
	TTTensor C = TTTensor::entrywise_product(A, B);
	TTOperator Co = TTOperator::entrywise_product(Ao, Bo);
	FullTensor Cf = FullTensor::entrywise_product(Af, Bf);
	
	TEST(frob_norm(Cf - FullTensor(Co))/frob_norm(Cf) < 1e-14);
	TEST(frob_norm(Cf - FullTensor(C))/frob_norm(Cf) < 1e-14);
	
	
	TTTensor D1 = TTTensor::entrywise_product(A, A);
	TTTensor D2(A);
	D2.entrywise_square();
	
	
	TTOperator Do1 = TTOperator::entrywise_product(Ao, Ao);
	TTOperator Do2(Ao);
	Do2.entrywise_square();
	
	FullTensor Df = FullTensor::entrywise_product(Af, Af);

	TEST(approx_equal(Df, FullTensor(D1), 1e-4));
	TEST(approx_equal(Df, FullTensor(D2), 1e-4));
	TEST(approx_equal(Df, FullTensor(Do1), 1e-4));
	TEST(approx_equal(Df, FullTensor(Do2), 1e-4));
)

UNIT_TEST(TT, soft_thresholding,
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 1.0);

    Index i,j,k;
    
	TTTensor A = TTTensor::construct_random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	TTTensor B = TTTensor::construct_random(std::vector<size_t>(10,2), std::vector<size_t>(9,2), rnd, dist);
	
	
	FullTensor Af(A);
	FullTensor Bf(B);
	
	TTOperator Ao(Af);
	TTOperator Bo(Bf);
	
	TTTensor C = TTTensor::entrywise_product(A, B);
	TTOperator Co = TTOperator::entrywise_product(Ao, Bo);
	FullTensor Cf = FullTensor::entrywise_product(Af, Bf);
	
	TEST(frob_norm(Cf - FullTensor(Co))/frob_norm(Cf) < 1e-14);
	TEST(frob_norm(Cf - FullTensor(C))/frob_norm(Cf) < 1e-14);
)
