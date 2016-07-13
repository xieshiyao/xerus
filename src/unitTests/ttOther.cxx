#include<xerus.h>

#include "../../include/xerus/test/test.h"
using namespace xerus;

static misc::UnitTest tt_entryprod("TT", "entrywise_product", [](){
    Index i,j,k;
    
	TTTensor A = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2));
	TTTensor B = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2));
	
	Tensor Af(A);
	Tensor Bf(B);
	
	TTOperator Ao(Af);
	TTOperator Bo(Bf);
	
	TTTensor C = entrywise_product(A, B);
	TTOperator Co = entrywise_product(Ao, Bo);
	Tensor Cf = entrywise_product(Af, Bf);
	
	MTEST(frob_norm(Cf - Tensor(Co))/frob_norm(Cf) < 1e-13, frob_norm(Cf - Tensor(Co))/frob_norm(Cf));
	MTEST(frob_norm(Cf - Tensor(C))/frob_norm(Cf) < 1e-13, frob_norm(Cf - Tensor(C))/frob_norm(Cf));
	
	
	TTTensor D1 = entrywise_product(A, A);
	
	TTOperator Do1 = entrywise_product(Ao, Ao);
	
	Tensor Df = entrywise_product(Af, Af);

	MTEST(approx_equal(Df, Tensor(D1), 1e-13), frob_norm(Df-Tensor(D1)));
	MTEST(approx_equal(Df, Tensor(Do1), 1e-13), frob_norm(Df- Tensor(Do1)));
});

static misc::UnitTest tt_soft("TT", "soft_thresholding", [](){
    Index i,j,k;
    
	TTTensor A = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2));
	TTTensor B = TTTensor::random(std::vector<size_t>(10,2), std::vector<size_t>(9,2));
	
	
	Tensor Af(A);
	Tensor Bf(B);
	
	TTOperator Ao(Af);
	TTOperator Bo(Bf);
	
	TTTensor C = entrywise_product(A, B);
	TTOperator Co = entrywise_product(Ao, Bo);
	Tensor Cf = entrywise_product(Af, Bf);
	
	TEST(frob_norm(Cf - Tensor(Co))/frob_norm(Cf) < 1e-14);
	TEST(frob_norm(Cf - Tensor(C))/frob_norm(Cf) < 1e-14);
});
