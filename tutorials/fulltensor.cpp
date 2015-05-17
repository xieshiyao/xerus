/**
 * @file 
 * @short the source code to the "getting started" tutorial
 */

#include <xerus.h>

int main() {
	// constructing a simple 3x3 tensor (aka matrix)
    xerus::FullTensor A({3,3});

	// assigning individual entries in the tensor
    A[{0,0}] = 1.0;
    A[{1,1}] = 1.0;
    A[{2,2}] = 1.0;

    xerus::FullTensor B({3,3});
    xerus::FullTensor C({3,3});
    
	// simple tensor arithmetic: scalar multiplication and addition / subtraction
    B = 2*A;
    C = B-A;

    std::cout << "The Frobenius norm of A-C equals " << frob_norm(A-C) << std::endl; 

	// with c++11 random number distributions random tensors can be constructed
	// every entry will be sampled from the given distribution
    std::mt19937_64 rnd;
    std::normal_distribution<double> dist (0.0, 10.0);
    xerus::FullTensor X = xerus::FullTensor::construct_random({3,3,3}, rnd, dist);

	// reserving the variable names i,j and k to denote indices
    xerus::Index i,j,k,l,m,n;
    xerus::FullTensor D;
    
	// more complex tensor arithmetic can be written in an Einstein-like notation
	// here: the simple matrix product A*B
    D(i,k) = A(i,j)*B(j,k);

	// and the matrix product A*transpose(B)
    D(i,k) = A(i,j)*B(k,j);

	// similarly contractions of tensors of higher degree
    D(l,m,n) = X(i,j,k)*A(i,l)*B(j,m)*C(k,n);

	// after these contractions D should be equal to 2*X
	// we can verify this with a simple textual output of X and D
    std::cout << "X:" << std::endl << X.to_string() << std::endl;
    std::cout << "D:" << std::endl << D.to_string() << std::endl;
}
