/**
 * @file 
 * @short the source code to the "TTTensors" tutorial
 */

#include <xerus.h>

int main() {
	xerus::Index i,j,k;
	
	// the order of the tensors in the following calculation
	const size_t d = 5;
	
	// the dimensions of the states and the operators
	const std::vector<size_t> stateDims(d, 4);
	const std::vector<size_t> operatorDims(2*d, 4);
	
	// create random ttTensors of rank 2
    xerus::TTTensor B = xerus::TTTensor::random(stateDims, 2);
	xerus::TTTensor X = xerus::TTTensor::random(stateDims, 2);
	
	// and set the TTOperator A to be the identity
	xerus::TTOperator A = xerus::TTOperator::identity(operatorDims);

    // solve the system A*X = B with the ALS (this should give us X == B)
	xerus::ALS(A, X, B);
	
	// verify that X is indeed equal to B by calculating the Frobenius Norm of the difference
	std::cout << "Frobenius norm of X-B is: " << frob_norm(X-B) << " this should be almost 0..." << std::endl;
	
	// replace the operator by a random rank-2 operator
	A = xerus::TTOperator::random(operatorDims, 2);
	
	// ensure that A is symmetric by calculating @f$ A\cdot A^T @f$
	// here i^d signifies, that i should represent a multi-index of dimension d
	// the TTOperator of order 2d is thus fully indexed by two indices of the form i^d, j^d
	A(i/2, k/2) = A(i/2, j/2) * A(k/2, j/2);
	
	// the rank of A increased in the last operation:
	using xerus::misc::operator<<;
	std::cout << "The rank of A*A^T is " << A.ranks() << std::endl;
	
	// create a performance data object to keep track of the current algorithm progress (and print it to cout)
	xerus::PerformanceData perfData(true);
	
	// apply the ALS algorithm to the new system A*X=B and try to converge up to a relative error of @f$ 10^{-4} @f$
	xerus::ALS_SPD(A, X, B, 1e-4, perfData);
	
	// as the ALS will not modify the rank of X, the residual will most likely not be zero in the end
	// here i&n denotes that i should be a multiindex spanning all but n indices of the given tensor
	// in this case j&0 simply denotes that j should span all indices of X and B
	std::cout << "Residual ||A*X-B|| = " << frob_norm(A(i^d, j^d)*X(j&0) - B(i&0))/frob_norm(B) << " this is likely not equal to 0..." << std::endl;
}

/**
 * examplary output of this program:
 * > Frobenius norm of X-B is: 1.12371e-13 this should be almost 0...
 * > The rank of A*A^T is { 4, 4, 4, 4 }
 * > Residual ||A*X-B|| = 8170.59 this is likely not equal to 0...
 */


