#include <xerus.h>

int main() {
	// construct the stiffness matrix A using a lambda function
	// and dividing it by h^2 = multiplying it with N^2
	xerus::Tensor A({512,512}, [](const std::vector<size_t> &idx){
		if (idx[0] == idx[1]) {
			return 2.0;
		} else if (idx[1] == idx[0]+1 || idx[1]+1 == idx[0]) {
			return -1.0;
		} else {
			return 0.0;
		}
	});
	
	A *= 512*512;
	
	// reinterpret the 512x512 tensor as a 2^18 tensor
	// and create (Q)TT decomposition of it
	A.reinterpret_dimensions(std::vector<size_t>(18, 2));
	xerus::TTOperator ttA(A);
	
	// and verify its rank
	using xerus::misc::operator<<;
	std::cout << "ttA ranks: " << ttA.ranks() << std::endl;
	
	// the right hand side of the equation both as Tensor and in (Q)TT format
	auto b = xerus::Tensor::ones(std::vector<size_t>(9, 2));
	auto ttb = xerus::TTTensor::ones(b.dimensions);
	
	// construct a random initial guess of rank 3 for the ALS algorithm
	xerus::TTTensor ttx = xerus::TTTensor::random(std::vector<size_t>(9, 2), std::vector<size_t>(8, 3));
	
	// and solve the system with the default ALS algorithm for symmetric positive operators
	xerus::ALS_SPD(ttA, ttx, ttb);
	
	// to perform arithmetic operations we need to define some indices
	xerus::Index i,j,k;
	
	// calculate the residual of the just solved system to evaluate its accuracy
	// here i^9 denotes a multiindex named i of dimension 9 (ie. spanning 9 indices of the respective tensors)
	double residual = frob_norm( ttA(i^9,j^9)*ttx(j^9) - ttb(i^9) );
	std::cout << "residual: " << residual << std::endl;
	
	// as an comparison solve the system exactly using the Tensor / operator
	xerus::Tensor x;
	x(j^9) = b(i^9) / A(i^9, j^9);
	
	// and calculate the Frobenius norm of the difference
	std::cout << "error: " << frob_norm(x - xerus::Tensor(ttx)) << std::endl;
}
