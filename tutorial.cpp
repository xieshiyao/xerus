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

// File containing the unit tests for the xerus library

#include "xerus.h"
    
int main() {
    xerus::FullTensor A({3,3});

    A[{0,0}] = 1.0;
    A[{1,1}] = 1.0;
    A[{2,2}] = 1.0;

    xerus::FullTensor B({3,3});
    xerus::FullTensor C({3,3});
    
    B = 2*A;
    C = B-A;

    std::cout << "The frobenius norm of A-C equals " << frob_norm(A-C) << std::endl; 

    std::mt19937_64 rnd;
    std::normal_distribution<double> dist (0.0, 10.0);
    xerus::FullTensor X = xerus::FullTensor::construct_random({3,3,3}, rnd, dist);

    xerus::Index i,j,k;
    xerus::FullTensor D;
    
    D(i,k) = A(i,j)*B(j,k);

    D(i,k) = A(i,j)*B(k,j);

    xerus::Index l,m,n;
    
    D(l,m,n) = X(i,j,k)*A(i,l)*B(j,m)*C(k,n);

    std::cout << "X:" << std::endl << X.to_string() << std::endl;
    
    std::cout << "D:" << std::endl << D.to_string() << std::endl;
}
