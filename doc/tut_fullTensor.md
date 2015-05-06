# Getting started

In this part of the tutorial the basic functionalities of xerus for FullTensors are introduced. 

In xerus the class FullTensor is used to represent all dense tensors, i.e. tensors that are not sparse and not decomposed into a tensor network. 
Lets start of by creating a simple 3x3 matrix, i.e. a tensor of order two and dimensions `{3,3}`.
\code{.cpp}
#include "xerus.h"

int main() {
    xerus::FullTensor A({3,3});
\endcode
By default xerus zero-initializes every entry of a created tensor. 
The elements of a tensor object can be accesed via the `[]` operator. For example 
\code{.cpp}
    A[{0,0}] = 1.0;
    A[{1,1}] = 1.0;
    A[{2,2}] = 1.0;
\endcode
makes `A` the 3x3 identity matrix. 
Xerus FullTensors supports all basic arithmetic like assignemnt, addition, subtraction and scalar multiplication through the canonical operators.
\code{.cpp}
    xerus::FullTensor B({3,3});
    xerus::FullTensor C({3,3});
    
    B = 2*A;
    C = B-A;
\endcode
Now `C` should also be the 3x3 identity, which we can check using the frobenius norm of the difference
\code{.cpp}
    std::cout << "The frobenius norm of A-C equals " << frob_norm(A-C) << std::endl; 
\endcode

To create random constructed tensors via `C++11` random facilities, we can use
\code{.cpp}
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);
    xerus::FullTensor X = xerus::FullTensor::construct_random({3,3,3}, rnd, dist);
\endcode
which creates a `3x3x3` tensor with independend gaussian distributed random entries.


Xerus allows to perform contractions between Tensors using Einstein-like index conventions. To calculate the matrix matrix product of A and B we can simply write
\code{.cpp}
    xerus::Index i,j,k;
    xerus::FullTensor D;
    
    D(i,k) = A(i,j)*B(j,k);
\endcode
The syntax is however much more powerfull, for example the matrix product @f$ A B^T @f$ we can simply write 
\code{.cpp}
    D(i,k) = A(i,j)*B(k,j);
\endcode
or to calculate the contraction @f$ D_{l,m,n} = \sum_{i,j,k} X_{i,j,k} A_{i,l} B_{j,m} C_{k,n} @f$ we can write
\code{.cpp}
    xerus::Index l,m,n;
    
    D(i,k) = X(i,j,k)*A(i,l)*B(j,m)*C(k,n);
\endcode

As A and B are the identity and C is only @f$ 2I @f$ D should be two times X, which we can check by printing both X and D

\code{.cpp}
    std::cout << "X:" << std::endl << X.to_string() << std::endl;
    
    std::cout << "D:" << std::endl << D.to_string() << std::endl;
\endcode

The full source code can be found in `tutorial.cpp`

