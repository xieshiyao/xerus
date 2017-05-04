# The Tensor Class

The basic building stone of this library and all Tensor Network methods is the Tensor, represented in the class `xerus::Tensor`.
To simplify the work with these objects, `xerus` contains a number of helper functions that allow the quick creation and 
modification of sparse and dense tensors. In the following we will list the most important ones but advise you to also read
the tutorials on [indices and equations]() and [decompositions]() to have the full toolset with which to work on individual 
tensors.

## Creation of Tensors
The most basic tensors can be created with the empty constructor
~~~.cpp
A = xerus::Tensor()
~~~
it is of degree 0 and represents the single number 0. Similarly the constructors that take either the degree or a vector of 
dimensions as input create (sparse) tensors that are equal to 0 everywhere
~~~.cpp
// creates a 1x1x1 tensor with entry 0
B = xerus::Tensor(3);
// creates a sparse 2x2x2 tensor without any entries
C = xerus::Tensor({2,2,2});
~~~
The latter of these can be forced to create a dense tensor instead which can either be initialized to 0 or uninitialized
~~~.cpp
// creates a dense 2x2x2 tensor with all entries set to 0
D = xerus::Tensor({2,2,2}, xerus::Tensor::Representation::Dense);
// creates a dense 2x2x2 tensor with uninitialized entries
E = xerus::Tensor({2,2,2}, xerus::Tensor::Representation::Dense, xerus::Tensor::Initialisation::None);
~~~

Other common tensors (apart from the 0 tensor) are available through named constructors:
~~~.cpp
// a 2x3x4 tensor with all entries = 1
xerus::Tensor::ones({2,3,4});
// an (3x4) x (3x4) identity operator
xerus::Tensor::identity({3,4,3,4});
// a 3x4x3x4 tensor with superdiagonal = 1 (where all 4 indices coincide) and = 0 otherwise
xerus::Tensor::kronecker({3,4,3,4});
// a 2x2x2 tensor with a 1 in position {1,1,1} and 0 everywhere else
xerus::Tensor::dirac({2,2,2}, {1,1,1});

// a 4x4x4 tensor with i.i.d. Gaussian random values
xerus::Tensor::random({4,4,4});
// a (4x4) x (4x4) random orthogonal operator drawn according to the Haar measure
xerus::Tensor::random_orthogonal({4,4},{4,4});
// a 4x4x4 sparse tensor with 10 random entries in uniformly distributed random positions
xerus::Tensor::random({4,4,4}, 10);
~~~

If the entries of the tensor should be calculated externally, it is possible to either pass the raw data directly (as
`std::unique_ptr<double>` or `std::shared_ptr<double>`) or use a callback / lambda function to populate the entries:
~~~.cpp
std::unique_ptr<double> ptr = foo();
// transfer ownership of the data to the Tensor object of size 2x2x2
// NOTE: make sure that the dimensions are correct as there is no way for xerus to check this!
F = xerus::Tensor({2,2,2}, ptr);

// create a dense 2x2x2 tensor with every entry populated by a callback (lambda) function
G = xerus::Tensor({2,2,2}, [](const std::vector<size_t> &_idx) -> double {
	// somehow derive the entry from the index positions in _idx
	double result = double( _idx[0] * _idx[1] * _idx[2] );
	return result;
});

// create a sparse 10x10x10 tensor with 5 entries determined by a callback (lambda) function
H = xerus::Tensor({10,10,10}, 5, [](size_t num, size_t max) -> std::pair<size_t, double> {
	// insert number 1/5, 2/5, 3/5, 4/5, 5/5
	// at positions 0 (= {0,0,0}), 17 (= {0,1,7}), 34 (= {0,3,4}) ... respectively
	return std::pair<size_t,double>(num*17, double(num)/double(max));
});
~~~

## Used Representation

## Operators and Modifications

## Advanced Use and Ownership of Data
