# The Tensor Class

The basic building stone of this library and all Tensor Network methods is the Tensor, represented in the class `xerus::Tensor`.
To simplify the work with these objects, `xerus` contains a number of helper functions that allow the quick creation and 
modification of sparse and dense tensors. In the following we will list the most important ones but advise you to also read
the tutorials on [indices and equations](\ref md_indices) and [decompositions]() to have the full toolset with which to work on individual 
tensors.

## Creation of Tensors
The most basic tensors can be created with the empty constructor
~~~.cpp
// creates a degree 0 tensor
A = xerus::Tensor();
~~~
~~~.py
# creates a degree 0 tensor
A = xerus.Tensor()
~~~
it is of degree 0 and represents the single number 0. Similarly the constructors that take either the degree or a vector of 
dimensions as input create (sparse) tensors that are equal to 0 everywhere
~~~.cpp
// creates a 1x1x1 tensor with entry 0
B = xerus::Tensor(3);
// creates a sparse 2x2x2 tensor without any entries
C = xerus::Tensor({2,2,2});
~~~
~~~.py
# creates a 1x1x1 tensor with entry 0
B = xerus.Tensor(3)
# creates a sparse 2x2x2 tensor without any entries
C = xerus.Tensor([2,2,2])
# equivalently: xerus.Tensor(dim=[2,2,2])
~~~
The latter of these can be forced to create a dense tensor instead which can either be initialized to 0 or uninitialized
~~~.cpp
// creates a dense 2x2x2 tensor with all entries set to 0
D = xerus::Tensor({2,2,2}, xerus::Tensor::Representation::Dense);
// creates a dense 2x2x2 tensor with uninitialized entries
E = xerus::Tensor({2,2,2}, xerus::Tensor::Representation::Dense, xerus::Tensor::Initialisation::None);
~~~
~~~.py
# creates a dense 2x2x2 tensor with all entries set to 0
D = xerus.Tensor(dim=[2,2,2], repr=xerus.Tensor.Representation.Dense)
# creates a dense 2x2x2 tensor with uninitialized entries
E = xerus.Tensor(dim=[2,2,2], repr=xerus.Tensor.Representation.Dense, init=xerus.Tensor.Initialisation.None)
~~~

Other commonly used tensors (apart from the 0 tensor) are available through named constructors:
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
// a 4x4x4 sparse tensor with 10 random entries in uniformly distributed random positions
xerus::Tensor::random({4,4,4}, 10);
// a (4x4) x (4x4) random orthogonal operator drawn according to the Haar measure
xerus::Tensor::random_orthogonal({4,4},{4,4});
~~~
~~~.py
# a 2x3x4 tensor with all entries = 1
xerus.Tensor.ones([2,3,4])
# an (3x4) x (3x4) identity operator
xerus.Tensor.identity([3,4,3,4])
# a 3x4x3x4 tensor with superdiagonal = 1 (where all 4 indices coincide) and = 0 otherwise
xerus.Tensor.kronecker([3,4,3,4])
# a 2x2x2 tensor with a 1 in position {1,1,1} and 0 everywhere else
xerus.Tensor.dirac([2,2,2], [1,1,1])
# equivalently xerus.Tensor.dirac(dim=[2,2,2], pos=[1,1,1])

# a 4x4x4 tensor with i.i.d. Gaussian random values
xerus.Tensor.random([4,4,4])
# a 4x4x4 sparse tensor with 10 random entries in uniformly distributed random positions
xerus.Tensor.random([4,4,4], n=10)
# a (4x4) x (4x4) random orthogonal operator drawn according to the Haar measure
xerus.Tensor.random_orthogonal([4,4],[4,4])
~~~

If the entries of the tensor should be calculated externally, it is possible in c++ to either pass the raw data directly (as
`std::unique_ptr<double>` or `std::shared_ptr<double>`, check section 'Advanced Use and Ownership of Data' for the latter!) 
or use a callback / lambda function to populate the entries:
~~~.cpp
std::unique_ptr<double> ptr = foo();
// transfer ownership of the data to the Tensor object of size 2x2x2
// NOTE: make sure that the dimensions are correct as there is no way for xerus to check this!
F = xerus::Tensor({2,2,2}, ptr);

// create a dense 3x3x3 tensor with every entry populated by a callback (lambda) function
G = xerus::Tensor({3,3,3}, [](const std::vector<size_t> &_idx) -> double {
	// somehow derive the entry from the index positions in _idx
	double result = double( _idx[0] * _idx[1] * _idx[2] );
	return result;
});

// create a sparse 16x16x16 tensor with 16 entries determined by a callback (lambda) function
H = xerus::Tensor({16,16,16}, 16, [](size_t num, size_t max) -> std::pair<size_t, double> {
	// insert number 1/5, 2/5, 3/5, 4/5, 5/5
	// at positions 0 (= {0,0,0}), 17 (= {0,1,1}), 34 (= {0,2,2}) ... respectively
	return std::pair<size_t,double>(num*17, double(num)/double(max));
});
~~~
~~~.py
# Transfering ownership of raw data directly is not possible from within python.

# create a dense 3x3x3 tensor with every entry populated by a callback (lambda) function
G = xerus.Tensor.from_function([3,3,3], lambda idx: idx[0]*idx[1]*idx[2])

# Creating a sparse tensor from a callback function in python is not supported by xerus.
# This behaviour can easily be achieved by setting the correspoding values as seen below though.
~~~

In python raw data structures are not directly compatible to those used in `xerus` internally. Tensors can be constructed from 
`numpy.ndarray` objects though. This function will also implicitely accept pythons native array objects.
~~~.py
# convert a numpy tensor (identity matrix) to a xerus.Tensor object
T = xerus.Tensor.from_ndarray(numpy.eye(2))
# alternatively the function also accepts pythons native arrays
U = xerus.Tensor.from_ndarray([[1,0], [0,1]])
~~~

Last but not least it is possible to populate the entries of a tensor by explicitely accessing them. During this process, an 
initially sparse 0 tensor as it is created by the default constructors will automatically be converted to a dense object as soon
as `xerus` deems this to be preferable.
~~~.cpp
// creating an identity matrix by explicitely setting non-zero entries
V = xerus::Tensor({2,2});
V[{0,0}] = 1.0; // equivalently: V[0] = 1.0;
V[{1,1}] = 1.0; // equivalently: V[3] = 1.0;
~~~
~~~.py
# creating an identity matrix by explicitely setting non-zero entries
V = xerus.Tensor([2,2])
V[[0,0]] = 1.0 # equivalently: V[0] = 1.0
V[[1,1]] = 1.0 # equivalently: V[3] = 1.0
~~~
Explicitely constructing a tensor similarly in a blockwise fashion will be covered in the tutorial on [indices and equations](\ref md_indices)


## Sparse and Dense Representations
The last example already mentioned, that `xerus` will dynamically convert sparse tensors do their dense counterpart when this
seems to be preferable. For this it does not matter whether the number of entries increased due to explicit access, summation or
contraction of tensors or as the result of decompositions.

This behaviour can be modified by changing the global setting
~~~.cpp
// tell xerus to convert sparse tensors to dense if 1 in 4 entries are non-zero
xerus::Tensor::sparsityFactor = 4;
~~~
~~~.py
# tell xerus to convert sparse tensors to dense if 1 in 4 entries are non-zero
xerus.Tensor.sparsityFactor = 4
~~~
in particular, setting the [sparsityFactor](\ref xerus::Tensor::sparsityFactor) to 0 will disable this feature.
~~~.cpp
// stop xerus from automatically converting sparse tensors to dense
xerus::Tensor::sparsityFactor = 0;
~~~
~~~.py
# stop xerus from automatically converting sparse tensors to dense
xerus.Tensor.sparsityFactor = 0
~~~
Note though, that calculations with non-sparse Tensors that are stored in a sparse representation are typically much slower than
in dense representation. You should thus manually convert overly full sparse Tensors to the dense representation.

To do this there are a number of ways to interact with the representation of `xerus::Tensor` objects. Above we already saw, that
the constructors can be used to explicitely construct sparse (default behaviour) or dense tensors. For already existing objects
you can use the member functions [.is_sparse()](\ref xerus::Tensor::is_sparse()) and [.is_dense()](\ref xerus::Tensor::is_dense()) to query their representation. To change representations call the 
member functions [.use_dense_representation()](\ref xerus::Tensor::use_dense_representation()) or [.use_sparse_representation()](\ref xerus::Tensor::use_sparse_representation()) to change it inplace or [.dense_copy()](\ref xerus::Tensor::dense_copy()) or 
[.sparse_copy()](\ref xerus::Tensor::sparse_copy()) to obtain new tensor objects with the desired representation.

To make more informed decisions about whether a conversion might be useful the tensor objects can be queried for the number of
defined entries with [.sparsity()](\ref xerus::Tensor::sparsity()) or for the number of non-zero entries with [.count_non_zero_entries()](\ref xerus::Tensor::count_non_zero_entries()) ().
~~~.cpp
// create a sparse tensor with 100 random entries
W = xerus::Tensor::random({100,100}, 100);
// query its sparsity. likely output: "100 100"
std::cout << W.sparsity() << ' ' << W.count_non_zero_entries() << std:endl;

// store an explicit 0 value in the sparse representation
W[{0,0}] = 0.0;
// query its sparsity. likely output: "101 100"
std::cout << W.sparsity() << ' ' << W.count_non_zero_entries() << std:endl;

// convert the tensor to dense representation
W.use_dense_representation();
// query its sparsity. likely output: "10000 100"
std::cout << W.sparsity() << ' ' << W.count_non_zero_entries() << std:endl;
~~~
~~~.py
# create a sparse tensor with 100 random entries
W = xerus.Tensor.random(dim=[100,100], n=100)
# query its sparsity. likely output: "100 100"
print(W.sparsity(), W.count_non_zero_entries())

# store an explicit 0 value in the sparse representation
W[[0,0]] = 0.0
# query its sparsity. likely output: "101 100"
print(W.sparsity(), W.count_non_zero_entries())

# convert the tensor to dense representation
W.use_dense_representation()
# query its sparsity. likely output: "10000 100"
print(W.sparsity(), W.count_non_zero_entries())
~~~


## Operators and Modifications
We have already seen the most basic method of modifying a tensor via the [operator[]](\ref xerus::Tensor::operator[]()). With it
and the index notation presented in the [indices and equations](\ref md_indices) tutorial, most desired manipulations can be 
represented. Some of them would still be cumbersome though, so `xerus` includes several helper functions to make your life easier.
The purpose of this section is to present the most important ones.

Naturally tensor objects can be used in arithmetic expressions whereever this is well defined: addition and subtraction of equally
sized tensors as well as multiplication and division by scalars. Note, that tensors in general do not have the cannonical
isomorphism to operators that matrices have. As such the matrix multiplication can not be generalized trivially to tensors. Instead
you will have to use indexed equations to express such contractions in xerus.
~~~.cpp
xerus::Tensor A = xerus::Tensor::random({100});
// normalizing the vector
A = A / frob_norm(A);
// adding a small pertubation
xerus::Tensor B = A + 1e-5 * xerus::Tensor::random({100});
// normalizing B
B /= frob_norm(B);
// determining the change from A to B
B -= A;
std::cout << "Distance on the unit sphere: " << B.frob_norm() << std::endl;
~~~
~~~.py
A = xerus.Tensor.random([100])
# normalizing the vector
A = A / xerus.frob_norm(A)
# adding a small pertubation
B = A + 1e-5 * xerus.Tensor.random([100])
# normalizing B
B /= B.frob_norm()
# determining the change from A to B
B -= A
print("Distance on the unit sphere: ", B.frob_norm())
~~~
Where we have already used [.frob_norm()](\ref xerus::Tensor::frob_norm()) to calculate the Frobenius norm of a Tensor.

Many theoretical results use flattenings or expansions (vectorizations or tensorizations) to reduce tensor problems to ones with
vectors and matrices. While this is not as common in practical applications it can still be useful at times. Because the data
that the computer uses internally does not need not be reordered, such a reinterpretation of the tensor objects are constant in
time.
~~~.cpp
A = xerus::Tensor::identity({2,2});
// flatten the identity matrix to the vector (1, 0, 0, 1)
A.reinterpret_dimensions({4});
~~~
~~~.py
A = xerus.Tensor.identity([2,2])
# flatten the identity matrix to the vector (1, 0, 0, 1)
A.reinterpret_dimensions([4])
~~~
This operation is obviously only possible when the total number of entries remains unchanged.

If you want to change the dimensions of a tensor such that the total size changes, you have to specify how to do this. `xerus`
provides three functions to help you in such a case: [.resize_mode()](\ref xerus::Tensor::resize_mode()) changes the dimension
of a single mode by adding zero slates or removing existing slates at a given position; [.fix_mode()](\ref xerus::Tensor::fix_mode())
reduces the tensor to an object of degree d-1 that corresponds to the slate, selected in the call to the function; finally
[.remove_slate()](\ref xerus::Tensor::remove_slate()) is a simplified version of `.resize_mode()` that removes a single slate
from the tensor, reducing the dimension of the specified mode by one.
~~~.cpp
// start with a 2x2 matrix filled with ones
A = xerus::Tensor::ones({2,2});
// insert another row (consisting of zeros) i.e. increase dimension of mode 0 to 3
A.resize_mode(0, 3, 1);
// select the first column i.e. reduce to fixed value 0 for mode 1
A.fix_mode(1, 0);
// expected output: "1.0 0.0 1.0"
std::cout << A.to_string() << std::endl;
~~~
~~~.py
# start with a 2x2 matrix filled with ones
A = xerus.Tensor.ones([2,2])
# insert another row (consisting of zeros) i.e. increase dimension of mode 0 to 3
A.resize_mode(mode=0, newDim=3, cutPos=1)
# select the first column i.e. reduce to fixed value 0 for mode 1
A.fix_mode(mode=1, value=0)
# expected output: "1.0 0.0 1.0"
print(A)
~~~

At the moment the Hadamard product is not available in a indexed notation (due to a lack of overloadable operators). Its
behaviour can instead be achieved with [entrywise_product()](\ref xerus::entrywise_product()).
~~~.cpp
// constructing a tensor with i.i.d. entries sampled from a Chi-Squared distribution by (Note: this is not the most efficient way to achieve this!)
// 1. constructing a tensor with Gaussian i.i.d. entries
A = xerus::Tensor::random({10,10,10,10});
// 2. performing an entrywise product of this tensor with itself
A = entrywise_product(A, A);
~~~
~~~.py
# constructing a tensor with i.i.d. entries sampled from a Chi-Squared distribution by (Note: this is not the most efficient way to achieve this!)
# 1. constructing a tensor with Gaussian i.i.d. entries
A = xerus.Tensor.random([10,10,10,10])
# 2. performing an entrywise product of this tensor with itself
A = xerus.entrywise_product(A, A)
~~~


## Output and Storing
In the above examples we have already seen two kind of queries to the `Tensor` objects: [.to_string()](\ref xerus::Tensor::to_string()) 
to obtain a human readable string representation of the Tensor and [.frob_norm()](\ref xerus::Tensor::frob_norm()) or equivalently
[frob_norm(Tensor))(\ref xerus::frob_norm(Tensor)) to obtain the Frobenius norm of a Tensor. Note for the former, that it is 
meant purely for debugging purposes, in particular of smaller objects, and it is not adequately possible to reconstruct the original
Tensor from this output.

Storing Tensors to files such that they can be reconstructed exactly from those is instead possible with [save_to_file()](\ref xerus::misc::save_to_file())
and respectively [load_from_file()](\ref xerus::misc::load_from_file()).
~~~.cpp

~~~



## Advanced Use and Ownership of Data
shared_ptr, factor, rest c++ only: accessing internal representations
