---
layout: post
title: "QTT Decomposition"
date: 1000-12-10
topic: "Examples"
section: "Examples"
---
__tabsInit
# QTT Decomposition

This guide can be used as a quick start into `xerus`. It will introduce some basic functionality of the library,
demonstrate the general layout and is enough for very basic applications. It is recommended to also have a look
at the more detailed guides for all classes one wishes to use though - or even have a look at the doxygen class documentation for details on all functions.

It is assumed that you have already obtained and compiled the library itself as well as know how to link against it.
If this is not the case, please refer to the [building xerus](/building_xerus) page.

In the following we will solve a FEM equation arising from the heat equation using a QTT decomposition and the ALS algorithm.

To start we create the stiffness matrix as a dense (ie. not sparse or decomposed) tensor.
To define the entries we pass a function to the constructor of the `Tensor` object that will be
called for every entry with a vector of size_t integers defining the indices of the current entry.
As it is simpler to think of the stiffness matrix in its original form rather than the QTT form we will 
set the dimensions to 512x512.

__tabsStart
~~~ cpp
xerus::Tensor A({512,512}, [](const std::vector<size_t> &idx){
	if (idx[0] == idx[1]) {
		return 2.0;
	} else if (idx[1] == idx[0]+1 || idx[1]+1 == idx[0]) {
		return -1.0;
	} else {
		return 0.0;
	}
});
~~~
__tabsMid
~~~ python
def A_fill(idx):
	if idx[0] == idx[1] :
		return 2
	elif idx[1] == idx[0]+1 or idx[1]+1 == idx[0] :
		return -1
	else :
		return 0

A = xerus.Tensor([512,512], A_fill)
~~~
__tabsEnd

To account for the $ h^2 $ factor that we have ignored so far we simply multipy the operator by $ N^2 $.

__tabsStart
~~~ cpp
A *= 512*512;
~~~
__tabsMid
~~~ python
A *= 512*512
~~~
__tabsEnd

By reinterpreting the dimension and thus effectively treating the tensor as a $ 2^{18} $ instead of a $ 512^2 $ tensor,
the decomposition into a `TTTensor` will give us the stiffness matrix in a QTT format.

__tabsStart
~~~ cpp
A.reinterpret_dimensions(std::vector<size_t>(18, 2));
xerus::TTOperator ttA(A);
~~~
__tabsMid
~~~ python
A.reinterpret_dimensions([2,]*18)
ttA = xerus.TTOperator(A)
~~~
__tabsEnd

As the Laplace operator is representable exactly in a low-rank QTT format, the rank of this `ttA` should be low after this construction.
We can verify this by printing the ranks:

__tabsStart
~~~ cpp
using xerus::misc::operator<<; // to be able to pipe vectors
std::cout << "ttA ranks: " << ttA.ranks() << std::endl;
~~~
__tabsMid
~~~ python
print("ttA ranks:", ttA.ranks())
~~~
__tabsEnd

For the right-hand-side we will take a simple tensor that is equal to 1 at every position. As this is a commonly used tensor,
we can simply use the named constructor provide dby `xerus`.

__tabsStart
~~~ cpp
auto b = xerus::Tensor::ones(std::vector<size_t>(9, 2));
auto ttb = xerus::TTTensor::ones(b.dimensions);
~~~
__tabsMid
~~~ python
b = xerus.Tensor.ones([2,]*9)
ttb = xerus.TTTensor.ones(b.dimensions)
~~~
__tabsEnd

To have an initial vector for the ALS algorithm we create a random TTTensor of the desired rank 
(3 in this case - note, that this is the exact rank of the solution).

__tabsStart
~~~ cpp
xerus::TTTensor ttx = xerus::TTTensor::random(std::vector<size_t>(9, 2), std::vector<size_t>(8, 3));
~~~
__tabsMid
~~~ python
ttx = xerus.TTTensor.random([2,]*9, [3,]*8)
~~~
__tabsEnd

With these three tensors (the operator `ttA`, the right-hand-side `ttb` and the initial guess `ttx`)
we can now perform the ALS algorithm to solve for `ttx` (note here, that the _SPD suffix chooses the variant of the ALS
that assumes that the given operator is symmetric positive definite)

__tabsStart
~~~ cpp
xerus::ALS_SPD(ttA, ttx, ttb);
~~~
__tabsMid
~~~ python
xerus.ALS_SPD(ttA, ttx, ttb)
~~~
__tabsEnd

To verify the calculation performed by the ALS we will need to perform some arithmetic operations. 
As these require the definition of (relative) index orderings in the tensors, we define some indices

__tabsStart
~~~ cpp
xerus::Index i,j,k;
~~~
__tabsMid
~~~ python
i,j,k = xerus.indices(3)
~~~
__tabsEnd

and use these in calculations like `A(i,j)*x(j) - b(i)`. Note though, that our tensors are of a higher
degree due to the QTT decomposition and we thus need to specify the corresponding dimension of the
multiindices i,j, and k with eg. `i^9` to denote a multiindex of dimension 9.

__tabsStart
~~~ cpp
double residual = frob_norm( ttA(i^9,j^9)*ttx(j^9) - ttb(i^9) );
std::cout << "residual: " << residual << std::endl;
~~~
__tabsMid
~~~ python
residual = xerus.frob_norm( ttA(i^9,j^9)*ttx(j^9) - ttb(i^9) )
print("residual:", residual)
~~~
__tabsEnd

To get the actual error of the ALS solution (rather than just the residual) we can calculate the exact solution
of the system using the FullTensors A, x and b and the / operator

__tabsStart
~~~ cpp
xerus::Tensor x;
x(j^9) = b(i^9) / A(i^9, j^9);
~~~
__tabsMid
~~~ python
x = xerus.Tensor()
x(j^9) << b(i^9) / A(i^9, j^9)
~~~
__tabsEnd

In the comparison of this exact solution `x` and the ALS solution `ttx`, we have to decide whether we want to cast
the TTTensor to a Tensor or vice versa to be able to subtract them.

__tabsStart
~~~ cpp
std::cout << "error: " << frob_norm(x - xerus::Tensor(ttx)) << std::endl;
~~~
__tabsMid
~~~ python
print("error:", xerus.frob_norm(x - xerus.Tensor(ttx)))
~~~
__tabsEnd

The expected output of this whole program now looks something like
~~~
ttA ranks: { 3, 3, 3, 3, 3, 3, 3, 3 }
residual: 4.93323e-11
error: 1.48729e-20
~~~

We could now also print the solution with `x.to_string()` or store in in a file to plot it with another program.
We could change our operator to define other FEM systems and many more things. As this is only a very short
introduction though, we will stop here and refer the interested reader to either the following more detailed guides or
to their own curiosity.


The full source code of this tutorial looks as follows

__tabsStart
~~~ cpp
{% include examples/qtt.cpp %}
~~~
__tabsMid
~~~ python
{% include examples/qtt.py %}
~~~
__tabsEnd

