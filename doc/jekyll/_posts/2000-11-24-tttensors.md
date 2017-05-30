---
layout: post
title: "TT-Tensors"
date: 2000-11-24
topic: "Basic Usage"
section: "Documentation"
---
__tabsInit
# Tensor Train / MPS Tensors and Operators

__warnStart

This Chapter is still **work in progress**.
Currently missing are sections on:
* introduction with short explanation of the TT format (both TTTensor and TTOperator)
* advanced usage: 
  * lazy evaluation
  * access to nodes
  * require_correct_format()
  * sparse components

__warnEnd

__infoStart

In the following, whenever something is said for the `TTTensor` class, it also holds true for the `TTOperator` class unless
otherwise stated. Everything said explicitely about the `TTOperator` class though does not hold for the `TTTensor` class.

__infoEnd


## Construction

The basic constructors for `TTTensor` and `TTOperator` objects work analogous to the respective `Tensor` constructors.

__tabsStart
~~~ cpp
// create a degree 0 TTTensor with entry 0
xerus::TTTensor();

// create a 1x1x1x1 TTOperator with entry 0
xerus::TTOperator(4);

// create a 2x3x4 TTTensor with all entries equal to 0
xerus::TTTensor({2,3,4});

// create a random TTTensor of dimensions 4x4x4 and ranks (2, 2)
xerus::TTTensor::random({4,4,4}, {2,2});
// and a random TTOperator of dimensions 5x5x5x5 and ranks (2)
xerus::TTOperator::random({5,5,5,5}, {2});

// construct an identity operator of dimensions (2x3) x (2x3)
xerus::TTOperator::identity({2,3,2,3});

// a 2x3x4 TTTensor with all entries = 1
xerus::TTTensor::ones({2,3,4});

// a 3x4x3x4 TTTensor with superdiagonal = 1 (where all 4 indices coincide) and = 0 otherwise
// note that this TTTensor is of rank r = min(dimension)
xerus::TTTensor::kronecker({3,4,3,4});

// a 2x2x2 TTTensor with a 1 in position {1,1,1} and 0 everywhere else
xerus::TTTensor::dirac({2,2,2}, {1,1,1});
~~~
__tabsMid
~~~ py
# create a degree 0 TTTensor with entry 0
xerus.TTTensor()

# create a 1x1x1x1 TTOperator with entry 0
xerus.TTOperator(4)

# create a 2x3x4 TTTensor with all entries equal to 0
xerus.TTTensor([2,3,4])

# create a random TTTensor of dimensions 4x4x4 and ranks (2, 2)
xerus.TTTensor.random([4,4,4], [2,2])
# and a random TTOperator of dimensions 5x5x5x5 and ranks (2)
xerus.TTOperator.random([5,5,5,5], [2])

# construct an identity operator of dimensions (2x3) x (2x3)
xerus.TTOperator.identity([2,3,2,3])

# a 2x3x4 TTTensor with all entries = 1
xerus.TTTensor.ones([2,3,4])

# a 3x4x3x4 TTTensor with superdiagonal = 1 (where all 4 indices coincide) and = 0 otherwise
# note that this TTTensor is of rank r = min(dimension)
xerus.TTTensor.kronecker([3,4,3,4])

# a 2x2x2 TTTensor with a 1 in position {1,1,1} and 0 everywhere else
xerus.TTTensor.dirac([2,2,2], [1,1,1]);
~~~
__tabsEnd

Note, that all of the above constructors, except for the identity operator, work both for `TTTensor`s and `TTOperator`s.

In case the TTTensor was derived mathematically, every component of it is likely available individually. In such cases the TTTensor
can also be constructed per component in xerus.

__tabsStart
~~~ cpp
// construct a TTTensor of the correct degree
xerus::TTTensor U(3);

// and set every component individually to previously calculated Tensors
U.set_component(0, c0);
U.set_component(1, c1);
U.set_component(2, c2);

// ensure that U is canonicalized
U.canonicalize_left();
~~~
__tabsMid
~~~ py
# construct a TTTensor of the correct degree
U = xerus.TTTensor(3)

# and set every component individually to previously calculated Tensors
U.set_component(0, c0)
U.set_component(1, c1)
U.set_component(2, c2)

# ensure that U is canonicalized
U.canonicalize_left()
~~~
__tabsEnd

Note two things here:
1. The components passed to [.set_component()](__doxyref(xerus::TTNetwork::set_component)) have to be of dimensions $r_{i-1}\times n_i\times r_i$ (for TTTensors) respectively $r_{i-1}\times m_i\times n_i\times r_i$ (for TTOperators). Neither ranks nor external dimensions have to coincide with the respective values the TTTensor had before the call to `set_component`, but before any other function is used, the ranks of neighboring components have to be consistent.
2. As `xerus` will not check for orthogonality of the components, it will assume that the thusly constructed TTTensor is not canonicalized (unless only the core was changed). To ensure numerical stability of the calcuations you should always finish your construction with a call to [.canonicalize_left()](__doxyref(xerus::TTNetwork::canonicalize_left)) or fix the canonicalization in another way (see below).


## Core Position and Canonicalization
Per default `xerus` will always construct TTTensors in which all but the first component are orthogonal (ie. canonical with core position 0 / left-canonical).
Unless you explicitely construct TTTensor via `.set_component()` or change their core position you thus don't have to worry about
this detail.

__warnStart

If two or more TTTensors are added, contracted (with TTOperators) or multiplied entrywise, that do not have the same core position (or if one of them is not canonicalized),
`xerus` will **not** attempt to canonicalize the result. As such calculations can become very unstable numerically you should
be very cautious when you do this.

__warnEnd

If you want to ensure that several TTTensors are compatible, you can check [.canonicalized](__doxyref(xerus::TTNetwork::canonicalized))
to see if the object is canonicalized and [.corePosition](__doxyref(xerus::TTNetwork::corePosition)) to see which component is not
orthogonal.

Canonicalization can be achieved by a call to [.move_core()](__doxyref(xerus::TTNetwork::move_core)), 
[.canonicalize_left()](__doxyref(xerus::TTNetwork::canonicalize_left)) (equivalent to `.move_core(0)`) or 
[.canonicalize_right()](__doxyref(xerus::TTNetwork::canonicalize_right)) (equivalent to `.move_core(.degree()-1)`. (There is also 
[.assume_core_position()](__doxyref(xerus::TTNetwork::assume_core_position)) if you really know what you are doing.)


## Casting from and to Tensors
Fully contracting a TTTensor network to a single tensor is simple, both theoretically and in `xerus`: simply cast the TTTensor to
a Tensor. The other direction is a bit more complicated. Per default, `xerus` will calculate the exact (up to machine precision)
rank revealing decomposition of the input tensor via successive SVDs. Optional arguments can modify this behaviour though to
instead calculate a quasi-best approximation with a given rank or precision.

__tabsStart
~~~ cpp
// construct a random tensor
auto A = xerus::Tensor::random({10,10,10,10});

// calculate a quasi-best approximation of A with ranks (8,8,8)
xerus::TTTensor ttA(A, {8,8,8});

// contract the TT network to a single tensor
xerus::Tensor A2(ttA);

// and output the error of the low-rank approximation
XERUS_LOG(info, "approximation quality: " << frob_norm(A2-A)/frob_norm(A));
~~~
__tabsMid
~~~ py
# construct a random tensor
A = xerus.Tensor.random([10,10,10,10])

# calculate a quasi-best approximation of A with ranks (8,8,8)
ttA = xerus.TTTensor(A, [8,8,8])

# contract the TT network to a single tensor
A2 = xerus.Tensor(ttA)

# and output the error of the low-rank approximation
print("approximation quality:", frob_norm(A2-A)/frob_norm(A))
~~~
__tabsEnd


## Ranks and Rounding
To decrease the rank of an already existing TTTensor, `xerus` offers the [.round()](__doxyref(xerus::TTNetwork::round)) method.
It can be called with either a single rank that will be used as a maximum for all positions, a maximal rank tuple or a precision
`eps` that determines how small singular values need to be relative to the largest singular valueto be discarded.

Determining the current rank can be done with [.ranks()](__doxyref(xerus::TTNetwork::ranks)) or with [.rank(i)](__doxyref(xerus::TTNetwork::rank)).

__tabsStart
~~~ cpp
// construct a random tensor
auto A = xerus::Tensor::random({5,5,5,5});

// calculate an exact TT decomposition of A
xerus::TTTensor ttA(A);

// we expect this to be of maximal rank (5, 25, 5)
XERUS_LOG(info, "ranks: " << ttA.ranks());

// from maximal rank (25) down to 1
for (size_t r = 25; r>0; --r) {
	// round ttA down to rank r
	ttA.round(r);
	
	// and determine the approximation error
	XERUS_LOG(info, "rank: " << ttA.rank(1) << " approximation error: " 
		<< frob_norm(xerus::Tensor(ttA)-A)/frob_norm(A));
}
~~~
__tabsMid
~~~ py
# construct a random tensor
A = xerus.Tensor.random([5,5,5,5])

# calculate an exact TT decomposition of A
ttA = xerus.TTTensor(A)

# we expect this to be of maximal rank (5, 25, 5)
print("ranks:", ttA.ranks())

# from maximal rank (25) down to 1
for r in reverse(xrange(1, 25)) :
	# round ttA down to rank r
	ttA.round(r)
	
	# and determine the approximation error
	print("rank:", ttA.rank(1), "approximation error:", 
		xerus.frob_norm(xerus.Tensor(ttA)-A)/A.frob_norm())
~~~
__tabsEnd

Rounding to a fixed rank tuple, ie. finding a quasi-best approximation on the set of TTTensors of this smaller rank, is a projection
that is often used in Riemannian algorithms and is necessary because many operations on TTTensors increase their rank. 
(For a more precise but also more time consuming alternative to `.round()` see the chapter on the [ALS algorithm](/als).)

<table class="table">
  <thead>
    <tr>
      <th style="text-align: right">Operation</th>
      <th style="border-left: 2px solid #ddd;">Rank of the Result</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">TTTensor + TTTensor</code></td>
      <td style="border-left: 2px solid #ddd;">$r_1 + r_2$</td>
    </tr>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">TTOperator * TTTensor</code></td>
      <td style="border-left: 2px solid #ddd;">$r_1 * r_2$</td>
    </tr>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">entrywise_product(TTTensor, TTTensor)</code></td>
      <td style="border-left: 2px solid #ddd;">$r_1 * r_2$</td>
    </tr>
  </tbody>
</table>

This change in ranks can also be used on purpose to increase the rank of a candidate solution.

__tabsStart
~~~ cpp
// increase the ranks of ttA by 1 (with high probability)
ttA += ttA.frob_norm() * 1e-10 
	* xerus::TTTensor::random(ttA.dimensions, std::vector<size_t>(ttA.degree()-1,1));
~~~
__tabsMid
~~~ py
# increase the ranks of ttA by 1 (with high probability)
ttA += ttA.frob_norm() * 1e-10 \
	* xerus.TTTensor.random(ttA.dimensions, [1]*(ttA.degree()-1));
~~~
__tabsEnd

Though above method will naturally only increase the ranks with high probability and up to the maximal possible ranks.

Finally we want to mention the soft thresholding operator [.soft_threshold()](__doxyref(xerus::TTNetwork::soft_threshold)) as it
can also decrease the rank of a TTTensor. Note that it does not provide a quasi-best approximation though.


## Output and Storing
In analogy to the `Tensor` class, TTTensors can be queried for their degree ([.degree()](__doxyref(xerus::TensorNetwork::degree))),
their dimensions ([.dimensions](__doxyref(xerus::TensorNetwork::dimensions))) and Frobenius norm ([.frob_norm()](__doxyref(xerus::TTNetwork::frob_norm))),
though not for their 1-norm (as it is hard to calculate).

Of particular interest are also the above mentioned ranks of a given TTTensor that can be obtained with [.ranks()](__doxyref(xerus::TTNetwork::ranks))
or with [.rank(i)](__doxyref(xerus::TTNetwork::rank)).

TTTensors themselves have no `.to_string` method as it is not quite obvious what that function should do. To obtain an output
similar to that of the Tensor class you can cast the object to Tensor though.

__tabsStart
~~~ cpp
// construct a random TTTensor of rank 2
auto U = xerus::TTTensor::random({4,4,4}, 2);

// allow piping of vectors to streams
using xerus::misc::operator<<;

// query for different information:
XERUS_LOG(info, "degree: " << U.degree());
XERUS_LOG(info, "dimensions: " << U.dimensions);
XERUS_LOG(info, "frob_norm: " << U.frob_norm());
XERUS_LOG(info, "ranks: " << U.ranks());
XERUS_LOG(info, "entries:\n" << xerus::Tensor(U).to_string());
~~~
__tabsMid
~~~ py
# construct a random TTTensor of rank 2
U = xerus.TTTensor.random({4,4,4}, 2);

# query for different information:
print("degree:", U.degree())
print("dimensions:", U.dimensions)
print("frob_norm:", U.frob_norm())
print("ranks:", U.ranks())
print("entries:\n", xerus.Tensor(U).to_string())
~~~
__tabsEnd

Storing and restoring TTTensorsto / from files uses the same functions as Tensors do: [save_to_file()](__doxyref(xerus::misc::save_to_file))
and [load_from_file()](__doxyref(xerus::misc::load_from_file)). These function perform no transformation to other formats.

__tabsStart
~~~ cpp
// construct a random 3x3x3 TTTensor
auto A = xerus::TTTensor::random({3, 3, 3}, 2);

// store the TTTensor to the file "tensor.dat"
xerus::misc::save_to_file(A, "tensor.dat");

// load the TTTensor from the file
auto B = xerus::misc::load_from_file<xerus::TTTensor>("tensor.dat");
// the following line would throw an exception
// auto C = xerus::misc::load_from_file<xerus::Tensor>("tensor.dat");

// check for correct reconstruction
XERUS_LOG(info, "restoration error: " << frob_norm(A-B));
~~~
__tabsMid
~~~ py
# construct a random 3x3x3 TTTensor
A = xerus.TTTensor.random([3, 3, 3], 2)

# store the TTTensor to the file "tensor.dat"
xerus.misc.save_to_file(A, "tensor.dat")

# load the TTTensor from the file
B = xerus.misc.load_from_file("tensor.dat");
# the type of B is determined automatically by the content of the file
print("loaded type:", type(B))

# check for correct reconstruction
print("restoration error:", xerus.frob_norm(A-B))
~~~
__tabsEnd


While the [operator[]](__doxyref(xerus::TensorNetwork::operator[])) cannot be used to set a single entry, it does provide a way
to efficiently calculate it. If you need to read many entries though, there might be more efficient means to do so. To calculate
(close to) all entries you should definitely cast the network to a single Tensor. If you need many entries, but it is inefficient
to calculate the contracted tensor due to its size, using a `SinglepointMeasurementSet` might be the most efficient way (see
the corresponding chapter [Tensor Completion](/completion) to see how to do that).



## Operators and Modifications
All well defined operators (`+`, `-`, `scalar *`, `/ scalar`, `+=`, `-=`, `/= scalar`, `*= scalar`) are implemented for TTTensors with their
canonical meaning.

In analogy to the Tensor class, the methods [.fix_mode()](__doxyref(xerus::TensorNetwork::fix_mode)), 
[.remove_slate()](__doxyref(xerus::TensorNetwork::remove_slate)) and [.resize_mode()](__doxyref(xerus::TensorNetwork::resize_mode)) can
be implemented efficiently and are thus provided by `xerus`. (For an example of these see the [Tensor](/tensor) chapter.)

Similarly the Hadamard product is implemented as [entrywise_product()](__doxyref(xerus::misc::entrywise_product)) though it necessarily
increases the rank as stated above.

Because the indexed equations for TTTensors are somewhat limited (see below), there are two more function provided explicitely for TTTensors:
[.transpose()](__doxygen(xerus::TTTensor::transpose)) to exchange the first half of the modes of a TTOperator with its second half; and
the dyadic product [dyadic_product()](__doxyref(xerus::misc::dyadic_product)) which accepts (either two or) a list of TTTensors to
join in a dyadic product, which is useful e.g. for operators of the common form 
$M_1 \otimes I \otimes I + M_2 \otimes A \otimes I + M_3 \otimes I \otimes A$ as shown in the following example:

__tabsStart
~~~ cpp
// "calculate" M_1, M_2, M3 and A
auto M1 = xerus::TTOperator::random({100,100}, {});
auto M2 = xerus::TTOperator::random({100,100}, {});
auto M3 = xerus::TTOperator::random({100,100}, {});
auto  A = xerus::TTOperator::random({10,10}, {});
auto  I = xerus::TTOperator::identity({10,10});

// construct total operator
auto Op = xerus::dyadic_product({M1, I, I})
	+ xerus::dyadic_product({M2, A, I})
	+ xerus::dyadic_product({M3, I, A});
~~~
__tabsMid
~~~ cpp
# "calculate" M_1, M_2, M3 and A
M1 = xerus.TTOperator.random([100,100], [])
M2 = xerus.TTOperator.random([100,100], [])
M3 = xerus.TTOperator.random([100,100], [])
A  = xerus.TTOperator.random([10,10], [])
I  = xerus.TTOperator.identity([10,10])

# construct total operator
Op = xerus.dyadic_product([M1, I, I]) \
	+ xerus.dyadic_product([M2, A, I]) \
	+ xerus.dyadic_product([M3, I, A])
~~~
__tabsEnd

This is only well defined for TTOperators (and TTTensors), because the index order for the operator functionality is fixed in the format. Because
the Tensor class is less unique, it is not possible to provide the same functionality for the Tensor class (it can be achieved there
via indexed equations though).



## Indexed Equations
Reordering the modes of a Tensor in the TT format is highly non-trivial. As a consequence indexed equations with tensors in a
TT format are somewhat limited - if you want to stay in the TT format. 

__warnStart

At the time of this writing (`xerus` version 3), `xerus` will only recognize that the result of an indexed equation 
is of a TT format if it was achieved by summation, contraction with a TTOperator or multiplications
with scalars - all without a change in index orders. In any other case it will result in a object of the more general TensorNetwork class
(see [Tensor Networks](/tensornetworks)).

__warnEnd

__tabsStart
~~~ cpp
auto U = xerus::TTTensor::random({10,10,10}, {2,2});

// fast index reshuffle
xerus::TensorNetwork A;
A(i,j,k) = U(k,j,i);
// A is not stored as TTTensor even though it would still be
// one from a mathematical point of view.

// slow index reshuffle
xerus::TTTensor B;
B(i,j,k) = U(k,j,i);
// xerus does not recognize the similarity in the networks and will thus
// fully contract U to a single tensor and then use SVDs to construct B
~~~
__tabsMid
~~~ cpp
U = xerus.TTTensor.random([10,10,10], [2,2])

# fast index reshuffle
A = xerus.TensorNetwork()
A(i,j,k) << U(k,j,i)
# A is not stored as TTTensor even though it would still be
# one from a mathematical point of view.

# slow index reshuffle
B = xerus.TTTensor()
B(i,j,k) << U(k,j,i)
# xerus does not recognize the similarity in the networks and will thus
# fully contract U to a single tensor and then use SVDs to construct B
~~~
__tabsEnd

Which method you should choose depends on your usecase. If you only want to contract other tensors to it, it is absolutely fine
to use `A`. As a general Tensor Network it does not define a summation though and calculating the Frobenius norm will be much 
slower than for `U` (for more details see [Tensor Networks](/tensornetworks)).
Calculating `B` was relatively slow, but as a one-time cost it might be acceptable in many usecases.
As the ideal solution performance-wise it would furthermore have been possible to construct `B` via 
[.set_component()](__doxyref(xerus::TTNetwork::set_component)) from the components of `U` ([.get_component()](__doxyref(xerus::TTNetwork::get_component)) )
though that in turn would have been more work to write and might be less understandable to read.

Typical, well-defined indexed equations that TTTensors and TTOperators were built for look as follows

__tabsStart
~~~ cpp
// well defined and efficient indexed equation with TTTensors and TTOperators:
// x = x + \alpha * A^T(Ax - b)
ttx(i&0) = ttx(i&0) + alpha * ttA(j/2, i/2) * ( ttA(j/2, k/2) * ttx(k&0) - ttb(j&0) );
~~~
__tabsMid
~~~ py
# well defined and efficient indexed equation with TTTensors and TTOperators:
# x = x + \alpha * A^T(Ax - b)
ttx(i&0) << ttx(i&0) + alpha * ttA(j/2, i/2) * ( ttA(j/2, k/2) * ttx(k&0) - ttb(j&0) )
~~~
__tabsEnd

