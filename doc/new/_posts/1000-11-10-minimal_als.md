---
layout: post
title: "The ALS Algorithm"
date: 1000-12-10
topic: "Examples"
section: "Examples"
---
__tabsInit
# The ALS Algorithm

Implementing the Alternating Least Squares (ALS) algorithm (also known as single-site DMRG) for the first time was the most 
important step for us to understand the TT format and its 
intricacies. We still think that it is a good point to start so we want to provide a simple implementation of the ALS
algorithm as an example. Using the `xerus` library this will be an efficient implementation using only about 100 lines of code (without comments).

## Introduction
The purpose of this page is not to give a full derivation of the ALS algorithm. The interested reader is instead refered to the
original publications on the matter. Let us just shortly recap the general idea of the algorithm to refresh your memory though.

Solving least squares problems of the form $\operatorname{argmin}_x \\|Ax - b\\|^2$ for large dimensions is a difficult endeavour. 
Even if $x$ and $b$ are given in the TT-Tensor and $A$ in the TT-Operator format with small ranks this is far from trivial.
There is a nice property of the TT format though, that we can use to construct a Gauss-Seidel-like iterative scheme: the linear
dependence of the represented tensor on all of its component tensors.
Due to this linearity, we can formulate a smaller subproblem: fixing all but one components of $x$, the resulting minimization
problem is again of the form of a least squares problem as above but with a projected $\hat b = Pb$ and with a smaller $\hat A=PAP^T$. 
In practice, these projections can be obtained by simply contracting all fixed components of $x$ to $A$ and $b$ (assuming all
fixed components are orthogonolized). The ALS algorithm will now simply iterate over the components of $x$ and solve these smaller subproblems.

There are a few things we should note before we start implementing this algorithm
* It is enough to restrict ourselves to the case of symmetric positive-semidefinite operators $A$. Any non-symmetric problem can be solved by setting $A'=A^TA$ and $b' = A^Tb$.
* We should always move our core of $x$ to the position currrently being optimized to make our lives easier (for several reasons...).
* Calculating the local operators $\hat A$ for components $i$ and $i+1$ is highly redundant. All components of $x$ up to the $i-1$'st have to be contracted with $A$ in both cases. Effectively this means, that we will keep stacks of $x^TAx$ contracted up to the current index ("left" of the current index) as well as contracted at all indices above the currrent one ("right" of it) and similarly for $x^T b$.


## Pseudo-Code
Let us start by writing down the algorithm in pseudo-code and then fill out the steps one by one.

__tabsStart
~~~ cpp
// while we are not done
	// for every position p = 0..degree(x) do
		// local operator = left stack(A) * p'th component of A * right stack(A)
		// local rhs = left stack(b) * p'th component of b * right stack(b)
		// p'th component of x = solution of the local least squares problem
		
		// remove top entry of the right stacks
		// add position p to left stacks
	
	// for every position p = degree(x)..0 do
		// same as above OR simply move core and update stacks
~~~
__tabsMid
~~~ py
# while we are not done
	# for every position p = 0..degree(x) do
		# local operator = left stack(A) * p'th component of A * right stack(A)
		# local rhs = left stack(b) * p'th component of b * right stack(b)
		# p'th component of x = solution of the local least squares problem
		
		# remove top entry of the right stacks
		# add position p to left stacks
	
	# for every position p = degree(x)..0 do
		# same as above OR simply move core and update stacks
~~~
__tabsEnd

## Helper Class
We want our main loop to resemble the above pseudo code as closely as possible, so we have to define some helper functions to
update the stacks. To ensure that all functions work on the same data without passing along references all the time, we will
define a small helper class, that holds all relevant variables: the degree `d` of our problem, the left and right stacks for `A`
and `b`, the TT tensors `A`, `x` and `b` themselves and the norm of `b`. As a parameter of the algorithm we will also store the
maximal number of iterations

__tabsStart
~~~ cpp
class InternalSolver {
	const size_t d;
	
	std::vector<Tensor> leftAStack;
	std::vector<Tensor> rightAStack;
	
	std::vector<Tensor> leftBStack;
	std::vector<Tensor> rightBStack;
	
	TTTensor& x;
	const TTOperator& A;
	const TTTensor& b;
	const double solutionsNorm;
public:
	size_t maxIterations;
	
	InternalSolver(const TTOperator& _A, TTTensor& _x, const TTTensor& _b) 
		: d(_x.degree()), x(_x), A(_A), b(_b), solutionsNorm(frob_norm(_b)), maxIterations(1000)
	{ 
		leftAStack.emplace_back(Tensor::ones({1,1,1}));
		rightAStack.emplace_back(Tensor::ones({1,1,1}));
		leftBStack.emplace_back(Tensor::ones({1,1}));
		rightBStack.emplace_back(Tensor::ones({1,1}));
	}
~~~
__tabsMid
~~~ py
class InternalSolver :
	def __init__(self, A, x, b):
		self.A = A
		self.b = b
		self.x = x
		self.d = x.degree()
		self.solutionsNorm = b.frob_norm()
		
		self.leftAStack = [ xe.Tensor.ones([1,1,1]) ]
		self.leftBStack = [ xe.Tensor.ones([1,1]) ]
		self.rightAStack = [ xe.Tensor.ones([1,1,1]) ]
		self.rightBStack = [ xe.Tensor.ones([1,1]) ]
		
		self.maxIterations = 1000
~~~
__tabsEnd

Note here, that we initialized the stacks with tensors of dimensions $1\times 1\times 1$ (respectively $1\times 1$). By doing this
we won't have to distinguish between the first, last or any other component that is being optimized. As the first and last component
of a TT Tensor already have an additional mode of dimension 1 in `xerus`, we can simply contract them to these tensors containing
a 1 entry as if they were any of the middle components.

## Updating the Stacks
To add the next entry to the left (right) stacks we have to contract the next components to the previous results. To increase
readability of the equations, we first store references to the $p$'th components in `Ai`, `xi` and `bi` before we contract them
using `xerus`'s indexed equations.

__tabsStart
~~~ cpp
	void push_left_stack(const size_t _position) {
		Index i1, i2, i3, j1 , j2, j3, k1, k2;
		const Tensor &xi = x.get_component(_position);
		const Tensor &Ai = A.get_component(_position);
		const Tensor &bi = b.get_component(_position);
		
		Tensor tmpA, tmpB;
		tmpA(i1, i2, i3) = leftAStack.back()(j1, j2, j3)
				*xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3);
		leftAStack.emplace_back(std::move(tmpA));
		tmpB(i1, i2) = leftBStack.back()(j1, j2)
				*xi(j1, k1, i1)*bi(j2, k1, i2);
		leftBStack.emplace_back(std::move(tmpB));
	}
	
	
	void push_right_stack(const size_t _position) {
		Index i1, i2, i3, j1 , j2, j3, k1, k2;
		const Tensor &xi = x.get_component(_position);
		const Tensor &Ai = A.get_component(_position);
		const Tensor &bi = b.get_component(_position);
		
		Tensor tmpA, tmpB;
		tmpA(i1, i2, i3) = xi(i1, k1, j1)*Ai(i2, k1, k2, j2)*xi(i3, k2, j3)
				*rightAStack.back()(j1, j2, j3);
		rightAStack.emplace_back(std::move(tmpA));
		tmpB(i1, i2) = xi(i1, k1, j1)*bi(i2, k1, j2)
				*rightBStack.back()(j1, j2);
		rightBStack.emplace_back(std::move(tmpB));
	}
~~~
__tabsMid
~~~ py
	def push_left_stack(self, pos) :
		i1,i2,i3, j1,j2,j3, k1,k2 = xe.indices(8)
		Ai = self.A.get_component(pos)
		xi = self.x.get_component(pos)
		bi = self.b.get_component(pos)
		
		tmpA = xe.Tensor()
		tmpB = xe.Tensor()
		tmpA(i1, i2, i3) << self.leftAStack[-1](j1, j2, j3)\
				*xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3)
		self.leftAStack.append(tmpA)
		tmpB(i1, i2) << self.leftBStack[-1](j1, j2)\
				*xi(j1, k1, i1)*bi(j2, k1, i2)
		self.leftBStack.append(tmpB)
	
	def push_right_stack(self, pos) :
		i1,i2,i3, j1,j2,j3, k1,k2 = xe.indices(8)
		Ai = self.A.get_component(pos)
		xi = self.x.get_component(pos)
		bi = self.b.get_component(pos)
		
		tmpA = xe.Tensor()
		tmpB = xe.Tensor()
		tmpA(j1, j2, j3) << xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3) \
				* self.rightAStack[-1](i1, i2, i3)
		self.rightAStack.append(tmpA)
		tmpB(j1, j2) << xi(j1, k1, i1)*bi(j2, k1, i2) \
				* self.rightBStack[-1](i1, i2)
		self.rightBStack.append(tmpB)
~~~
__tabsEnd

We will also use a small helper function to calculate the residual norm of the current solution:

__tabsStart
~~~ cpp
	double calc_residual_norm() {
		Index i,j;
		return frob_norm(A(i/2, j/2)*x(j&0) - b(i&0)) / solutionsNorm;
	}
~~~
__tabsMid
~~~ py
	def calc_residual_norm(self) :
		i,j = xe.indices(2)
		return xe.frob_norm(self.A(i/2, j/2)*self.x(j&0) - self.b(i&0)) / self.solutionsNorm
~~~
__tabsEnd


## The Main Loop
With these helper functions we can now start to write the main optimization loop. Before we start with the optimization though
we have to ensrure that the core of `x` is at the right position and create the full right stacks - with the above helper functions 
this has become very easy:

__tabsStart
~~~ cpp
	void solve() {
		x.move_core(0, true);
		for (size_t pos = d-1; pos > 0; --pos) {
			push_right_stack(pos);
		}
		// ...
	}
~~~
__tabsMid
~~~ py
	def solve(self) :
		self.x.move_core(0, True)
		for pos in reversed(xrange(1, self.d)) :
			self.push_right_stack(pos)
		
		# ...
~~~
__tabsEnd

The pseudo code now starts with "while we are not done". We have already created a `maxIterations` variable to decide whether
we are "done", but we typically also want an endcriterion based on the current residual. For the purpose of this example we will
simply check, whether the residual decreased by at least 1% within the last 10 sweeps. To do this, we have to keep a record
of the residuals after every sweep:

__tabsStart
~~~ cpp
	void solve() {
		// ...
		std::vector<double> residuals(10, 1000.0);
		
		for (size_t itr = 0; itr < maxIterations; ++itr) {
			residuals.push_back(calc_residual_norm());
			
			if (residuals.back()/residuals[residuals.size()-10] > 0.99) {
				return; 
			}
			
			// ...
		}
	}
~~~
__tabsMid
~~~ py
	def solve(self) :
		# ...
		residuals = [1000]*10
		
		for itr in xrange(self.maxIterations) :
			residuals.append(self.calc_residual_norm())
			
			if residuals[-1]/residuals[-10] > 0.99 :
				return
			
			# ...
~~~
__tabsEnd


The individual iterations now consist of a sweep from left to right and one back to the left. For every position we will construct
the local operator and right-hand-side via indexed equations, solve this local problem and store the result in the corresponding
component of `x`

__tabsStart
~~~ cpp
for (size_t corePosition = 0; corePosition < d; ++corePosition) {
	Tensor op, rhs;
	
	const Tensor &Ai = A.get_component(corePosition);
	const Tensor &bi = b.get_component(corePosition);
	
	op(i1, i2, i3, j1, j2, j3) = leftAStack.back()(i1, k1, j1)
	                             * Ai(k1, i2, j2, k2)
	                             * rightAStack.back()(i3, k2, j3);
	rhs(i1, i2, i3) = leftBStack.back()(i1, k1) 
	                  * bi(k1, i2, k2) 
	                  * rightBStack.back()(i3, k2);
	
	xerus::solve(x.component(corePosition), op, rhs, 0);
	
	// move core and update stacks
	// ...
}


~~~
__tabsMid
~~~ py
for pos in xrange(self.d):
	op = xe.Tensor()
	rhs = xe.Tensor()
	
	Ai = self.A.get_component(pos)
	bi = self.b.get_component(pos)
	
	op(i1, i2, i3, j1, j2, j3) << self.leftAStack[-1](i1, k1, j1) \
	                              * Ai(k1, i2, j2, k2) \
	                              * self.rightAStack[-1](i3, k2, j3)
	rhs(i1, i2, i3) << self.leftBStack[-1](i1, k1) \
	                   * bi(k1, i2, k2) \
	                   * self.rightBStack[-1](i3, k2)
	
	tmp = xe.Tensor()
	tmp(i1&0) << rhs(j1&0) / op(j1/2, i1/2)
	self.x.set_component(pos, tmp)
	
	# move core and update stacks
	# ...
~~~
__tabsEnd

To prepare for the next position we simply have to move the core, pop the top element of the right stacks and add the new element
to the left stacks (unless this already was the last position of course!)

__tabsStart
~~~ cpp
// move core and update stacks
if (corePosition+1 < d) {
	x.move_core(corePosition+1, true);
	push_left_stack(corePosition);
	rightAStack.pop_back();
	rightBStack.pop_back();
}
~~~
__tabsMid
~~~ py
# move core and update stacks
if pos+1 < self.d :
	self.x.move_core(pos+1, True)
	self.push_left_stack(pos)
	self.rightAStack.pop()
	self.rightBStack.pop()
~~~
__tabsEnd


For the sweep from right to left we could now write a similar loop, for the sake of simplicity we will just move the core and 
stacks back to the left though (it seems to depend on the example which of these two variants is faster). So the full iteration looks
as follows:

__tabsStart
~~~ cpp
for (size_t itr = 0; itr < maxIterations; ++itr) {
	residuals.push_back(calc_residual_norm());
	
	if (residuals.back()/residuals[residuals.size()-10] > 0.99) {
		return; // We are done!
	}
	
	// Sweep Left -> Right
	for (size_t corePosition = 0; corePosition < d; ++corePosition) {
		Tensor op, rhs;
		
		const Tensor &Ai = A.get_component(corePosition);
		const Tensor &bi = b.get_component(corePosition);
		
		op(i1, i2, i3, j1, j2, j3) = leftAStack.back()(i1, k1, j1)
		                             * Ai(k1, i2, j2, k2)
		                             * rightAStack.back()(i3, k2, j3);
		rhs(i1, i2, i3) = leftBStack.back()(i1, k1) 
		                  * bi(k1, i2, k2) 
		                  * rightBStack.back()(i3, k2);
		
		xerus::solve(x.component(corePosition), op, rhs);
		
		if (corePosition+1 < d) {
			x.move_core(corePosition+1, true);
			push_left_stack(corePosition);
			rightAStack.pop_back();
			rightBStack.pop_back();
		}
	}
	
	
	// Sweep Right -> Left : only move core and update stacks
	x.move_core(0, true);
	for (size_t corePosition = d-1; corePosition > 0; --corePosition) {
		push_right_stack(corePosition);
		leftAStack.pop_back();
		leftBStack.pop_back();
	}
	
}
~~~
__tabsMid
~~~ py
for itr in xrange(self.maxIterations) :
	residuals.append(self.calc_residual_norm())
	if residuals[-1]/residuals[-10] > 0.99 :
		return
	
	# sweep left -> right
	for pos in xrange(self.d):
		op = xe.Tensor()
		rhs = xe.Tensor()
		
		Ai = self.A.get_component(pos)
		bi = self.b.get_component(pos)
		
		op(i1, i2, i3, j1, j2, j3) << self.leftAStack[-1](i1, k1, j1) \
		                              * Ai(k1, i2, j2, k2) \
		                              * self.rightAStack[-1](i3, k2, j3)
		rhs(i1, i2, i3) << self.leftBStack[-1](i1, k1) \
		                   * bi(k1, i2, k2) \
		                   * self.rightBStack[-1](i3, k2)
		
		tmp = xe.Tensor()
		tmp(i1&0) << rhs(j1&0) / op(j1/2, i1/2)
		self.x.set_component(pos, tmp)
		
		if pos+1 < self.d :
			self.x.move_core(pos+1, True)
			self.push_left_stack(pos)
			self.rightAStack.pop()
			self.rightBStack.pop()
	
	
	# right -> left, only move core and update stack
	self.x.move_core(0, True)
	for pos in reversed(xrange(1,self.d)) :
		self.push_right_stack(pos)
		self.leftAStack.pop()
		self.leftBStack.pop()
~~~
__tabsEnd


## Finishing Touches
All we need now is a smal function that creates an object of our helper class and calls the `.solve()` method

__tabsStart
~~~ cpp
void simpleALS(const TTOperator& _A, TTTensor& _x, const TTTensor& _b)  {
	InternalSolver solver(_A, _x, _b);
	return solver.solve();
}
~~~
__tabsMid
~~~ py
def simpleALS(A, x, b) :
	solver = InternalSolver(A, x, b)
	solver.solve()
~~~
__tabsEnd

To see what is happening in the algorithm we can also add some output to the iterations along the lines of

__tabsStart
~~~ cpp
XERUS_LOG(simpleALS, "Iteration: " << itr << " Residual: " << residuals.back());
~~~
__tabsMid
~~~ py
print("Iteration:",itr, "Residual:", residuals[-1])
~~~
__tabsEnd


## Testing our Implementation
To test our freshly implemented ALS, we can use a random TT-Operator (after symmetrizing it via $A'=A^T A$). By constructing a
random low-rank solution $s$ and providing the right-hand-side $b=A' s$ to the algorithm, we can ensure that a low-rank solution
actually exists. As a starting point we will furthermore provide a random initial TT tensor `x` of the correct rank.

As we know the correct solution via this construction and the random operator `A` is invertible with high probability, we can then
calculate the residual and the actual error to verify that our algorithm works as intended. (We have no proof, that the algorithm
always finds the correct solution in this setting, but empirically this is the case...)

__tabsStart
~~~ cpp
int main() {
	Index i,j,k;
	
	auto A = TTOperator::random(std::vector<size_t>(16, 4), std::vector<size_t>(7,2));
	A(i/2,j/2) = A(i/2, k/2) * A(j/2, k/2);
	
	auto solution = TTTensor::random(std::vector<size_t>(8, 4), std::vector<size_t>(7,3));
	TTTensor b;
	b(i&0) = A(i/2, j/2) * solution(j&0);
	
	auto x = TTTensor::random(std::vector<size_t>(8, 4), std::vector<size_t>(7,3));
	simpleALS(A, x, b);
	
	XERUS_LOG(info, "Residual: " << frob_norm(A(i/2, j/2) * x(j&0) - b(i&0))/frob_norm(b));
	XERUS_LOG(info, "Error: " << frob_norm(solution-x)/frob_norm(x));
}
~~~
__tabsMid
~~~ py
if __name__ == "__main__":
	i,j,k = xe.indices(3)
	
	A = xe.TTOperator.random([4]*16, [2]*7)
	A(i/2,j/2) << A(i/2, k/2) * A(j/2, k/2)
	
	solution = xe.TTTensor.random([4]*8, [3]*7)
	b = xe.TTTensor()
	b(i&0) << A(i/2, j/2) * solution(j&0)
	
	x = xe.TTTensor.random([4]*8, [3]*7)
	simpleALS(A, x, b)
	
	print("Residual:", xe.frob_norm(A(i/2, j/2) * x(j&0) - b(i&0))/xe.frob_norm(b))
	print("Error:", xe.frob_norm(solution-x)/xe.frob_norm(x))
~~~
__tabsEnd


## Complete Sourcecode
The full source code of this example looks as follows

__tabsStart
~~~ cpp
{% include examples/als.cpp %}
~~~
__tabsMid
~~~ python
{% include examples/als.py %}
~~~
__tabsEnd


