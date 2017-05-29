---
layout: post
title: "Decompositions and Solve"
date: 2000-11-26
topic: "Basic Usage"
section: "Documentation"
---
__tabsInit
# Decompositions and Solve

In matrix calculus the decomposition of a matrix into the matrix product of several matrices with special properties (eg.
into an orthogonal and a triangular matrix (QR) or orthogonal matrices and diagonal of singular values (SVD)) are among the
most powerful tools to devise numerical algorithms. In the case of tensors of higher degree it is necessary to indicate along
which modes the decomposition is supposed to happen, so `xerus` uses the notation of indexed equations explained in the previous
chapter ([Indices and Equations](/indices)).

## QR Decompositions
To provide an intuitive approach to decompositions, `xerus` uses the assignment of multiple tensors with a single operator to 
denote them. Here `(Q, R) = QR(A)` reads as "Q and R are defined as the QR-Decomposition of A", even though we naturally have
to provide indices to make this line well defined:

__tabsStart
~~~ cpp
// perform QR decomposition of A and store result in Q and R
(Q(i,r), R(r,j)) = xerus.QR(A(i,j));
~~~
__tabsMid
~~~ py
# perform QR decomposition of A and store result in Q and R
(Q(i,r), R(r,j)) << xerus.QR(A(i,j))
~~~
__tabsEnd

In these decompositions we can distribute the modes of the original tensor as we please, but to be well defined, `Q` and `R` must
only share one and exactly one index.

__tabsStart
~~~ cpp
// well defined QR decomposition of a degree 4 tensor
(Q(i,r,k), R(l,r,j)) = xerus.QR(A(i,j,k,l));
// invalid: (Q(i,r,s), R(r,s,j)) = xerus.QR(A(i,j));
~~~
__tabsMid
~~~ py
# well defined QR decomposition of a degree 4 tensor
(Q(i,r,k), R(l,r,j)) << xerus.QR(A(i,j,k,l))
# invalid: (Q(i,r,s), R(r,s,j)) << xerus.QR(A(i,j))
~~~
__tabsEnd

For convenience `xerus` defines four variants of the QR decomposition. Assuming the input is of size $m\times n$, $min = \operatorname{min}(m,n)$
and $r$ the rank of the input we have following resulting objects:

<table class="table">
  <thead>
    <tr>
      <th style="border: none;"></th>
      <th style="border-bottom: none; border-left: 2px solid #ddd;">Left-Hand-Side</th>
      <th style="border: none;"></th>
      <th style="border-bottom: none; border-left: 2px solid #ddd; ">Right-Hand-Side</th>
      <th style="border: none;"></th>
    </tr>
    <tr>
      <th style="border-top: none;">Decomposition</th>
      <th style="border-top: none; border-left: 2px solid #ddd;">Property</th>
      <th style="border-top: none;">Dimensions</th>
      <th style="border-top: none; border-left: 2px solid #ddd;">Property</th>
      <th style="border-top: none;">Dimensions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">xerus.QR</code></td>
      <td style="border-left: 2px solid #ddd;">orthogonal</td>
      <td>$m\times min$</td>
      <td style="border-left: 2px solid #ddd;">upper triangular</td>
      <td>$min\times n$</td>
    </tr>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">xerus.RQ</code></td>
      <td style="border-left: 2px solid #ddd;">upper triangular</td>
      <td>$m\times min$</td>
      <td style="border-left: 2px solid #ddd;">orthogonal</td>
      <td>$min\times n$</td>
    </tr>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">xerus.QC</code></td>
      <td style="border-left: 2px solid #ddd;">orthogonal</td>
      <td>$m\times r$</td>
      <td style="border-left: 2px solid #ddd;">upper or lower triangular</td>
      <td>$r\times n$</td>
    </tr>
    <tr>
      <td style="text-align: right"><code class="highlighter-rouge">xerus.CQ</code></td>
      <td style="border-left: 2px solid #ddd;">upper or lower triangular</td>
      <td>$m\times r$</td>
      <td style="border-left: 2px solid #ddd;">orthogonal</td>
      <td>$r\times n$</td>
    </tr>
  </tbody>
</table>


## Singular Value Decompositions

The Singular Value Decomposition in `xerus` is called very much like the `QR` decomposition:

__tabsStart
~~~ cpp
// calculate the SVD of A and store the resulting matrices in U, S and Vt
(U(i,r1), S(r1,r2), Vt(r2,j)) = xerus.SVD(A(i,j));
~~~
__tabsMid
~~~ py
# calculate the SVD of A and store the resulting matrices in U, S and Vt
(U(i,r1), S(r1,r2), Vt(r2,j)) << xerus.SVD(A(i,j))
~~~
__tabsEnd

In this form it is rank-revealing (so `S` is of dimensions $r\times r$ instead of $\operatorname{min}(m,n)\times\operatorname{min}(m,n)$)
and exact, but it is possible to pass optional arguments to use it as a truncated SVD.

__tabsStart
~~~ cpp
// calculate the SVD, truncated to at most 5 singular values
size_t numSingularVectors = 5;
// or until a singular value is smaller than 0.01 times the maximal singular value
double epsilon = 0.01;
(U(i,r1), S(r1,r2), Vt(r2,j)) = xerus.SVD(A(i,j), numSingularVectors, epsilon);
~~~
__tabsMid
~~~ py
# calculate the SVD, truncated to at most 5 singular values
numSingularVectors = 5
# or until a singular value is smaller than 0.01 times the maximal singular value
epsilon = 0.01
(U(i,r1), S(r1,r2), Vt(r2,j)) = xerus.SVD(A(i,j), maxRank=numSingularVectors, eps=epsilon);
~~~
__tabsEnd


## Solving Linear Equations

A common application for matrix decompositions is to solve matrix equations of the form $A\cdot x = b$ for $x$ via QR, LU,
LDL$^T$ or Cholesky decompositions. In xerus this is again provided in indexed notation via the `operator/`.

__tabsStart
~~~ cpp
// solve A(i,j)*x(j) = b(i) for x
x(j) = b(i) / A(i,j);
~~~
__tabsMid
~~~ py
# solve A(i,j)*x(j) = b(i) for x
x(j) << b(i) / A(i,j)
~~~
__tabsEnd

Depending on the representation and some properties of `A`, `xerus` will automatically choose one of the above mentioned decompositions
to solve the system.
