---
layout: post
title: "Nomenclature"
date: 2000-12-15 16:25:06 +0000
topic: "Introduction"
section: "Documentation"
---


# Nomenclature

Tensor (Network) methods were developed independently in several different fields. As such there is a large variety of different
names for the same concepts and at times even several concepts for the same name, depending on the context. To avoid confusion 
we want to explain most terms as they are used throughout this library.

It is not strictly speaking necessary to read this chapter to successfully use the library, but it might very well help you find 
the functions you were looking for if you are used to slightly different notation.

## Tensors

For us, a **tensor** is always a multidimensional array of real numbers. With the regular tensor product $ \otimes $ this means 
that for any tensor $ T $:

$$ T \in \mathbb{R}^{n_1 \times n_2 \times \cdots \times n_d} = \mathbb{R}^{n_1} \otimes \mathbb{R}^{n_2} \otimes \cdots \otimes \mathbb{R}^{n_d} $$

Here the number $ d $ of subspaces that were joined with the tensor product is called the **degree** (or sometimes **order**) of $ T $.
Just as $ n_i $ was the dimension of the respective subspace $ \mathbb{R}^{n_i} $, it is also the **dimension of the $i$-th mode**
of $T$. The full **dimensions** of $T$ (note the plural), or equivalently the dimensions of the **tensor space** to which $T$ belongs, are given by the ordered $d$-tuple $(n_1, n_2, \dots, n_d)$.

With discrete sets $[n] = \{1,2,\dots,n\}$ we can alternatively define a tensor entrywise as:

$$ T[i_1, i_2, \dots, i_d] \in \mathbb{R}\quad\text{for}\quad i_1\in [n_1], i_2\in [n_2], \dots, i_d\in [n_d]  $$

We call the notation $T[i_1, i_2, \dots, i_d]$ an **indexed tensor** with $d$ **indices** $i_1,\dots,i_d$. The index $i_1$ indexes the first
**mode** of $T$, $i_2$ respectively indexes the second mode and so on. The **dimension of the $j$-th index** $i_j$ is equal to $n_j$.

Instead of indexing a tensor with individual indices, one or more **multiindices** can be used. Every multiindex can be represented
by an ordered tuple. E.g. we could write the last definition of tensors as

$$ T[\mathbb{i}] = T[\mathbb{i}_1, \mathbb{i}_2, \dots] \in \mathbb{R}\quad\text{for}\quad \mathbb{i} \in [n_1]\otimes [n_2] \otimes \cdots \otimes [n_d] $$

Here the **span** (i.e. the degree of the indexspace, here $[n_1]\otimes [n_2] \otimes \cdots \otimes [n_d]$) of the index $\mathbb{i}$
is equal to the degree of $T$. The tensor is thus fully indexed by $\mathbb{i}$. It could alternatively be indexed by two multiindices
of **span** $d/2$ or one with span $2$ and another one with span $d-2$ or ...

An individual value stored in a tensor, e.g.

$$ U[3, 5, 7] = 3.141 $$

is called an **entry** of the tensor. In the example the entry has the **position** $(3, 5, 7)$. If a tensor is equal to $0$
in most positions, it can be stored efficiently in a sparse **representation**. 

**Fixing a mode** to a single value
we would receive a row (e.g. $b[i] = A[2,i]$) or a column (e.g. $b[i] = A[i,5]$) in the matrix case. In the general tensor case 
(e.g. $S[i,j,k] = T[i,3,j,k]$) we call the resulting tensor of degree $d-1$ a **slate** of $T$.

For every permutation $p:[d]\rightarrow[d]$ there is a **reshuffling** (or **reordering**) $R$ such that

$$ 
\begin{align}
R: \mathbb{R}^{n_1 \times n_2 \times \cdots \times n_d} &\rightarrow \mathbb{R}^{n_{p(1)} \times n_{p(2)} \times \cdots \times n_{p(d)}} \\
T&\mapsto S\;:\;S[i_1, i_2, \dots, i_d] = T[i_{p(1)}, i_{p(2)}, \dots, i_{p(d)}]
\end{align}
$$

e.g. for matrices the transposition is a reordering or for tensors of degree $3$ the following is one of five possible (nontrivial) reorderings

$$
\begin{align}
\mathbb{R}^{n_1 \times n_2 \times n_3} &\rightarrow \mathbb{R}^{n_1 \times n_3 \times n2} \\
T&\mapsto S\;:\;S[i,j,k] = T[i,k,j]
\end{align}
$$

The **contraction** of a mode of a tensor with a mode of another tenser is equal to the sum over all tensor products of the slates
of those modes. E.g. any matrix-matrix product or defining a tensor $S$ entrywise as

$$ S[i,j] = \sum_k \sum_l T[k,i,l] \cdot U[j,k,l] $$

In the **Einstein notation** it is customary to perform such sums over all indices that appear exactly twice in a product implicitely,
i.e.

$$ S[i,j] = T[k,i,l] \cdot U[j,k,l] = \sum_k \sum_l T[k,i,l] \cdot U[j,k,l] $$


## Tensor Networks

A **tensor network** is a set of tensors together with a set of contractions between them. It is called a network because
it can be represented by a graph where every **node** represents a tensor and the **links** between them indicate contractions.
Any mode of any tensor in the network that is not part of any contractions is an **external mode** (or **external link**) of the
network.

A very common set of tensor networks are the **Tensor Train Tensors** (or short **TT-Tensors**, also known as Matrix Product States / MPS).
They consist of a linear row of tensors with one external mode each and internal links only between neighbors. In a slight modification
with two external modes per tensor they are called **Tensor Train Operators** (or short **TT-Operator**, also known as Matrix Product
 Operator / MPO) as they are often used in conjunction with TT-Tensors in an operator functionality.

As there is an obvious ordering for the nodes of a TT-Tensor, we also say that the tensor with the first external mode is the 
**first component**, and so on. It its **cannonical form**, all but one component are orthogonalized. The remaining non-orthogonalized
component is called the **core** of the TT-Tensor. If the **core position** is $0$, i.e. the $0$-th component is the core, the
TT-Tensor is in its **left-cannonical** form, respectively **right-cannonical** with core-position $d-1$.

The ordered tuple of dimensions of the shared modes between first and second, second and third... nodes is called the **rank**
of the TT-Tensor. Via truncated SVD decompositions these can be reduced. As this looses some precision in the representation of
the original tensor such a process is called **rounding**.

