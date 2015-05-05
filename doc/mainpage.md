# xerus - a general purpose tensor library    {#mainpage}

## Introduction

The `xerus` library is a general purpose library for numerical calculations with higher order tensors, Tensor-Train Decompositions / Matrix Product States and other Tensor Networks.
The focus of development was the simple usability and adaptibility to any setting that requires higher order tensors or decompositions thereof. It thus achieves a level of 
adaptibility that is unmatched by any other open-source tensor library available.

Some key features include:
* Modern code and concepts incorporating many features of the new `c++11` standard.
* Calculation with tensors of arbitrary orders using an intuitive Einstein-like notation `A(i,j) = B(i,k,l) * C(k,j,l);`.
* Full implementation of the Tensor-Train decompositions (MPS) with all neccessary capabilities (including Algorithms like ALS and DMRG).
* Capabilites to handle arbitrary Tensor Networks and the possibility for the user to define own decomposition formats.
* Lazy evaluation of contractions featuring heuristics to find good contraction orders.
* Fast calculations thanks to wide-spread backends on the `blas` basis (in particular `openblas` or `atlas`) and the `suiteSparse` sparse capabilities.
* Thread-safe up to 256 threads.

## Version History

We have just released our first stable version 1.0 in May 2015. It can be obtained via [git](https://git.hemio.de/xerus/xerus/tree/master) or as an archived download via the same link.
The current development version is also available in the same git repository (branch 'development') and it might include some more features than the latest stable release (especially
at these early stages of development), but be adviced that these development versions are particularly bad documented and might change drastically from one day to the next.

## Getting Started

There are a number of tutorials to get you started using the `xerus` library.
* [Building xerus](@ref tutorial_build) - instruction on how to build the library iteself and your first own program using it.
* [First Calculations](@ref md_tut_fullTensor) - defining and using tensors of arbitrary degree and shorthand index notations.
* [TT Tensors](@ref md_tut_ttTensor) - using the MPS or Tensor-Train decomposition.
* [Debugging](@ref md_tut_debugging) - using `xerus`'s capabilities to debug your own application.

## Issues

Should you have any problems with the library do not hesitate to contact us at [contact[at]libxerus.org] or describe your problem in the [issuetracker](https://git.hemio.de/xerus/xerus/issues).

