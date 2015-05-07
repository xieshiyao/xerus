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
* [Building xerus](@ref md_tut_build) - instruction on how to build the library iteself and your first own program using it.
* [First Calculations](@ref md_tut_fullTensor) - defining and using tensors of arbitrary degree and shorthand index notations.
* [TT Tensors](@ref md_tut_ttTensor) - using the MPS or Tensor-Train decomposition.
* [Debugging](@ref md_tut_debugging) - using `xerus`'s capabilities to debug your own application.

## Issues

Should you have any problems with the library do not hesitate to contact us at [contact[at]libxerus.org](mailto:contact[at]libxerus.org) or describe your problem in the [issuetracker](https://git.hemio.de/xerus/xerus/issues).


## Licence and Attribution

The `xerus` library is published under the AGPL v3.0 public licence. With proper attribution you are thus allowed to freely use the source code but if you distribute your software including the `xerus`
library (or you provide a service based on this library) you have to provide the full source code under a compatible licence. 
For more information see the [AGPL licence](http://www.gnu.org/licenses/agpl-3.0.html).

If this library proved useful in your scientific research we would be grateful for contributions to the source code and / or a citation. As there is no related article yet, a corresponding bibtex 
entry might look as follows
~~~
@Misc{xerus,
	author = {Huber, Benjamin and Wolf, Sebastian},
	title = {Xerus - A General Purpose Tensor Library},
	howpublished = {\url{https://libxerus.org/}},
	year = {2014--2015}
}
~~~
