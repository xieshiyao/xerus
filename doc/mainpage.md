# xerus - a general purpose tensor library    {#mainpage}

## Introduction

The `xerus` library is a general purpose library for numerical calculations with higher order tensors, Tensor-Train Decompositions / Matrix Product States and general Tensor Networks.
The focus of development was the simple usability and adaptibility to any setting that requires higher order tensors or decompositions thereof. 

The key features include:
* Modern code and concepts incorporating many features of the new `C++11` standard.
* Calculation with tensors of arbitrary orders using an intuitive Einstein-like notation `A(i,j) = B(i,k,l) * C(k,j,l);`.
* Full implementation of the Tensor-Train decompositions (MPS) with all neccessary capabilities (including Algorithms like ALS, ADF and CG).
* Lazy evaluation of multiple tensor contractions featuring heuristics to find the most effective contraction order.
* Direct integration of the `blas` and `lapack`, as high performance linear algebra backends.
* Fast sparse tensor calculation by usage of the `suiteSparse` sparse matrix capabilities.
* Capabilites to handle arbitrary Tensor Networks.

## Version History

We released our first stable version 1.0 in May 2015 and are currently at version 2.2.1. It can be obtained via [git](https://git.hemio.de/xerus/xerus) or as an archived download via the same link.
The current development version is also available in the same git repository (branch 'development') and it might include more features than the latest stable release, but be adviced that these development versions are particularly bad documented and might change drastically from one day to the next.

## Getting Started

There are a number of tutorials to get you started using the `xerus` library.
* [Building xerus](@ref md_building_xerus) - instruction on how to build the library iteself and your first own program using it.
* [Quick-Start guide](_quick-_start-example.html) - a short introduction into the basic `xerus` functionality.
* [TT Tensors](_t_t-_tensors_01_07_m_p_s_08-example.html) - using the MPS or Tensor-Train decomposition.
* [Debugging](@ref md_tut_debugging) - using `xerus`'s capabilities to debug your own application.

## Issues

Should you have any problems with the library do not hesitate to contact us at [contact[at]libxerus.org](mailto:contact[at]libxerus.org) or describe your problem in the [issuetracker](https://git.hemio.de/xerus/xerus/issues).


## Licence and Attribution

The `xerus` library is published under the AGPL v3.0. With proper attribution you are thus allowed to freely use and modify the source code but if you distribute your software including the `xerus`
library (or you provide a service based on this library) you have to provide the full source code under a compatible licence. For more information see the [AGPL](http://www.gnu.org/licenses/agpl-3.0.html).

If this library proved useful in your scientific research we would be grateful for any contributions to the source code and citations. As there is no related article yet, a corresponding bibtex 
entry might look as follows
~~~
@Misc{xerus,
	author = {Huber, Benjamin and Wolf, Sebastian},
	title = {Xerus - A General Purpose Tensor Library},
	howpublished = {\url{https://libxerus.org/}},
	year = {2014--2017}
}
~~~
