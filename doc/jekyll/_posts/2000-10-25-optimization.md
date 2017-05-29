---
layout: post
title: "Optimization"
date: 2000-10-25 16:25:06 +0000
topic: "Advanced Usage"
section: "Documentation"
---


# Optimizations

If you are like us, you want to get the fastest possible version of your numerical code to run as many samples as possible and
solve the largest systems possible. To this end there are a number of possible optimizations already provided for you by the 
`xerus` library. The following list expands on the most relevant of them in roughly the order of effectiveness.

## Disabling Runtime Checks
The library contains many runtime checks for out-of-bounds access, other invalid inputs (like illegal contractions), consistency
and even to check the correct behaviour of internal structures. Depending on the complexity your code and the time spent inside
`xerus` (and not one of the libraries it uses) you can expect a large performance gain by disabling these checks in the `config.mk`
file during compilation of xerus.

It is not advisable to do this while developing, as it will be much more difficult to detect errors in your calls to `xerus`
functions, but once you have established, that your code works as expected you might want to try replacing the `libxerus.so` object
used by your project with one compiled with the `-D XERUS_DISABLE_RUNTIME_CHECKS` flag.

## Use c++ instead of Python
The interface between the languages python and c++ makes it necessary to perform operations for the sole purpose of compatibility
between the otherwise incompatible languages. Often this includes copies of vectors of integers (whenever dimensions are specified
or queried) and sometimes even deep copies of whole tensors (`.from_ndarray()` and `.to_ndarray()`). The only way to get rid of
this overhead is to write your appliation in c++ instead of python. Most instructions that `xerus` offers for python look very
similar in c++, so a transition might be simpler than you think. Simply check out the rest of the tutorials to compare the code
snippets.

This transition is particularly useful, if you wrote your own numerical algorithms in python. As an example consider the simple
ALS implementation in the example section ([ALS](minimal_als))), where the c++ implementation is faster by about a factor of two.
If most of the runtime is spend inside one of `xerus`'s own algorithms like the `ALS`, it is likely not worth much though.

## Compiling Xerus with High Optimizations
Per default the library already compiles with high optimization settings (corresponding basically to `-O3`) as there is rarely
any reason to use lower settings for numerical code. If you want to spend a significant amount of cpu hours in numerical code
using the `xerus` library though, it might be worthwile to go even further.

The most significant change in runtime speed gains due to compiler settings at this point will come from link-time optimizations
(for `c++`projects using `xerus`).
To make use of them you will need a sufficiently recent versions of the `g++` compiler and `ar` archiver. After compiling the
`libxerus.so` object with the `USE_LTO = TRUE` flag you can then enable `-flto` in your own compilation process. The optimizations
that will be used then extending more than a single compilation unit and might thus use significant system resources during 
compilation.

If link-time optimization is not an option (or not sufficient) it is also possible to replace the high optimizations flag in your
`config.mk` file with the `DANGEROUS_OPTIMIZATION = TRUE` flag. This will enable non-IEEE-conform optimizations that should
typically only change floating point results in the least significant bit but might lead to undefined behaviour in case a `NaN`
or overflow is encountered during runtime. (It is rumored that there is an even higher optimization setting available for `xerus`
for those who know how to find it and want to get even the last 1% of speedup...)


## Avoiding Indexed Expressions
The comfort of being able to write Einstein-notation-like equations in the source code of the form `A(i,k) = B(i,j)*C(j,k);` 
comes with the price of a certain overhead during runtime. It is in the low single-digit percent range for typical applications
but can become significant when very small tensors are being used and the time for the actual contraction thus becomes negligible.

In such cases it can be useful to replace such equations (especially ones that are as simple as above) with the explicit statement 
of contractions and reshuffels. For above equation that would simply be
~~~ cpp
// equivalent to A(i,k) = B(i,j^2)*C(j^2,k)
contract(A, B, false, C, false, 2);
~~~
i.e. read as: contract two tensors and store the result in A, left hand side B, not transposed, right hand side C, not transposed, contract two modes.

If it is necessary to reshuffle a tensor to be able to contract it in such a way. This can be done
with the `reshuffle` function.
~~~ cpp
// equivalent to: A(i,j,k) = B(i,k,j)
reshuffle(A, B, {0,2,1});
~~~

Decompositions similarly have their low(er) level calls. They require properly reshuffled tensors and you have to provide a
`splitPosition`, i.e. the number of modes that will be represented by the left-hand-side of the result.
~~~ cpp
// equivalent to: (Q(i,j,r), R(r, k)) = xerus.QR(A(i,j,k))
calculate_qr(Q, R, A, 2);
~~~

It is our opinion that code written with these functions instead of indexed espressions are often much harder to understand
and the speedup is typically small... but just in case you really want to, you now have the option to use them.

