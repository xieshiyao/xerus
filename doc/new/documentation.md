---
layout: post
section: "Documentation"
---

# Documentation

This is the semi-complete documentation of the `xerus` library. It does not provide you with precise function declarations or
class hierarchies (check out the [doxygen documentation](doxygen) for those) but instead focuses on small working code snippets
to demonstrate `xerus`'s capabilities. 

The documentation is categorized into chapters that each focus on individual features of the library. In particular the chapters in the "Basic
Usage" section only build upon the knowledge of the previous chapters so that they can easily be read in order. They are still
useful as lookup resources later on though.

If you have not already done so, you will most likely want to start by downloading and [building xerus](building_xerus). If
you are uncertain or confused by our nomenclature at any point you can reference the [nomenclature](nomenclature) chapter which
should hopefully allow you to deduce the meaning of all tensor specific terms we use.

The "Basic Usage" section starts out with the creation and modification of (sparse or dense) tensors ([The Tensor Class](tensor)),
how to use them in indexed equations to denote contractions, summations etc. ([Indices and Equations](indices)) and finally 
how to denote generalizations of matrix decompositions ([Decompositions and Solve](decompositions)) before explainig any more
elaborate tensor format. The focus of `xerus` so far clearly lies on the Tensor Train decomposition ([TT-Tensors](tttensors)) and
algorithms for those to solve least squares problems ([Alternating Algorithms](als), [Riemannian Algorithms](riemannian)) or
tensor recovery and completion problems ([Tensor Completion / Recovery](completion)). It is possible to use `xerus`'s capabilities
to construct and contract arbitrary tensor networks though ([General Tensor Networks](tensornetworks)).

In the "Advanced Usage" section you find instructions on how to optimize your usage of the `xerus` library to gain those last
10-20% of speedup ([Optimization](optimization)) and explanations of the debugging tools that `xerus` gives you ([Debugging](debugging)).
