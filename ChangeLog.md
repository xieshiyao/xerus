## Changelog

Potentially breaking changes are marked with an exclamation point '!' at the begin of their description.

* 2015-10-14 v1.4.1
 * Fixed a bug with the recannonicalization of products of TTOperators.
 * Many improvements for the ADF algorithm for tensor completion.
 * Added TTTensor::ones and TTOperator::ones to create TTTensors with all entries equal to 1.

* 2015-09-09 v1.4.0
 * Added an implementation of the ADF algorithm.
 * Added a custom allocator to better handle the small allocations of Indices etc. (requires to set the corresponding option in the config.mk). This allocator is highly recommended for use with the ADF algorithm.
 * Reduced overhead (time not spent in blas) of TensorNetwork contractions.
 * Tensor decomposition now allow arbitrary index orders in the target as well as the source.
 * Added TN::measure(...) to perform a batch of evaluations of the network.
 * Added TN::fix_slate(...) to fix a single index to a specific value.
 * Added TN::reduce_representation() to trace out any full-rank edge.
 * Added TN::transfer_core(...) and TN::round_edge(...) analogously to TT rounding / core movement operations.
 * ! renamed class::construct_random to class::random and TTOperator::construct_identity to TTOperator::identity.
 * Added several named constructors:
  * Tensor::ones(...) to construct tensors filled with all entries equal to 1.
  * Tensor::identity(...) to construct a (sparse) tensor equal to the TTOperator::identity.
  * Tensor::dirac(...) to construct a (sparse) tensor with a single 1 at a given position and 0 everywhere else.
  * Tensor::kronecker(...) to construct a tensor that is 1 when all indices are equal and 0 otherwise.
 * ! misc::approx_equal now uses a relative difference (for consistency with the tensor variant xerus::approx_equal).
 * Several bugfixes for SparseTensors.
 * Several minor bugfixes (TN::draw no longer uses temporary files, span-zero indices in assignment, ...).

* 2015-07-10 v1.3.0
 * Added several retractions to allow Riemannian algorithms on TT manifolds.
 * Added the basic (Riemannian) steepest descent algorithm.
 * Added the (so far untested) CG algorithm for TT tensors.
 * ! New PerformanceData class to log stepcount, time and residual of algorithm runs.
 * New algorithm to determine the largest entry of a TT tensor.
 * Fully rewrote the contraction heuristics, now including a local search for more optimal contraction orders.
 * Several bugfixes for order-zero tensors, TT tensors of more than full rank and other rare cases.

* 2015-06-15 v1.2.1
 * entrywise_product for FullTensor and TTNetwork.
 * TensorNetworks now without extraneous factor.

* 2015-06-15 v1.2.0
 * Rank revealing QR decomposition in TTNetwork:move_core. 
 * ALS performance improvement by skipping full-rank nodes.
 * Added TensorNetwork::draw to get a graphical representation of a given tensor network.
 * Added algorithm DecompositionALS to approximate a FullTensor by a TTTensor.
 * New Quick-Start guide.
 * Several bugfixes, more consistent operator availability and much more documentation.

* 2015-06-02 v1.1.1
 * Complete documentation for Tensor, FullTensor and SparseTensor
 * Added implementation of the projection ALS
 * Fixed performance bug in logging system
 * ALS now uses energy functional for convergence testing
 * Several minor fixes

* 2015-06-01 v1.1.0
 * ! New representation for the TT tensors: all component tensors are now of degree 3 / 4 (operators).
 * TT::get_component and ::set_component functions that handle the meta data.
 * More strict ISO c++ conformity: no names starting with '__'.
 * More elaborate performance profiling.

* 2015-05-27 v1.0.2
 * Added support for Clang++
 * ! Cleanup of namespace structure
 * Several documentation fixes

* 2015-05-20 v1.0.1
 * fixed a major bug with the handling of shared resources in TensorNetworks
 * added TTTensor/TTOperator::chop() function
 * some minor fixes

* 2015-05-19 v1.0
 * first stable release
