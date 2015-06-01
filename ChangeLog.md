## Changelog

* 2015-06-01 v1.1.0
 * New representation for the TT tensors: all component tensors are now of degree 3 / 4 (operators).
 * TT::get_component and ::set_component functions that handle the meta data.
 * More strict ISO c++ conformity: no names starting with '__'.
 * More elaborate performance profiling.

* 2015-05-27 v1.0.2
 * Added support for Clang++
 * Cleanup of namespace structure
 * Several documentation fixes

* 2015-05-20 v1.0.1
 * fixed a major bug with the handling of shared resources in TensorNetworks
 * added TTTensor/TTOperator::chop() function
 * some minor fixes

* 2015-05-19 v1.0
 * first stable release
