// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
// 
// Xerus is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// Xerus is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with Xerus. If not, see <http://www.gnu.org/licenses/>.
//
// For further information on Xerus visit https://libXerus.org 
// or contact us at contact@libXerus.org.

/**
 * @file
 * @brief Definition of the python bindings.
 */

#ifdef XERUS_EXPERIMENTAL_PYTHON_WRAPPER

#include <boost/function.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include "xerus.h"

using namespace boost::python;



template<class t_type, class... args>
static xerus::internal::IndexedTensor<t_type>* indexing_wrapper(t_type &_this, args... _indices) {
	return new xerus::internal::IndexedTensor<t_type>(std::move(_this(_indices...)));
}


template<typename containedType>
struct custom_vector_to_list{
	static PyObject* convert(const std::vector<containedType>& v){
		list ret;
		for(const containedType& e : v) {
			ret.append(e);
		}
		return incref(ret.ptr());
	}
};
template<typename containedType>
struct custom_vector_from_seq{
	custom_vector_from_seq(){ converter::registry::push_back(&convertible,&construct,type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, converter::rvalue_from_python_stage1_data* data){
		void* storage=(reinterpret_cast<converter::rvalue_from_python_storage<std::vector<size_t> >*>(data))->storage.bytes;
		new (storage) std::vector<size_t>();
		std::vector<size_t>* v=reinterpret_cast<std::vector<size_t>*>(storage);
		Py_ssize_t l=PySequence_Size(obj_ptr); if(l<0) abort();
		v->reserve(size_t(l)); 
		for (Py_ssize_t i=0; i<l; i++) { 
			v->push_back(extract<size_t>(PySequence_GetItem(obj_ptr,i))); 
		}
		data->convertible=storage;
	}
};


// TODO get_copy() wrapper in python correct manner
BOOST_PYTHON_MODULE(libxerus) {
	using namespace xerus;
	
	bool show_user_defined = true;
	bool show_signatures = false;
	docstring_options doc_options(show_user_defined, show_signatures);
	
#define parametersDocstr "\n\nParameters\n----------\n"
#define	returnsDocstr "\n\nReturns\n-------\n"
	
	#define VECTOR_TO_PY(type, name) class_<std::vector<type>, boost::noncopyable>(name, no_init); \
		custom_vector_from_seq<type>(); \
		to_python_converter<std::vector<type>, custom_vector_to_list<type>>(); void(0)
	
	VECTOR_TO_PY(size_t, "IntegerVector");
	VECTOR_TO_PY(double, "DoubleVector");
	
	// --------------------------------------------------------------- index
	class_<Index>("Index",
		"helper class to define objects to be used in indexed expressions"
	)
		.def(init<int64_t>())
		.def("__pow__", &Index::operator^, "i**d changes the index i to span d indices in the current expression")
		.def("__xor__", &Index::operator^, "i^d changes the index i to span d indices in the current expression")
		.def("__div__", &Index::operator/, "i/n changes the index i to span 1/n of all the indices of the current object")
		.def("__and__", &Index::operator&, "i&d changes the index i to span all but d indices of the current object")
		.def("__str__", static_cast<std::string (*)(const Index &)>(&misc::to_string<Index>))
	;
	implicitly_convertible<int64_t, Index>();
	def("indices", +[](int _n)->boost::python::list{
		list res;
		for (int i = 0; i < _n; ++i) {
			res.append(Index());
		}
		return res;
	}, args("n"), "Creates n distinct indices."
		parametersDocstr "n : int"
		returnsDocstr "list of Index"
	);
	
	implicitly_convertible<internal::IndexedTensorReadOnly<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	implicitly_convertible<internal::IndexedTensorWritable<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	implicitly_convertible<internal::IndexedTensorMoveable<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	implicitly_convertible<internal::IndexedTensor<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	// NOTE in the following all __mul__ variants are defined for the ReadOnly indexed Tensors, even if they are meant for
	//		the moveable indexed tensors. boost will take care of the proper matching that way. if IndexedTensorMoveable
	//		defined an __mul__ function on its own it would overwrite all overloaded variants of the readonly indexed tensors
	//		and thus loose a lot of functionality.
	// ---------------------------------------------- indexedTensor<TN>
	using namespace internal;
#define ADD_MOVE_AND_RESULT_PTR(name, op, lhs_type, rhs_type, res_type) \
	.def(name, \
			+[](lhs_type &_l, rhs_type &_r) -> res_type* { \
				LOG(pydebug, "python wrapper: " name);\
				return new res_type(std::move(_l) op std::move(_r)); \
			}, return_value_policy<manage_new_object>())
	
	class_<internal::IndexedTensorReadOnly<TensorNetwork>, boost::noncopyable>("IndexedTensorNetworkReadOnly", no_init)
		ADD_MOVE_AND_RESULT_PTR("__add__", +, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__add__", +, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__sub__", -, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__sub__", -, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorMoveable<TensorNetwork>, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorMoveable<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<TensorNetwork>)
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, value_t _r) -> internal::IndexedTensorReadOnly<TensorNetwork>* {
				LOG(pydebug, "mul TN ro * scalar");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(_l) * _r);
			}, return_value_policy<manage_new_object>())
		.def("__rmul__", 
			+[](value_t _r, internal::IndexedTensorReadOnly<TensorNetwork> &_l) -> internal::IndexedTensorReadOnly<TensorNetwork>* {
				LOG(pydebug, "mul TN scalar * ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(_l) * _r);
			}, return_value_policy<manage_new_object>())
		.def("__div__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, value_t _r) -> internal::IndexedTensorReadOnly<TensorNetwork>* {
				LOG(pydebug, "div TN ro / scalar");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(_l) / _r);
			}, return_value_policy<manage_new_object>())
		.def("frob_norm", static_cast<value_t (*)(const IndexedTensorReadOnly<TensorNetwork> &)>(&frob_norm<TensorNetwork>))
		.def(float_(self)) // cast to double
	;
	
	class_<internal::IndexedTensorWritable<TensorNetwork>, boost::noncopyable, bases<internal::IndexedTensorReadOnly<TensorNetwork>>>("IndexedTensorNetworkWriteable", no_init)
	;
	class_<internal::IndexedTensorMoveable<TensorNetwork>, boost::noncopyable, bases<internal::IndexedTensorWritable<TensorNetwork>>>("IndexedTensorNetworkMoveable", no_init)
	;
	class_<internal::IndexedTensor<TensorNetwork>, boost::noncopyable, bases<internal::IndexedTensorWritable<TensorNetwork>>>("IndexedTensorNetwork", no_init)
		.def("__lshift__", 
			+[](internal::IndexedTensor<TensorNetwork> &_lhs, internal::IndexedTensorReadOnly<Tensor> &_rhs) {
				std::move(_lhs) = std::move(_rhs);
			})
		.def("__lshift__", 
			+[](internal::IndexedTensor<TensorNetwork> &_lhs, internal::IndexedTensorReadOnly<TensorNetwork> &_rhs) {
				std::move(_lhs) = std::move(_rhs);
			})
	;
	
	
	// --------------------------------------------- indexedTensor<Tensor>
	
	class_<internal::IndexedTensorReadOnly<Tensor>, boost::noncopyable>("IndexedTensorReadOnly", no_init)
		ADD_MOVE_AND_RESULT_PTR("__add__", +, IndexedTensorReadOnly<Tensor>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<Tensor>)
		ADD_MOVE_AND_RESULT_PTR("__sub__", -, IndexedTensorReadOnly<Tensor>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<Tensor>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorReadOnly<Tensor>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__mul__", *, IndexedTensorReadOnly<Tensor>, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
		ADD_MOVE_AND_RESULT_PTR("__div__", /, IndexedTensorReadOnly<Tensor>, IndexedTensorReadOnly<Tensor>, IndexedTensorMoveable<Tensor>)
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, value_t _r) -> internal::IndexedTensorReadOnly<Tensor>* {
				LOG(pydebug, "mul ro * scalar");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(_l) * _r);
			}, return_value_policy<manage_new_object>())
		.def("__rmul__", 
			+[](value_t _r, internal::IndexedTensorReadOnly<Tensor> &_l) -> internal::IndexedTensorReadOnly<Tensor>* {
				LOG(pydebug, "mul scalar * ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(_l) * _r);
			}, return_value_policy<manage_new_object>())
		.def("__div__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, value_t _r) -> internal::IndexedTensorReadOnly<Tensor>* {
				LOG(pydebug, "div ro / scalar");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(_l) / _r);
			}, return_value_policy<manage_new_object>())
		.def("frob_norm", static_cast<value_t (*)(const IndexedTensorReadOnly<Tensor> &)>(&frob_norm<Tensor>))
		.def(float_(self)) // cast to double
	;
	class_<internal::IndexedTensorWritable<Tensor>, boost::noncopyable, bases<internal::IndexedTensorReadOnly<Tensor>>>("IndexedTensorWriteable", no_init)
	;
	class_<internal::IndexedTensorMoveable<Tensor>, boost::noncopyable, bases<internal::IndexedTensorWritable<Tensor>>>("IndexedTensorMoveable", no_init)
	;
	class_<internal::IndexedTensor<Tensor>, boost::noncopyable, bases<internal::IndexedTensorWritable<Tensor>>>("IndexedTensor", no_init)
		.def("__lshift__", 
			+[](internal::IndexedTensor<Tensor> &_lhs, internal::IndexedTensorReadOnly<Tensor> &_rhs) {
				std::move(_lhs) = std::move(_rhs);
			})
		.def("__lshift__", 
			+[](internal::IndexedTensor<Tensor> &_lhs, internal::IndexedTensorReadOnly<TensorNetwork> &_rhs) {
				std::move(_lhs) = std::move(_rhs);
			})
	;
	
	
	// ----------------------------------------------------------- decompositions
	
	// ----------------------------------------------------------- Tensor
	{ scope Tensor_scope = 
		class_<Tensor>("Tensor",
			"a non-decomposed Tensor in either sparse or dense representation"
		)
			.def(init<const Tensor::DimensionTuple&>(args("dimensions"), "constructs a Tensor with the given dimensions"))
			.def(init<const TensorNetwork&>())
	// 		.def(init<const Tensor::DimensionTuple&, const boost::function<value_t()> &>()) // could define as argument boost::python::object _f, but maybe rather create a new converter?
			.add_property("dimensions", +[](Tensor &_A) {
				return _A.dimensions;
			})
			.def("degree", &Tensor::degree)
			.def_readonly("factor", &Tensor::factor)
			.def_readonly("size", &Tensor::size)
			.def("frob_norm", &Tensor::frob_norm)
			.def("random", 
				+[](std::vector<size_t> _dim) {
					static std::random_device rd;
					std::mt19937_64 rnd(rd());
					std::normal_distribution<double> dist(0.0, 1.0);
					return xerus::Tensor::random(_dim, rnd, dist);
				}).staticmethod("random")
			.def("ones", &Tensor::ones, args("dimensions"), 
				 "Constructs a Tensor of given dimensions that is equal to 1 everywhere."
				  parametersDocstr "dimensions : list or tuple of int"
			).staticmethod("ones")
			.def("identity", &Tensor::identity).staticmethod("identity")
			.def("kronecker", &Tensor::kronecker).staticmethod("kronecker")
			.def("dirac", static_cast<Tensor (*)(Tensor::DimensionTuple, const Tensor::MultiIndex&)>(&Tensor::dirac))
			.def("dirac", static_cast<Tensor (*)(Tensor::DimensionTuple, const size_t)>(&Tensor::dirac)).staticmethod("dirac")
			.def("has_factor", &Tensor::has_factor)
			.def("is_dense", &Tensor::is_dense)
			.def("is_sparse", &Tensor::is_sparse)
			.def("sparsity", &Tensor::sparsity)
	// 		.def("reinterpret_dimensions", &Tensor::reinterpret_dimensions)
			.def("__call__", &indexing_wrapper<Tensor>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<Tensor,Index, Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__str__", &Tensor::to_string)
			.def(self * other<value_t>())
			.def(other<value_t>() * self)
			.def(self / other<value_t>())
			.def(self + self)
			.def(self - self)
			.def("__getitem__", +[](Tensor &_this, size_t _i) {
				// TODO Note: for loops expect that an IndexError will be raised for illegal indexes to allow proper detection of the end of the sequence.
				return _this[_i];
			})
			.def("__getitem__", +[](Tensor &_this, std::vector<size_t> _idx) {
				return _this[_idx];
			})
			.def("__setitem__", +[](Tensor &_this, size_t _i, value_t _val) {
				_this[_i] = _val;
			})
			.def("__setitem__", +[](Tensor &_this, std::vector<size_t> _i, value_t _val) {
				_this[_i] = _val;
			})
		;
		enum_<Tensor::Representation>("Representation", "Possible representations of Tensor objects.")
			.value("Dense", Tensor::Representation::Dense)
			.value("Sparse", Tensor::Representation::Sparse)
// 			.export_values() // would define Tensor.Sparse = Tensor.Representation.Sparse etc.
		;
	} // close Tensor_scope
	
	
	// ------------------------------------------------------------- TensorNetwork
	VECTOR_TO_PY(TensorNetwork::Link, "LinkVector");
	VECTOR_TO_PY(TensorNetwork::Node, "TensorNetworkNodeVector");
	{ scope TN_scope =
		class_<TensorNetwork>("TensorNetwork")
			.def(init<Tensor>())
			.add_property("dimensions", +[](TensorNetwork &_A) {
				return _A.dimensions;
			})
			.def("degree", &TensorNetwork::degree)
			.def("datasize", &TensorNetwork::datasize)
			.add_property("nodes", +[](TensorNetwork &_this){
				return _this.nodes;
			})
			.def("node", +[](TensorNetwork &_this, size_t _i) {
				return _this.nodes[_i];
			})
			.add_property("externalLinks", +[](TensorNetwork &_this){
				return _this.externalLinks;
			})
			.def("__call__", &indexing_wrapper<TensorNetwork>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def("__call__", &indexing_wrapper<TensorNetwork,Index, Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
			.def(self * other<value_t>())
			.def(other<value_t>() * self)
			.def(self / other<value_t>())
			.def("__getitem__", +[](TensorNetwork &_this, size_t _i) {
				// TODO Note: for loops expect that an IndexError will be raised for illegal indexes to allow proper detection of the end of the sequence.
				return _this[_i];
			})
			.def("__getitem__", +[](TensorNetwork &_this, std::vector<size_t> _idx) {
				return _this[_idx];
			})
// 			.def("reshuffle_nodes", +[](TensorNetwork &_this, boost::python::object _f) { //TODO
// 				_this.reshuffle_nodes(_f);
// 			})
			.def("require_valid_network", +[](TensorNetwork &_this) {
				_this.require_valid_network();
			})
			.def("require_correct_format", &TensorNetwork::require_correct_format)
			.def("swap_external_links", &TensorNetwork::swap_external_links)
			.def("round_edge", &TensorNetwork::round_edge)
			.def("transfer_core", &TensorNetwork::transfer_core)
			.def("transfer_core", +[](TensorNetwork &_this, size_t _from, size_t _to) {
				_this.transfer_core(_from, _to);
			})
			.def("reduce_representation", &TensorNetwork::reduce_representation)
			//TODO find_common_edge (wrapper that returns python::tuple
			.def("sanitize", &TensorNetwork::sanitize)
			.def("fix_slate", &TensorNetwork::fix_slate) //TODO rename
			.def("remove_slate", &TensorNetwork::remove_slate)
			.def("resize_dimension", &TensorNetwork::resize_dimension)
			.def("resize_dimension", +[](TensorNetwork &_this, size_t _mode,  size_t _newDim) {
				_this.resize_dimension(_mode, _newDim);
			})
			.def("contract", static_cast<void (TensorNetwork::*)(size_t, size_t)>(&TensorNetwork::contract))
			.def("contract", static_cast<size_t (TensorNetwork::*)(const std::set<size_t>&)>(&TensorNetwork::contract)) //TODO write converter
			.def("contraction_cost", &TensorNetwork::contraction_cost)
			.def("draw", &TensorNetwork::draw)
			.def("frob_norm", &TensorNetwork::frob_norm)
		;
		
		class_<TensorNetwork::TensorNode>("TensorNode")
			.def("size", &TensorNetwork::TensorNode::size)
			.def("degree", &TensorNetwork::TensorNode::degree)
// 			.def_readonly("erased", &TensorNetwork::TensorNode::erased) // internal
// 			.def("erase", &TensorNetwork::TensorNode::erase) // internal
			.def_readonly("tensorObject", &TensorNetwork::TensorNode::tensorObject) // TODO should return an optional
			.add_property("neighbors", +[](TensorNetwork::TensorNode &_this){
				return _this.neighbors;
			})
		;
		
		class_<TensorNetwork::Link>("TensorNetworkLink")
			.def_readonly("other", &TensorNetwork::Link::other)
			.def_readonly("indexPosition", &TensorNetwork::Link::indexPosition)
			.def_readonly("dimension", &TensorNetwork::Link::dimension)
			.def_readonly("external", &TensorNetwork::Link::external)
			.def("links", &TensorNetwork::Link::links)
		;
	} // closes TN_scope

	
	// ------------------------------------------------------------- TTNetwork	
	class_<TTTensor, bases<TensorNetwork>>("TTTensor")
		.def(init<const Tensor&, optional<value_t, size_t>>())
		.def(init<const Tensor&, value_t, TensorNetwork::RankTuple>())
		.def(init<Tensor::DimensionTuple>())
		.def(init<size_t>())
		.def("get_component", &TTTensor::get_component, return_value_policy<copy_const_reference>())
		.def("set_component", &TTTensor::set_component)
		.def_readonly("cannonicalized", &TTTensor::cannonicalized)
		.def_readonly("corePosition", &TTTensor::corePosition)
		.def("ranks", &TTTensor::ranks)
		.def("rank", &TTTensor::rank)
// 		.def("frob_norm", &TTTensor::frob_norm) // NOTE unneccessary because correct call is inherited
		.def("random", 
			+[](std::vector<size_t> _dim, std::vector<size_t> _rank) {
				static std::random_device rd;
				std::mt19937_64 rnd(rd());
				std::normal_distribution<double> dist(0.0, 1.0);
				return xerus::TTTensor::random(_dim, _rank, rnd, dist);
			}).staticmethod("random")
		.def("ones", &TTTensor::ones).staticmethod("ones")
// 		.def("identity", &TTTensor::identity).staticmethod("identity") // exists for  TTOperator only
// 		.def("kronecker", &TTTensor::kronecker).staticmethod("kronecker") //TODO
// 		.def("dirac", static_cast<TTTensor (*)(Tensor::DimensionTuple, const Tensor::MultiIndex&)>(&TTTensor::dirac)) //TODO
// 		.def("dirac", static_cast<TTTensor (*)(Tensor::DimensionTuple, const size_t)>(&TTTensor::dirac)).staticmethod("dirac") //TODO
		
		.def("reduce_to_maximal_ranks", &TTTensor::reduce_to_maximal_ranks).staticmethod("reduce_to_maximal_ranks")
		.def("dyadic_product", static_cast<TTTensor (*)(const std::vector<TTTensor> &)>(&TTTensor::dyadic_product))
			.staticmethod("dyadic_product") //  TODO should not be in TTTensor scope
		.def("entrywise_square", &TTTensor::entrywise_square) // TODO should not be a member method?
		.def("find_largest_entry", &TTTensor::find_largest_entry) // TODO should not be a member method!
		// TODO chop wrapper
		
		.def("round", +[](TTTensor &_this, std::vector<size_t> _rank){
			_this.round(_rank);
		})
		.def("round", static_cast<void (TTTensor::*)(const std::vector<size_t>&, double)>(&TTTensor::round))
		.def("round", static_cast<void (TTTensor::*)(double)>(&TTTensor::round))
		.def("round", static_cast<void (TTTensor::*)(size_t)>(&TTTensor::round))
		
		.def("soft_threshold", +[](TTTensor &_this, double _tau) {
			_this.soft_threshold(_tau);
		})
		.def("soft_threshold", +[](TTTensor &_this, std::vector<double> _taus) {
			_this.soft_threshold(_taus);
		})
		.def("soft_threshold", static_cast<void (TTTensor::*)(const double, const bool)>(&TTTensor::soft_threshold))
		.def("soft_threshold", static_cast<void (TTTensor::*)(const std::vector<double>&, const bool)>(&TTTensor::soft_threshold))
		
		.def("move_core", +[](TTTensor &_this, size_t _pos){
			_this.move_core(_pos);
		})
		.def("move_core", &TTTensor::move_core)
		
		.def("assume_core_position", &TTTensor::assume_core_position)
		.def("cannonicalize_left", &TTTensor::cannonicalize_left)
		.def("cannonicalize_right", &TTTensor::cannonicalize_right)
		.def(self + self)
		.def(self - self)
		.def(self * other<value_t>())
		.def(other<value_t>() * self)
		.def(self / other<value_t>())
	;
	def("entrywise_product", static_cast<TTTensor (*)(const TTTensor&, const TTTensor&)>(&entrywise_product));
	
	class_<TTOperator, bases<TensorNetwork>>("TTOperator")
		.def(init<const Tensor&, optional<value_t, size_t>>())
		.def(init<const Tensor&, value_t, TensorNetwork::RankTuple>())
		.def(init<Tensor::DimensionTuple>())
		.def(init<size_t>())
		.def("get_component", &TTOperator::get_component, return_value_policy<copy_const_reference>())
		.def("set_component", &TTOperator::set_component)
		.def_readonly("cannonicalized", &TTOperator::cannonicalized)
		.def_readonly("corePosition", &TTOperator::corePosition)
		.def("ranks", &TTOperator::ranks)
		.def("rank", &TTOperator::rank)
// 		.def("frob_norm", &TTOperator::frob_norm) // NOTE unneccessary because correct call is inherited
		.def("random", 
			+[](std::vector<size_t> _dim, std::vector<size_t> _rank) {
				static std::random_device rd;
				std::mt19937_64 rnd(rd());
				std::normal_distribution<double> dist(0.0, 1.0);
				return xerus::TTOperator::random(_dim, _rank, rnd, dist);
			}).staticmethod("random")
		.def("ones", &TTOperator::ones).staticmethod("ones")
// 		.def("kronecker", &TTOperator::kronecker).staticmethod("kronecker") //TODO
// 		.def("dirac", static_cast<TTOperator (*)(Tensor::DimensionTuple, const Tensor::MultiIndex&)>(&TTOperator::dirac)) //TODO
// 		.def("dirac", static_cast<TTOperator (*)(Tensor::DimensionTuple, const size_t)>(&TTOperator::dirac)).staticmethod("dirac") //TODO
		
		.def("reduce_to_maximal_ranks", &TTOperator::reduce_to_maximal_ranks).staticmethod("reduce_to_maximal_ranks")
		.def("dyadic_product", static_cast<TTOperator (*)(const std::vector<TTOperator> &)>(&TTOperator::dyadic_product))
			.staticmethod("dyadic_product") //  TODO should not be in TTOperator scope
		.def("entrywise_square", &TTOperator::entrywise_square) // TODO should not be a member method?
		.def("find_largest_entry", &TTOperator::find_largest_entry) // TODO should not be a member method!
		// TODO chop wrapper
		
		.def("round", +[](TTOperator &_this, std::vector<size_t> _rank){
			_this.round(_rank);
		})
		.def("round", static_cast<void (TTOperator::*)(const std::vector<size_t>&, double)>(&TTOperator::round))
		.def("round", static_cast<void (TTOperator::*)(double)>(&TTOperator::round))
		.def("round", static_cast<void (TTOperator::*)(size_t)>(&TTOperator::round))
		
		.def("soft_threshold", +[](TTOperator &_this, double _tau) {
			_this.soft_threshold(_tau);
		})
		.def("soft_threshold", +[](TTOperator &_this, std::vector<double> _taus) {
			_this.soft_threshold(_taus);
		})
		.def("soft_threshold", static_cast<void (TTOperator::*)(const double, const bool)>(&TTOperator::soft_threshold))
		.def("soft_threshold", static_cast<void (TTOperator::*)(const std::vector<double>&, const bool)>(&TTOperator::soft_threshold))
		
		.def("move_core", +[](TTOperator &_this, size_t _pos){
			_this.move_core(_pos);
		})
		.def("move_core", &TTOperator::move_core)
		
		.def("assume_core_position", &TTOperator::assume_core_position)
		.def("cannonicalize_left", &TTOperator::cannonicalize_left)
		.def("cannonicalize_right", &TTOperator::cannonicalize_right)
		.def(self + self)
		.def(self - self)
		.def(self * other<value_t>())
		.def(other<value_t>() * self)
		.def(self / other<value_t>())
		
		
		// for  TTOperator only:
		.def("identity", &TTOperator::identity<true>).staticmethod("identity") 
		.def("transpose", &TTOperator::transpose<true>)
	;
	def("entrywise_product", static_cast<TTOperator (*)(const TTOperator&, const TTOperator&)>(&entrywise_product));
	
	// ------------------------------------------------------------- Algorithms
	VECTOR_TO_PY(PerformanceData::DataPoint, "PerfDataVector");
	{ scope perfdata =
		class_<PerformanceData>("PerformanceData")
			.def_readwrite("printProgress", &PerformanceData::printProgress)
			.def_readwrite("startTime", &PerformanceData::startTime)
			.def_readwrite("stopTime", &PerformanceData::stopTime)
			.def_readwrite("additionalInformation", &PerformanceData::additionalInformation)
			.add_property("data", +[](PerformanceData &_this){
				return _this.data;
			}, +[](PerformanceData &_this, std::vector<PerformanceData::DataPoint> &_newData){
				_this.data = _newData;
			}, )
			// TODO convert the errFunc
			.add_property("errorFunction", +[](PerformanceData &_this){
				return _this.errorFunction;
			}, +[](PerformanceData &_this, boost::python::object _f){
				_this.errorFunction = _f;
			}, )
			.def(init<bool>())
			.def("start", &PerformanceData::start)
			.def("stop_timer", &PerformanceData::stop_timer)
			.def("continue_timer", &PerformanceData::continue_timer)
			.def("reset", &PerformanceData::reset)
			.def("get_elapsed_time", &PerformanceData::get_elapsed_time)
			.def("get_runtime", &PerformanceData::get_runtime)
			.def("add", +[](PerformanceData &_this, size_t _itr,  value_t _res){
				_this.add(_itr, _res);
			})
			.def("add", +[](PerformanceData &_this, size_t _itr,  value_t _res, const TTTensor &_x){
				_this.add(_itr, _res, _x);
			})
			.def("add", +[](PerformanceData &_this, size_t _itr,  value_t _res, const TTTensor &_x, size_t _flags){
				_this.add(_itr, _res, _x, _flags);
			})
			.def("add", +[](PerformanceData &_this, value_t _res, const TTTensor &_x, size_t _flags){
				_this.add(_res, _x, _flags);
			})
			.def("__nonzero__", +[](PerformanceData &_this){ return bool(_this); })
			.def("dump_to_file", &PerformanceData::dump_to_file)
			// TODO histogram
		;
		
		class_<PerformanceData::DataPoint>("DataPoint", no_init)
			.def_readonly("iterationCount", &PerformanceData::DataPoint::iterationCount)
			.def_readonly("elapsedTime", &PerformanceData::DataPoint::elapsedTime)
			.def_readonly("residual", &PerformanceData::DataPoint::residual)
			.def_readonly("error", &PerformanceData::DataPoint::error)
			.def_readonly("ranks", &PerformanceData::DataPoint::ranks)
			.def_readonly("flags", &PerformanceData::DataPoint::flags)
		;
	}
	
	class_<ALSVariant>("ALSVariant", init<uint, size_t, ALSVariant::LocalSolver, bool, optional<bool>>())
		.def_readwrite("sites", &ALSVariant::sites)
		.def_readwrite("numHalfSweeps", &ALSVariant::numHalfSweeps)
		.def_readwrite("convergenceEpsilon", &ALSVariant::convergenceEpsilon)
		.def_readwrite("useResidualForEndCriterion", &ALSVariant::useResidualForEndCriterion)
		.def_readwrite("preserveCorePosition", &ALSVariant::preserveCorePosition)
		.def_readwrite("assumeSPD", &ALSVariant::assumeSPD)
		.def_readwrite("localSolver", &ALSVariant::localSolver) // TODO we need some wrapper here
		
		// TODO add variants with perfdata
		// TODO test the following!
		// NOTE it is functionally enough to export the (A,x,b) and (x,b) variant without _eps or _numHalfSweeps
		//      as it is expected that there will be problems with ambiguity in python we thus leave out the _eps / _numHalfSweeps variants
		.def("__call__", +[](ALSVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b) {
			_this(_A, _x, _b);
		})
		.def("__call__", +[](ALSVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _eps) {
			_this(_A, _x, _b, _eps);
		})
		.def("__call__", +[](ALSVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps) {
			_this(_A, _x, _b, _numHalfSweeps);
		})
		.def("__call__", +[](ALSVariant &_this, TTTensor &_x, const TTTensor &_b) {
			_this(_x, _b);
		})
		.def("__call__", +[](ALSVariant &_this, TTTensor &_x, const TTTensor &_b, value_t _eps) {
			_this(_x, _b, _eps);
		})
		.def("__call__", +[](ALSVariant &_this, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps) {
			_this(_x, _b, _numHalfSweeps);
		})
	;
// 	add_property("ALS", +[](){ return ALS; });
	def("ALS", +[](TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
			ALS(_x, _b, _pd);
		})
	// ------------------------------------------------------------- misc
	def("approx_equal", static_cast<bool (*)(const TensorNetwork&, const TensorNetwork&, double)>(&approx_equal));
	def("approx_equal", static_cast<bool (*)(const Tensor&, const TensorNetwork&, double)>(&approx_equal));
	def("approx_equal", static_cast<bool (*)(const TensorNetwork&, const Tensor&, double)>(&approx_equal));
	def("approx_equal", static_cast<bool (*)(const Tensor&, const Tensor&, double)>(&approx_equal));
	def("approx_equal", +[](const Tensor& _l, const Tensor& _r) {
		return approx_equal(_l, _r);
	});
	def("approx_equal", +[](const Tensor& _l, const TensorNetwork& _r) {
		return approx_equal(_l, _r);
	});
	def("approx_equal", +[](const TensorNetwork& _l, const TensorNetwork& _r) {
		return approx_equal(_l, _r);
	});
	def("approx_equal", +[](const TensorNetwork& _l, const TensorNetwork& _r) {
		return approx_equal(_l, _r);
	});
	
	def("log", +[](std::string _msg){
		LOG_SHORT(info, _msg);
	});
	
	// TODO streamwriter stuff
	
	// identity returns the cpp name to a python object
// 	def("identity", identity_);
	
	// the following is probably not necessary because generic_error inherits from std::exception
// 	register_exception_translator<misc::generic_error>([](const misc::generic_error &_e){
// 		LOG(pydebug, "custom exception handler called with " << _e.what());
// 		PyErr_SetString(PyExc_UserWarning, _e.what());
// 	});
}

#endif // XERUS_EXPERIMENTAL_PYTHON_WRAPPER
