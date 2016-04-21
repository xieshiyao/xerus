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

BOOST_PYTHON_MODULE(libxerus) {
	using namespace xerus;
	
	class_<std::vector<size_t>, boost::noncopyable>("IntegerVector", no_init); // just to give it a human readable name in python
	class_<std::vector<double>, boost::noncopyable>("DoubleVector", no_init);
	custom_vector_from_seq<size_t>();
	to_python_converter<std::vector<size_t>, custom_vector_to_list<size_t>>();
	custom_vector_from_seq<double>();
	to_python_converter<std::vector<double>, custom_vector_to_list<double>>();
	
	// --------------------------------------------------------------- index
	class_<Index>("Index")
		.def(init<int64_t>())
		.def("__pow__", &Index::operator^)
		.def("__xor__", &Index::operator^)
		.def("__div__", &Index::operator/)
		.def("__and__", &Index::operator&)
		.def("__str__", static_cast<std::string (*)(const Index &)>(&misc::to_string<Index>))
	;
	implicitly_convertible<int64_t, Index>();
	def("indices", +[](int _n)->boost::python::list{
		list res;
		for (int i = 0; i < _n; ++i) {
			res.append(Index());
		}
		return res;
	});
	
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
	class_<Tensor>("Tensor")
		.def(init<const Tensor::DimensionTuple&>())
		.def(init<const TensorNetwork&>())
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
		.def("ones", &Tensor::ones).staticmethod("ones")
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
	
	
	// ------------------------------------------------------------- TTNetwork
	class_<TensorNetwork>("TensorNetwork")
		.def(init<Tensor>())
		.add_property("dimensions", +[](TensorNetwork &_A) {
			return _A.dimensions;
		})
		.def("degree", &TensorNetwork::degree)
		.def("datasize", &TensorNetwork::datasize)
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
		.def("frob_norm", &TensorNetwork::frob_norm)
	;
	
	class_<TTTensor, bases<TensorNetwork>>("TTTensor")
		.def(init<const Tensor&>())
		.def(init<const Tensor&, value_t>())
		.def(init<const Tensor&, value_t, size_t>())
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
	
	// ------------------------------------------------------------- Algorithms
// 	def("ALS", +[](){
// 		
// 	});
// 	def("ALS_SPD", &ALS_SPD.operator());
// 	def("DMRG", &DMRG.operator());
// 	def("DMRG_SPD", &DMRG_SPD.operator());
// 	def("ASD", &ASD.operator());
// 	def("ASD_SPD", &ASD_SPD.operator());
	
	// ------------------------------------------------------------- misc
	def("log", +[](std::string _msg){
		LOG_SHORT(info, _msg);
	});
	
	// the following is probably not necessary because we inherit from std::exception
// 	register_exception_translator<misc::generic_error>([](const misc::generic_error &_e){
// 		LOG(pydebug, "custom exception handler called with " << _e.what());
// 		PyErr_SetString(PyExc_UserWarning, _e.what());
// 	});
}

#endif // XERUS_EXPERIMENTAL_PYTHON_WRAPPER
