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

template<typename T>
static std::vector<T> to_std_vector(const object& iterable) {
	return std::vector<T>(stl_input_iterator<T>(iterable),
							stl_input_iterator<T>());
}

template<class... args>
static xerus::internal::IndexedTensor<xerus::Tensor>* call_operator_tensor(xerus::Tensor &_this, args... _indices) {
	return new xerus::internal::IndexedTensor<xerus::Tensor>(std::move(_this(_indices...)));
}


BOOST_PYTHON_MODULE(libxerus) {
	using namespace xerus;
	
	class_<std::vector<size_t>, boost::noncopyable>("IntegerVector", no_init); // just to give it a human readable name in python
	
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
		ADD_MOVE_AND_RESULT_PTR("__sub__", -, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorReadOnly<TensorNetwork>, IndexedTensorMoveable<TensorNetwork>)
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
		.def_readonly("dimensions", &Tensor::dimensions)
		.def("degree", &Tensor::degree)
		.def_readonly("factor", &Tensor::factor)
		.def_readonly("size", &Tensor::size)
		.def("frob_norm", &Tensor::frob_norm)
		.def("random", 
			+[](boost::python::list _dim) {
				static std::random_device rd;
				std::mt19937_64 rnd(rd());
				std::normal_distribution<double> dist(0.0, 1.0);
				return xerus::Tensor::random(to_std_vector<size_t>(_dim), rnd, dist);
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
		.def("__call__", &call_operator_tensor<>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__str__", &Tensor::to_string)
		.def(self * other<value_t>())
		.def(other<value_t>() * self)
		.def(self / other<value_t>())
		.def(self + self)
		.def(self - self)
	;
	
	
	// ------------------------------------------------------------- TTNetwork
	class_<TensorNetwork>("TensorNetwork")
		.def(init<Tensor>())
		.def_readonly("dimensions", &TensorNetwork::dimensions)
// 		.def("__call__", &call_operator_tensor<>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index, Index>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index, Index, Index>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index, Index, Index, Index>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
// 		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def(self * other<value_t>())
		.def(other<value_t>() * self)
		.def(self / other<value_t>())
		.def("frob_norm", &TensorNetwork::frob_norm)
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
}

#endif // XERUS_EXPERIMENTAL_PYTHON_WRAPPER
