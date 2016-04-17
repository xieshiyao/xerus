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
	
	
	// ---------------------------------------------- indexedTensor<TN>
	
	class_<internal::IndexedTensorReadOnly<TensorNetwork>, boost::noncopyable>("IndexedTensorNetworkReadOnly", no_init)
		.def("__add__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, internal::IndexedTensorReadOnly<TensorNetwork> &_r) {
				LOG(pydebug, "add TN ro + ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) + std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__sub__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, internal::IndexedTensorReadOnly<TensorNetwork> &_r){
				LOG(pydebug, "sub TN ro - ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, internal::IndexedTensorReadOnly<TensorNetwork> &_r){
				LOG(pydebug, "mul TN ro * ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, internal::IndexedTensorMoveable<TensorNetwork> &_r){
				LOG(pydebug, "mul TN ro * mv");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, value_t _r){
				LOG(pydebug, "mul TN ro * scalar");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * _r));
			}, return_value_policy<manage_new_object>())
		.def("__rmul__", 
			+[](value_t _r, internal::IndexedTensorReadOnly<TensorNetwork> &_l){
				LOG(pydebug, "mul TN scalar * ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * _r));
			}, return_value_policy<manage_new_object>())
		.def("__div__", 
			+[](internal::IndexedTensorReadOnly<TensorNetwork> &_l, value_t _r){
				LOG(pydebug, "div TN ro / scalar");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) / _r));
			}, return_value_policy<manage_new_object>())
// 		.def("frob_norm", &internal::IndexedTensorReadOnly<TensorNetwork>::frob_norm)
	;
	class_<internal::IndexedTensorWritable<TensorNetwork>, boost::noncopyable, bases<internal::IndexedTensorReadOnly<TensorNetwork>>>("IndexedTensorNetworkWriteable", no_init)
	;
	class_<internal::IndexedTensorMoveable<TensorNetwork>, boost::noncopyable, bases<internal::IndexedTensorWritable<TensorNetwork>>>("IndexedTensorNetworkMoveable", no_init)
		.def("__mul__", 
			+[](internal::IndexedTensorMoveable<TensorNetwork> &_l, internal::IndexedTensorReadOnly<TensorNetwork> &_r){
				LOG(pydebug, "mul TN mv * ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * std::move(_r)));
			}, return_value_policy<manage_new_object>())
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
		.def("__add__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorReadOnly<Tensor> &_r) {
				LOG(pydebug, "add ro + ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) + std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__add__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorMoveable<Tensor> &_r) {
				LOG(pydebug, "add ro + mv");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) + std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__sub__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorReadOnly<Tensor> &_r){
				LOG(pydebug, "sub ro - ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__sub__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorMoveable<Tensor> &_r){
				LOG(pydebug, "sub ro - mv");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorReadOnly<Tensor> &_r){
				LOG(pydebug, "mul ro * ro");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorMoveable<Tensor> &_r){
				LOG(pydebug, "mul ro * mv");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, value_t _r){
				LOG(pydebug, "mul ro * scalar");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) * _r));
			}, return_value_policy<manage_new_object>())
		.def("__rmul__", 
			+[](value_t _r, internal::IndexedTensorReadOnly<Tensor> &_l){
				LOG(pydebug, "mul scalar * ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) * _r));
			}, return_value_policy<manage_new_object>())
		.def("__div__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, value_t _r){
				LOG(pydebug, "div ro / scalar");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) / _r));
			}, return_value_policy<manage_new_object>())
		.def("__div__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorReadOnly<Tensor> &_r){
				LOG(pydebug, "div ro / ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) / std::move(_r)));
			}, return_value_policy<manage_new_object>())
// 		.def("frob_norm", &internal::IndexedTensorReadOnly<Tensor>::frob_norm)
		
	;
	class_<internal::IndexedTensorWritable<Tensor>, boost::noncopyable, bases<internal::IndexedTensorReadOnly<Tensor>>>("IndexedTensorWriteable", no_init)
	;
	class_<internal::IndexedTensorMoveable<Tensor>, boost::noncopyable, bases<internal::IndexedTensorWritable<Tensor>>>("IndexedTensorMoveable", no_init)
		.def("__add__", 
			+[](internal::IndexedTensorMoveable<Tensor> &_l, internal::IndexedTensorReadOnly<Tensor> &_r) {
				LOG(pydebug, "add mv + ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) + std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__sub__", 
			+[](internal::IndexedTensorMoveable<Tensor> &_l, internal::IndexedTensorReadOnly<Tensor> &_r){
				LOG(pydebug, "sub mv - ro");
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__mul__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorMoveable<Tensor> &_r){
				LOG(pydebug, "mul ro * mv");
				return new internal::IndexedTensorMoveable<TensorNetwork>(std::move(std::move(_l) * std::move(_r)));
			}, return_value_policy<manage_new_object>())
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
	
	implicitly_convertible<internal::IndexedTensorReadOnly<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	implicitly_convertible<internal::IndexedTensorWritable<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	implicitly_convertible<internal::IndexedTensorMoveable<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	implicitly_convertible<internal::IndexedTensor<Tensor>, internal::IndexedTensorMoveable<TensorNetwork>>();
	
	// ----------------------------------------------------------- decompositions
	
	
	// ----------------------------------------------------------- Tensor
	class_<Tensor>("Tensor")
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
// 		.def("dirac", static_cast<Tensor (*)(Tensor::DimensionTuple, Tensor::MultiIndex)>(&Tensor::dirac)).staticmethod("dirac")
// 		.def("dirac", static_cast<Tensor (*)(Tensor::DimensionTuple, size_t)>(&Tensor::dirac)).staticmethod("dirac")
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
		.def("frob_norm", &Tensor::frob_norm)
	;
	
	
	// ------------------------------------------------------------- TTNetwork
	
	
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
