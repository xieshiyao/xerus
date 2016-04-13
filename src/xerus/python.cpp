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

static xerus::Tensor random_tensor(boost::python::list _dim) {
	static std::random_device rd;
	std::mt19937_64 rnd(rd());
	std::normal_distribution<double> dist(0.0, 1.0);
	return xerus::Tensor::random(to_std_vector<size_t>(_dim), rnd, dist);
}

template<class... args>
static xerus::internal::IndexedTensor<xerus::Tensor>* call_operator_tensor(xerus::Tensor &_this, args... _indices) {
	return new xerus::internal::IndexedTensor<xerus::Tensor>(std::move(_this(_indices...)));
}


BOOST_PYTHON_MODULE(libxerus) {
	using namespace xerus;
	class_<Index>("Index")
		.def(init<int64_t>())
		.def("__pow__", &Index::operator^)
		.def("__xor__", &Index::operator^)
		.def("__div__", &Index::operator/)
		.def("__and__", &Index::operator&)
		.def("__str__", static_cast<std::string (*)(const Index &)>(&misc::to_string<Index>))
	;
	implicitly_convertible<int64_t, Index>();
	
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
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
		.def("__sub__", 
			+[](internal::IndexedTensorReadOnly<Tensor> &_l, internal::IndexedTensorMoveable<Tensor> &_r){
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
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
				return new internal::IndexedTensorMoveable<Tensor>(std::move(std::move(_l) - std::move(_r)));
			}, return_value_policy<manage_new_object>())
	;
	class_<internal::IndexedTensor<Tensor>, boost::noncopyable, bases<internal::IndexedTensorWritable<Tensor>>>("IndexedTensor", no_init)
		.def("__lshift__", 
			+[](xerus::internal::IndexedTensor<xerus::Tensor> &_lhs, xerus::internal::IndexedTensorReadOnly<xerus::Tensor> &_rhs) {
				std::move(_lhs) = std::move(_rhs);
			})
	;
	
	class_<Tensor>("Tensor")
		.def("random", &random_tensor).staticmethod("random")
		.def("__call__", &call_operator_tensor<>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__call__", &call_operator_tensor<Index, Index, Index, Index, Index, Index, Index>, return_value_policy<manage_new_object>())
		.def("__str__", &Tensor::to_string)
	;
	
	def("log", +[](std::string _msg){
		LOG_SHORT(info, _msg);
	});
}

#endif // XERUS_EXPERIMENTAL_PYTHON_WRAPPER
