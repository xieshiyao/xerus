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
/*
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

static xerus::internal::IndexedTensor<xerus::Tensor>* wrapper_test(xerus::Tensor &_this) {
	xerus::internal::IndexedTensor<xerus::Tensor> *result = new xerus::internal::IndexedTensor<xerus::Tensor>(std::move(_this()));
	return result;
}

static std::string logg(std::string _s) {
	LOG(testing, "python string?! = '" << _s << "'");
	return _s + " horst";
}

BOOST_PYTHON_MODULE(libxerus) {
	class_<xerus::Index>("Index");
	class_<xerus::internal::IndexedTensor<xerus::Tensor>, boost::noncopyable>("IndexedTensor", no_init);
	class_<xerus::Tensor>("Tensor")
		.def("random", &random_tensor).staticmethod("random")
		.def("__call__", &wrapper_test, return_value_policy<manage_new_object>());
	def("log", &logg);
}*/

