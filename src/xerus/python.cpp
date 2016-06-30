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
#include <boost/python/call.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <numpy/ndarrayobject.h>
#pragma GCC diagnostic pop
#include "xerus.h"
#include "xerus/misc/internal.h"

using namespace boost::python;


void variable_argument_member_to_tuple_wrapper(const std::string &_name, const std::string &_tmpName = "new_fn") {
	boost::python::str pyCode;
	pyCode += "def patch_call_fn():\n";
	pyCode += std::string("  original = ")+_name+"\n";
	pyCode += std::string("  def ")+_tmpName+"( this, *args ):\n";
	pyCode += "    return original( this, args )\n";
	pyCode += std::string("  return ")+_tmpName+"\n";
	pyCode += _name + " = patch_call_fn()\n";
	boost::python::exec(pyCode, scope().attr("__dict__"));
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
		void* storage=(reinterpret_cast<converter::rvalue_from_python_storage<std::vector<containedType> >*>(data))->storage.bytes;
		new (storage) std::vector<containedType>();
		std::vector<containedType>* v=reinterpret_cast<std::vector<containedType>*>(storage);
		Py_ssize_t l=PySequence_Size(obj_ptr); if(l<0) abort();
		v->reserve(size_t(l)); 
		for (Py_ssize_t i=0; i<l; i++) { 
			v->push_back(extract<containedType>(PySequence_GetItem(obj_ptr,i))); 
		}
		data->convertible=storage;
	}
};


// TODO get_copy() wrapper in python correct manner
BOOST_PYTHON_MODULE(libxerus) {
	using namespace xerus;
	
	import_array(); // for numpy
	
	bool show_user_defined = true;
	bool show_signatures = false;
	docstring_options doc_options(show_user_defined, show_signatures);
	
	// xerus version
	scope().attr("VERSION_MAJOR") = VERSION_MAJOR;
	scope().attr("VERSION_MINOR") = VERSION_MINOR;
	scope().attr("VERSION_REVISION") = VERSION_REVISION;
	scope().attr("VERSION_COMMIT") = VERSION_COMMIT;
	
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
	exec(
		"def indices(n=1):\n"
		"  \"\"\"Create n distinct indices.\"\"\"\n"
		"  i = 0\n"
		"  while i<n:\n"
		"    yield Index()\n"
		"    i += 1\n"
	, scope().attr("__dict__"));
	VECTOR_TO_PY(Index, "IndexVector");
	
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
	VECTOR_TO_PY(IndexedTensor<Tensor>*, "IndexedTensorList");
	
	
	// ----------------------------------------------------------- factorizations
	class_<TensorFactorisation, boost::noncopyable>("TensorFactorisation", boost::python::no_init)
		.def("__rlshift__", +[](TensorFactorisation &_rhs, object &_lhs){
			std::vector<IndexedTensor<Tensor>*> tmp = extract<std::vector<IndexedTensor<Tensor>*>>(_lhs);
			_rhs(tmp);
		})
	;
	class_<SVD, bases<TensorFactorisation>, boost::noncopyable>("SVD_temporary", boost::python::no_init);
	def("SVD", +[](IndexedTensor<Tensor> &_rhs)->TensorFactorisation*{
		return new SVD(std::move(_rhs));
	}, return_value_policy<manage_new_object,  // result is treated as a new object
	   with_custodian_and_ward_postcall<0,1>>()); // but the argument will not be destroyed before the result is destroyed
	
	class_<QR, bases<TensorFactorisation>, boost::noncopyable>("QR_temporary", boost::python::no_init);
	def("QR", +[](IndexedTensor<Tensor> &_rhs)->TensorFactorisation*{
		return new QR(std::move(_rhs));
	}, return_value_policy<manage_new_object,  // result is treated as a new object
	   with_custodian_and_ward_postcall<0,1>>()); // but the argument will not be destroyed before the result is destroyed
	
	class_<RQ, bases<TensorFactorisation>, boost::noncopyable>("RQ_temporary", boost::python::no_init);
	def("RQ", +[](IndexedTensor<Tensor> &_rhs)->TensorFactorisation*{
		return new RQ(std::move(_rhs));
	}, return_value_policy<manage_new_object,  // result is treated as a new object
	   with_custodian_and_ward_postcall<0,1>>()); // but the argument will not be destroyed before the result is destroyed
	
	class_<QC, bases<TensorFactorisation>, boost::noncopyable>("QC_temporary", boost::python::no_init);
	def("QC", +[](IndexedTensor<Tensor> &_rhs)->TensorFactorisation*{
		return new QC(std::move(_rhs));
	}, return_value_policy<manage_new_object,  // result is treated as a new object
	   with_custodian_and_ward_postcall<0,1>>()); // but the argument will not be destroyed before the result is destroyed
	
	class_<CQ, bases<TensorFactorisation>, boost::noncopyable>("CQ_temporary", boost::python::no_init);
	def("CQ", +[](IndexedTensor<Tensor> &_rhs)->TensorFactorisation*{
		return new CQ(std::move(_rhs));
	}, return_value_policy<manage_new_object,  // result is treated as a new object
	   with_custodian_and_ward_postcall<0,1>>()); // but the argument will not be destroyed before the result is destroyed
	
	
	// ----------------------------------------------------------- Tensor
	{ scope Tensor_scope = 
		class_<Tensor>("Tensor",
			"a non-decomposed Tensor in either sparse or dense representation"
		)
			.def(init<const Tensor::DimensionTuple&>(args("dimensions"), "constructs a Tensor with the given dimensions"))
			.def(init<const TensorNetwork&>())
			.def(init<const Tensor &>())
			.def("from_function", +[](const Tensor::DimensionTuple& _dim, PyObject *_f){
				return Tensor(_dim, [&](std::vector<size_t> _p) {
					return boost::python::call<double>(_f, _p);
				});
			}).staticmethod("from_function") 
			.def("from_ndarray", +[](PyObject *_npObject){
				//TODO check for dangling pointers!
				#pragma GCC diagnostic push
				#pragma GCC diagnostic ignored "-Wuseless-cast"
				#pragma GCC diagnostic ignored "-Wold-style-cast"
				#pragma GCC diagnostic ignored "-Wcast-qual"
				#pragma GCC diagnostic ignored "-Wpedantic"
				PyArrayObject *npa = reinterpret_cast<PyArrayObject*>(PyArray_FROM_OTF(_npObject, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY));
				#pragma GCC diagnostic pop
				int deg = PyArray_NDIM(npa);
				std::vector<size_t> dims;
				dims.resize(size_t(deg));
				if (PyArray_ISCONTIGUOUS(npa)) {
					for (int i=0; i<deg; ++i) {
						dims[size_t(i)] = size_t(PyArray_DIMS(npa)[i]);
					}
					Tensor result(dims, Tensor::Representation::Dense, Tensor::Initialisation::None);
					misc::copy(result.get_unsanitized_dense_data(), static_cast<double*>(PyArray_DATA(npa)), result.size);
					Py_DECREF(npa);
					return object(result);
				} else if (PyArray_ISFORTRAN(npa)) {
					std::vector<size_t> shuffle(dims);
					for (int i=0; i<deg; ++i) {
						dims[size_t(deg-i-1)] = size_t(PyArray_DIMS(npa)[i]);
						shuffle[size_t(deg-i-1)] = size_t(i);
					}
					Tensor result(dims, Tensor::Representation::Dense, Tensor::Initialisation::None);
					// TODO reduce number of copies
					Tensor tmp(dims, Tensor::Representation::Dense, Tensor::Initialisation::None);
					misc::copy(tmp.get_unsanitized_dense_data(), static_cast<double*>(PyArray_DATA(npa)), result.size);
					reshuffle(result, tmp, shuffle);
					Py_DECREF(npa);
					return object(result);
				} else {
					LOG(error, "could not convert ndarray of neither c nor fortran striding");
					Py_DECREF(npa);
					return object();
				}
			}).staticmethod("from_ndarray")
			.def("to_ndarray", +[](Tensor &_this){
				std::vector<npy_intp> dimensions;
				for (size_t d : _this.dimensions) {
					dimensions.emplace_back(npy_intp(d));
				}
				
				#pragma GCC diagnostic push
				#pragma GCC diagnostic ignored "-Wuseless-cast"
				#pragma GCC diagnostic ignored "-Wold-style-cast"
				#pragma GCC diagnostic ignored "-Wcast-qual"
				#pragma GCC diagnostic ignored "-Wpedantic"
				PyObject *pyObj; 
				if (_this.is_dense()) {
					pyObj = PyArray_SimpleNewFromData(int(_this.degree()), &dimensions[0], NPY_DOUBLE, _this.get_dense_data());
				} else {
					Tensor cpy(_this); // NOTE leaves _this as a sparse tensor
					//TODO reduce the numer of copies in this
					pyObj = PyArray_SimpleNewFromData(int(_this.degree()), &dimensions[0], NPY_DOUBLE, cpy.get_dense_data());
				}
				PyObject * res = PyArray_Copy(reinterpret_cast<PyArrayObject*>(pyObj));  // copy due to lifetime issues (copy is owned by numpy instead of us)
				Py_DECREF(pyObj);
				return res;
				#pragma GCC diagnostic pop
			})
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
			.def("all_entries_valid", &Tensor::all_entries_valid)
			.def("reorder_cost", &Tensor::reorder_cost)
			.def("reinterpret_dimensions", &Tensor::reinterpret_dimensions)
			.def("resize_mode", &Tensor::resize_mode,
				(arg("mode"), arg("newDimension"), arg("cutPosition")=~0ul)
			)
			.def("fix_mode", &Tensor::fix_mode)
			.def("remove_slate", &Tensor::remove_slate)
			.def("perform_trace", &Tensor::perform_trace)
			.def("offset_add", &Tensor::offset_add)
			.def("use_dense_representation", &Tensor::use_dense_representation)
			.def("use_sparse_representation", &Tensor::use_sparse_representation,
				(arg("epsilon")=EPSILON)
			)
			.def("ensure_own_data", &Tensor::ensure_own_data)
			.def("ensure_own_data_no_copy", &Tensor::ensure_own_data_no_copy)
			.def("apply_factor", &Tensor::apply_factor)
			.def("ensure_own_data_and_apply_factor", &Tensor::ensure_own_data_and_apply_factor)
			.def("multiIndex_to_position", &Tensor::multiIndex_to_position).staticmethod("multiIndex_to_position")
			.def("position_to_multiIndex", &Tensor::position_to_multiIndex).staticmethod("position_to_multiIndex")
			.def("__call__", +[](Tensor &_this, const std::vector<Index> &_idx){
				return  new xerus::internal::IndexedTensor<Tensor>(std::move(_this(_idx)));
			}, return_value_policy<manage_new_object>())
			.def("__str__", &Tensor::to_string)
			.def(self * other<value_t>())
			.def(other<value_t>() * self)
			.def(self / other<value_t>())
			.def(self + self)
			.def(self - self)
			.def(self += self)
			.def(self -= self)
			.def("__getitem__", +[](Tensor &_this, size_t _i) {
				if (_i >= _this.size) {
					PyErr_SetString(PyExc_IndexError, "Index out of range");
					throw_error_already_set();
				}
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
	variable_argument_member_to_tuple_wrapper("Tensor.__call__", "TensorCallOperator");
	def("reshuffle", static_cast<Tensor(*)(const Tensor&, const std::vector<size_t>&)>(&reshuffle));
	def("contract", static_cast<Tensor(*)(const Tensor&, bool, const Tensor&, bool, size_t)>(&contract));
	def("pseudo_inverse", static_cast<Tensor(*)(const Tensor&, size_t)>(&pseudo_inverse));
	def("entrywise_product", static_cast<Tensor(*)(const Tensor&, const Tensor&)>(&entrywise_product));
	def("approx_entrywise_equal", static_cast<bool(*)(const Tensor&, const Tensor&, value_t)>(&approx_entrywise_equal));
	
	
	
	// ------------------------------------------------------------- TensorNetwork
	VECTOR_TO_PY(TensorNetwork::Link, "LinkVector");
	VECTOR_TO_PY(TensorNetwork::TensorNode, "TensorNetworkNodeVector");
	// TODO allow modification of the network
	{ scope TN_scope =
		class_<TensorNetwork>("TensorNetwork")
			.def(init<Tensor>())
			.def(init<const TensorNetwork &>())
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
			.def("__call__", +[](TensorNetwork &_this, const std::vector<Index> &_idx){
				return  new xerus::internal::IndexedTensor<TensorNetwork>(std::move(_this(_idx)));
			}, return_value_policy<manage_new_object>())
			.def(self * other<value_t>())
			.def(other<value_t>() * self)
			.def(self / other<value_t>())
			.def("__getitem__", +[](TensorNetwork &_this, size_t _i) {
				if (_i >= misc::product(_this.dimensions)) {
					PyErr_SetString(PyExc_IndexError, "Index out of range");
					throw_error_already_set();
				}
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
			.def("measure", &TensorNetwork::measure)
			.def("swap_external_links", &TensorNetwork::swap_external_links)
			.def("round_edge", &TensorNetwork::round_edge)
			.def("transfer_core", &TensorNetwork::transfer_core,
				(arg("from"), arg("to"), arg("allowRankReduction")=true)
			)
			.def("reduce_representation", &TensorNetwork::reduce_representation)
			//TODO find_common_edge (wrapper that returns python::tuple
			.def("sanitize", &TensorNetwork::sanitize)
			.def("fix_mode", &TensorNetwork::fix_mode)
			.def("remove_slate", &TensorNetwork::remove_slate)
			.def("resize_mode", &TensorNetwork::resize_mode,
				(arg("mode"), arg("newDimension"), arg("cutPosition")=~0ul)
			)
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
			.add_property("tensorObject", +[](TensorNetwork::TensorNode &_this)->object{
				if (_this.tensorObject) {
					return object(_this.tensorObject.get());
				} else {
					return object(); // aka none
				}
			})
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
	variable_argument_member_to_tuple_wrapper("TensorNetwork.__call__", "TensorNetworkCallOperator");
	
	// ------------------------------------------------------------- TTNetwork	
	class_<TTTensor, bases<TensorNetwork>>("TTTensor")
		.def(init<const Tensor&, optional<value_t, size_t>>())
		.def(init<const Tensor&, value_t, TensorNetwork::RankTuple>())
		.def(init<Tensor::DimensionTuple>())
		.def(init<size_t>())
		.def(init<const TTTensor &>())
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
// 		.def("kronecker", &TTTensor::kronecker).staticmethod("kronecker") //TODO
// 		.def("dirac", static_cast<TTTensor (*)(Tensor::DimensionTuple, const Tensor::MultiIndex&)>(&TTTensor::dirac)) //TODO
// 		.def("dirac", static_cast<TTTensor (*)(Tensor::DimensionTuple, const size_t)>(&TTTensor::dirac)).staticmethod("dirac") //TODO
		
		.def("reduce_to_maximal_ranks", &TTTensor::reduce_to_maximal_ranks).staticmethod("reduce_to_maximal_ranks")
// 		.def("degrees_of_freedom", static_cast<size_t (TTTensor::*)()>(&TTTensor::degrees_of_freedom))
		.def("degrees_of_freedom", static_cast<size_t (*)(const std::vector<size_t>&, const std::vector<size_t>&)>(&TTTensor::degrees_of_freedom)).staticmethod("degrees_of_freedom")
		.def("dyadic_product", static_cast<TTTensor (*)(const std::vector<TTTensor> &)>(&TTTensor::dyadic_product))
			.staticmethod("dyadic_product") //  TODO should not be in TTTensor scope
		.def("entrywise_square", &TTTensor::entrywise_square) // TODO should not be a member method?
		// TODO chop wrapper
		
		.def("round", static_cast<void (TTTensor::*)(const std::vector<size_t>&, double)>(&TTTensor::round),
			(arg("ranks"), arg("epsilon")=EPSILON)
		)
		.def("round", static_cast<void (TTTensor::*)(double)>(&TTTensor::round))
		.def("round", static_cast<void (TTTensor::*)(size_t)>(&TTTensor::round))
		
		.def("soft_threshold", static_cast<void (TTTensor::*)(const double, const bool)>(&TTTensor::soft_threshold),
			(arg("tau"), arg("preventZero")=false)
		)
		.def("soft_threshold", static_cast<void (TTTensor::*)(const std::vector<double>&, const bool)>(&TTTensor::soft_threshold),
			(arg("tau"), arg("preventZero")=false)
		)
		
		.def("move_core", &TTTensor::move_core,
			(arg("position"), arg("keepRank")=false)
		)
		
		.def("assume_core_position", &TTTensor::assume_core_position)
		.def("cannonicalize_left", &TTTensor::cannonicalize_left)
		.def("cannonicalize_right", &TTTensor::cannonicalize_right)
		.def(self + self)
		.def(self - self)
		.def(self * other<value_t>())
		.def(other<value_t>() * self)
		.def(self / other<value_t>())
		.def(self += self)
		.def(self -= self)
	;
	def("entrywise_product", static_cast<TTTensor (*)(const TTTensor&, const TTTensor&)>(&entrywise_product));
	def("find_largest_entry", static_cast<size_t (*)(const TTTensor&, value_t, value_t)>(&find_largest_entry));
	
	class_<TTOperator, bases<TensorNetwork>>("TTOperator")
		.def(init<const Tensor&, optional<value_t, size_t>>())
		.def(init<const Tensor&, value_t, TensorNetwork::RankTuple>())
		.def(init<Tensor::DimensionTuple>())
		.def(init<size_t>())
		.def(init<const TTOperator &>())
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
// 		.def("degrees_of_freedom", static_cast<size_t (TTOperator::*)()>(&TTOperator::degrees_of_freedom))
		.def("degrees_of_freedom", static_cast<size_t (*)(const std::vector<size_t>&, const std::vector<size_t>&)>(&TTOperator::degrees_of_freedom)).staticmethod("degrees_of_freedom")
		.def("dyadic_product", static_cast<TTOperator (*)(const std::vector<TTOperator> &)>(&TTOperator::dyadic_product))
			.staticmethod("dyadic_product") //  TODO should not be in TTOperator scope
		.def("entrywise_square", &TTOperator::entrywise_square) // TODO should not be a member method!
		// TODO chop wrapper
		
		.def("round", static_cast<void (TTOperator::*)(const std::vector<size_t>&, double)>(&TTOperator::round),
			(arg("ranks"), arg("epsilon")=EPSILON)
		)
		.def("round", static_cast<void (TTOperator::*)(double)>(&TTOperator::round))
		.def("round", static_cast<void (TTOperator::*)(size_t)>(&TTOperator::round))
		
		.def("soft_threshold", static_cast<void (TTOperator::*)(const double, const bool)>(&TTOperator::soft_threshold),
			(arg("tau"), arg("preventZero")=false)
		)
		.def("soft_threshold", static_cast<void (TTOperator::*)(const std::vector<double>&, const bool)>(&TTOperator::soft_threshold),
			(arg("tau"), arg("preventZero")=false)
		)
		
		.def("move_core", &TTOperator::move_core,
			(arg("position"), arg("keepRank")=false)
		)
		
		.def("assume_core_position", &TTOperator::assume_core_position)
		.def("cannonicalize_left", &TTOperator::cannonicalize_left)
		.def("cannonicalize_right", &TTOperator::cannonicalize_right)
		.def(self + self)
		.def(self - self)
		.def(self += self)
		.def(self -= self)
		.def(self * other<value_t>())
		.def(other<value_t>() * self)
		.def(self / other<value_t>())
		
		
		// for  TTOperator only:
		.def("identity", &TTOperator::identity<>).staticmethod("identity") 
		.def("transpose", &TTOperator::transpose<>)
	;
	def("entrywise_product", static_cast<TTOperator (*)(const TTOperator&, const TTOperator&)>(&entrywise_product));
	def("find_largest_entry", static_cast<size_t (*)(const TTOperator&, value_t, value_t)>(&find_largest_entry));
	
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
			})
			.add_property("errorFunction", 
						  +[](PerformanceData &_this){ return _this.errorFunction; }, 
						  +[](PerformanceData &_this, PyObject *_f){ 
							  // TODO increase ref count for _f? also decrease it on overwrite?!
							  _this.errorFunction = [_f](const TTTensor &_x)->double{
								  return call<double>(_f, _x);
							}; 
						})
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
			.def("add", +[](PerformanceData &_this, size_t _itr,  value_t _res, const TTTensor &_x, size_t _flags){
				_this.add(_itr, _res, _x, _flags);
			}, (arg("iterationCount"), arg("residual"), arg("x"), arg("flags")=0) )
			.def("add", +[](PerformanceData &_this, value_t _res, const TTTensor &_x, size_t _flags){
				_this.add(_res, _x, _flags);
			}, (arg("residual"), arg("x"), arg("flags")=0) )
			.def("__nonzero__", +[](PerformanceData &_this){ return bool(_this); })
			.def("dump_to_file", &PerformanceData::dump_to_file)
			.def("__iadd__", +[](PerformanceData &_this, const std::string &_s){
				_this << _s;
			})
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
	
	class_<TTRetractionI>("TTRetractionI", init<const TTRetractionI &>());
	scope().attr("ALSRetractionI") = object(TTRetractionI(&ALSRetractionI));
	scope().attr("SubmanifoldRetractionI") = object(TTRetractionI(&SubmanifoldRetractionI));
	scope().attr("HOSVDRetractionI") = object(TTRetractionI(&HOSVDRetractionI));
	
	class_<TTRetractionII>("TTRetractionII", init<const TTRetractionII &>());
	scope().attr("ALSRetractionII") = object(TTRetractionII(&ALSRetractionII));
	scope().attr("SubmanifoldRetractionII") = object(TTRetractionII(&SubmanifoldRetractionII));
	scope().attr("HOSVDRetractionII") = object(TTRetractionII(&HOSVDRetractionII));
	
	class_<TTVectorTransport>("TTVectorTransport", init<const TTVectorTransport &>());
	scope().attr("ProjectiveVectorTransport") = object(TTVectorTransport(&ProjectiveVectorTransport));
	
	{ scope als_scope = 
		class_<ALSVariant>("ALSVariant", init<uint, size_t, ALSVariant::LocalSolver, bool, optional<bool>>())
			.def(init<const ALSVariant&>())
			.def_readwrite("sites", &ALSVariant::sites)
			.def_readwrite("numHalfSweeps", &ALSVariant::numHalfSweeps)
			.def_readwrite("convergenceEpsilon", &ALSVariant::convergenceEpsilon)
			.def_readwrite("useResidualForEndCriterion", &ALSVariant::useResidualForEndCriterion)
			.def_readwrite("preserveCorePosition", &ALSVariant::preserveCorePosition)
			.def_readwrite("assumeSPD", &ALSVariant::assumeSPD)
			.add_property("localSolver", 
						  +[](ALSVariant &_this){ return _this.localSolver; },
						  +[](ALSVariant &_this, ALSVariant::LocalSolver _s){ _this.localSolver = _s; })
			
			.def("__call__", +[](GeometricCGVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
				_this(_A, _x, _b, _pd);
			}, (arg("A"), arg("x"), arg("b"), arg("perfData")=NoPerfData) )
			
			.def("__call__", +[](GeometricCGVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _eps, PerformanceData &_pd) {
				_this(_A, _x, _b, _eps, _pd);
			}, (arg("A"), arg("x"), arg("b"), arg("epsilon"), arg("perfData")=NoPerfData) )
			
			.def("__call__", +[](GeometricCGVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, PerformanceData &_pd) {
				_this(_A, _x, _b, _numHalfSweeps, _pd);
			}, (arg("A"), arg("x"), arg("b"), arg("numHalfSweeps"), arg("perfData")=NoPerfData) )
			
			.def("__call__", +[](GeometricCGVariant &_this, TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
				_this(_x, _b, _pd);
			}, (arg("x"), arg("b"), arg("perfData")=NoPerfData) )
			
			.def("__call__", +[](GeometricCGVariant &_this, TTTensor &_x, const TTTensor &_b, value_t _eps, PerformanceData &_pd) {
				_this(_x, _b, _eps, _pd);
			}, (arg("x"), arg("b"), arg("epsilon"), arg("perfData")=NoPerfData) )
			
			.def("__call__", +[](GeometricCGVariant &_this, TTTensor &_x, const TTTensor &_b, size_t _numHalfSweeps, PerformanceData &_pd) {
				_this(_x, _b, _numHalfSweeps, _pd);
			}, (arg("x"), arg("b"), arg("numHalfSweeps"), arg("perfData")=NoPerfData) )
		;
		class_<ALSVariant::LocalSolver>("LocalSolver", boost::python::no_init);
		als_scope.attr("lapack_solver") = object(ALSVariant::LocalSolver(&ALSVariant::lapack_solver));
		als_scope.attr("ASD_solver") = object(ALSVariant::LocalSolver(&ALSVariant::ASD_solver));
	}
	scope().attr("ALS") = object(ptr(&ALS));
	scope().attr("ALS_SPD") = object(ptr(&ALS_SPD));
	scope().attr("DMRG") = object(ptr(&DMRG));
	scope().attr("DMRG_SPD") = object(ptr(&DMRG_SPD));
	scope().attr("ASD") = object(ptr(&ASD));
	scope().attr("ASD_SPD") = object(ptr(&ASD_SPD));
	
	def("decomposition_als", &decomposition_als, (arg("x"), arg("b"), arg("epsilon")=EPSILON, arg("maxIterations")=1000));
	
	class_<GeometricCGVariant>("GeometricCGVariant", init<size_t, value_t, bool, TTRetractionI, TTVectorTransport>())
		.def(init<const GeometricCGVariant&>())
		.def_readwrite("numSteps", &GeometricCGVariant::numSteps)
		.def_readwrite("convergenceEpsilon", &GeometricCGVariant::convergenceEpsilon)
		.def_readwrite("assumeSymmetricPositiveDefiniteOperator", &GeometricCGVariant::assumeSymmetricPositiveDefiniteOperator)
		.add_property("retraction", 
					  +[](GeometricCGVariant &_this){ return _this.retraction; }, 
					  +[](GeometricCGVariant &_this, TTRetractionI _r){ _this.retraction = _r; })
		.add_property("vectorTransport", 
					  +[](GeometricCGVariant &_this){ return _this.vectorTransport; }, 
					  +[](GeometricCGVariant &_this, TTVectorTransport _transp){ _this.vectorTransport = _transp; })
		
		.def("__call__", +[](GeometricCGVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
			_this(_A, _x, _b, _pd);
		}, (arg("A"), arg("x"), arg("b"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](GeometricCGVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _eps, PerformanceData &_pd) {
			_this(_A, _x, _b, _eps, _pd);
		}, (arg("A"), arg("x"), arg("b"), arg("epsilon"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](GeometricCGVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numSteps, PerformanceData &_pd) {
			_this(_A, _x, _b, _numSteps, _pd);
		}, (arg("A"), arg("x"), arg("b"), arg("numSteps"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](GeometricCGVariant &_this, TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
			_this(_x, _b, _pd);
		}, (arg("x"), arg("b"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](GeometricCGVariant &_this, TTTensor &_x, const TTTensor &_b, value_t _eps, PerformanceData &_pd) {
			_this(_x, _b, _eps, _pd);
		}, (arg("x"), arg("b"), arg("epsilon"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](GeometricCGVariant &_this, TTTensor &_x, const TTTensor &_b, size_t _numSteps, PerformanceData &_pd) {
			_this(_x, _b, _numSteps, _pd);
		}, (arg("x"), arg("b"), arg("numSteps"), arg("perfData")=NoPerfData) )
	;
	scope().attr("GeometricCG") = object(ptr(&GeometricCG));
	
	class_<SteepestDescentVariant>("SteepestDescentVariant", init<size_t, value_t, bool, TTRetractionII>())
		.def(init<const SteepestDescentVariant&>())
		.def_readwrite("numSteps", &SteepestDescentVariant::numSteps)
		.def_readwrite("convergenceEpsilon", &SteepestDescentVariant::convergenceEpsilon)
		.def_readwrite("assumeSymmetricPositiveDefiniteOperator", &SteepestDescentVariant::assumeSymmetricPositiveDefiniteOperator)
		.add_property("retraction", 
					  +[](SteepestDescentVariant &_this){ return _this.retraction; }, 
					  +[](SteepestDescentVariant &_this, TTRetractionII _r){ _this.retraction = _r; })
		
		.def("__call__", +[](SteepestDescentVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
			_this(_A, _x, _b, _pd);
		}, (arg("A"), arg("x"), arg("b"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](SteepestDescentVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _eps, PerformanceData &_pd) {
			_this(_A, _x, _b, _eps, _pd);
		}, (arg("A"), arg("x"), arg("b"), arg("epsilon"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](SteepestDescentVariant &_this, const TTOperator &_A, TTTensor &_x, const TTTensor &_b, size_t _numSteps, PerformanceData &_pd) {
			_this(_A, _x, _b, _numSteps, _pd);
		}, (arg("A"), arg("x"), arg("b"), arg("numSteps"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](SteepestDescentVariant &_this, TTTensor &_x, const TTTensor &_b, PerformanceData &_pd) {
			_this(_x, _b, _pd);
		}, (arg("x"), arg("b"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](SteepestDescentVariant &_this, TTTensor &_x, const TTTensor &_b, value_t _eps, PerformanceData &_pd) {
			_this(_x, _b, _eps, _pd);
		}, (arg("x"), arg("b"), arg("epsilon"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](SteepestDescentVariant &_this, TTTensor &_x, const TTTensor &_b, size_t _numSteps, PerformanceData &_pd) {
			_this(_x, _b, _numSteps, _pd);
		}, (arg("x"), arg("b"), arg("numSteps"), arg("perfData")=NoPerfData) )
	;
	scope().attr("SteepestDescent") = object(ptr(&SteepestDescent));
	
	// ------------------------------------------------------------- measurements
	
	class_<SinglePointMeasurementSet>("SinglePointMeasurementSet")
		.def(init<const SinglePointMeasurementSet&>())
		.def("get_position", +[](SinglePointMeasurementSet &_this, size_t _i){
			return _this.positions[_i];
		})
		.def("set_position", +[](SinglePointMeasurementSet &_this, size_t _i, std::vector<size_t> _pos){
			_this.positions[_i] = _pos;
		})
		.def("get_measuredValue", +[](SinglePointMeasurementSet &_this, size_t _i){
			return _this.measuredValues[_i];
		})
		.def("set_measuredValue", +[](SinglePointMeasurementSet &_this, size_t _i, value_t _val){
			_this.measuredValues[_i] = _val;
		})
		.def("add", &SinglePointMeasurementSet::add)
		.def("size", &SinglePointMeasurementSet::size)
		.def("degree", &SinglePointMeasurementSet::degree)
		.def("test_solution", &SinglePointMeasurementSet::test_solution)
		
		.def("random", +[](const std::vector<size_t> &_dim, size_t _numMeas){
			static std::random_device rd;
			std::mt19937_64 rnd(rd());
			return SinglePointMeasurementSet::random(_dim, _numMeas, rnd);
		}).staticmethod("random")
	;
	def("sort", static_cast<void (*)(SinglePointMeasurementSet&, size_t)>(&xerus::sort), (arg("measurements"), arg("splitPosition")=~0ul) );
	def("IHT", &IHT, (arg("x"), arg("measurements"), arg("perfData")=NoPerfData) );
	
	
	VECTOR_TO_PY(Tensor, "TensorVector");
	
	class_<RankOneMeasurementSet>("RankOneMeasurementSet")
		.def(init<const RankOneMeasurementSet&>())
		.def("get_position", +[](RankOneMeasurementSet &_this, size_t _i){
			return _this.positions[_i];
		})
		.def("set_position", +[](RankOneMeasurementSet &_this, size_t _i, std::vector<Tensor> _pos){
			_this.positions[_i] = _pos;
		})
		.def("get_measuredValue", +[](RankOneMeasurementSet &_this, size_t _i){
			return _this.measuredValues[_i];
		})
		.def("set_measuredValue", +[](RankOneMeasurementSet &_this, size_t _i, value_t _val){
			_this.measuredValues[_i] = _val;
		})
		.def("add", &RankOneMeasurementSet::add)
		.def("size", &RankOneMeasurementSet::size)
		.def("degree", &RankOneMeasurementSet::degree)
		.def("test_solution", &RankOneMeasurementSet::test_solution)
	;
	def("sort", static_cast<void (*)(RankOneMeasurementSet&, size_t)>(&xerus::sort), (arg("measurements"), arg("splitPosition")=~0ul) );
	
	
	class_<ADFVariant>("ADFVariant", init<size_t, double, double>())
		.def(init<ADFVariant>())
		.def_readwrite("maxIterations", &ADFVariant::maxIterations)
		.def_readwrite("targetResidualNorm", &ADFVariant::targetResidualNorm)
		.def_readwrite("minimalResidualNormDecrease", &ADFVariant::minimalResidualNormDecrease)
		
		.def("__call__", +[](ADFVariant &_this, TTTensor& _x, const SinglePointMeasurementSet& _meas, PerformanceData& _pd){
			return _this(_x, _meas, _pd);
		}, (arg("x"), arg("measurements"), arg("perfData")=NoPerfData) )
		.def("__call__", +[](ADFVariant &_this, TTTensor& _x, const SinglePointMeasurementSet& _meas, const std::vector<size_t>& _maxRanks, PerformanceData& _pd){
			return _this(_x, _meas, _maxRanks, _pd);
		}, (arg("x"), arg("measurements"), arg("maxRanks"), arg("perfData")=NoPerfData) )
		
		.def("__call__", +[](ADFVariant &_this, TTTensor& _x, const RankOneMeasurementSet& _meas, PerformanceData& _pd){
			return _this(_x, _meas, _pd);
		}, (arg("x"), arg("measurements"), arg("perfData")=NoPerfData) )
		.def("__call__", +[](ADFVariant &_this, TTTensor& _x, const RankOneMeasurementSet& _meas, const std::vector<size_t>& _maxRanks, PerformanceData& _pd){
			return _this(_x, _meas, _maxRanks, _pd);
		}, (arg("x"), arg("measurements"), arg("maxRanks"), arg("perfData")=NoPerfData) )
	;
	scope().attr("ADF") = object(ptr(&ADF));
	
	
	// ------------------------------------------------------------- misc
	def("frob_norm", +[](const Tensor& _x){ return _x.frob_norm(); });
	def("frob_norm", +[](const TensorNetwork& _x){ return _x.frob_norm(); });
	def("frob_norm", static_cast<value_t (*)(const IndexedTensorReadOnly<Tensor>&)>(&frob_norm));
	def("frob_norm", static_cast<value_t (*)(const IndexedTensorReadOnly<TensorNetwork>&)>(&frob_norm));
	
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
	def("approx_equal", +[](const TensorNetwork& _l, const Tensor& _r) {
		return approx_equal(_l, _r);
	});
	def("approx_equal", +[](const TensorNetwork& _l, const TensorNetwork& _r) {
		return approx_equal(_l, _r);
	});
	
	def("log", +[](std::string _msg){
		LOG_SHORT(info, _msg);
	});
	
	enum_<misc::FileFormat>("FileFormat")
		.value("BINARY", misc::FileFormat::BINARY)
		.value("TSV", misc::FileFormat::TSV)
	;
	
	def("save_to_file", +[](const Tensor &_obj, const std::string &_filename, misc::FileFormat _format){
		misc::save_to_file(_obj, _filename, _format);
	}, (arg("object"), arg("filename"), arg("format")=misc::FileFormat::BINARY) );
	
	def("save_to_file", +[](const TensorNetwork &_obj, const std::string &_filename, misc::FileFormat _format){
		misc::save_to_file(_obj, _filename, _format);
	}, (arg("object"), arg("filename"), arg("format")=misc::FileFormat::BINARY) );
	
	def("save_to_file", +[](const TTTensor &_obj, const std::string &_filename, misc::FileFormat _format){
		misc::save_to_file(_obj, _filename, _format);
	}, (arg("object"), arg("filename"), arg("format")=misc::FileFormat::BINARY) );
	
	def("save_to_file", +[](const TTOperator &_obj, const std::string &_filename, misc::FileFormat _format){
		misc::save_to_file(_obj, _filename, _format);
	}, (arg("object"), arg("filename"), arg("format")=misc::FileFormat::BINARY) );
	
	def("load_from_file", +[](std::string _filename){
		// determine type stored in the file
		std::ifstream in(_filename);
		if (!in) {
			return object();
		}
		std::string classname;
		in >> classname; // "Xerus"
		in >> classname;
		in.close();
		if (classname == "xerus::Tensor") {
			return object(misc::load_from_file<Tensor>(_filename));
		}
		if (classname == "xerus::TensorNetwork") {
			return object(misc::load_from_file<TensorNetwork>(_filename));
		}
		if (classname == "xerus::TTNetwork<false>") {
			return object(misc::load_from_file<TTTensor>(_filename));
		}
		if (classname == "xerus::TTNetwork<true>") {
			return object(misc::load_from_file<TTOperator>(_filename));
		}
		LOG(warning, "unknown class type '" << classname << "' in file '" << _filename << "'");
		return object();
	});
	
	// identity returns the cpp name to a python object
// 	def("identity", identity_);
	
	// the following is probably not necessary because generic_error inherits from std::exception
// 	register_exception_translator<misc::generic_error>([](const misc::generic_error &_e){
// 		LOG(pydebug, "custom exception handler called with " << _e.what());
// 		PyErr_SetString(PyExc_UserWarning, _e.what());
// 	});
}

#endif // XERUS_EXPERIMENTAL_PYTHON_WRAPPER
