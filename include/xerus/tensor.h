// // Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2015 Benjamin Huber and Sebastian Wolf. 
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
 * @brief Header file for the Tensor class.
 */

#pragma once

#include "basic.h"
#include "misc/sfinae.h"
#include <string>
#include <limits>
#include <memory>
#include "indexedTensorReadOnly.h"
#include "indexedTensor.h"
 
namespace xerus {
	class FullTensor;
	class SparseTensor;
	
	/// @brief Class that handles simple (non-decomposed) tensors in a dense or sparse representation.
	class Tensor {
	public:
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Auxiliary types- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		enum class Representation : bool { Dense, Sparse};
		
		enum class Initialisation : bool { Nothing, Zero};
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/// @brief Vector containing the individual dimensions of the tensor.
		std::vector<size_t> dimensions;
		
		/// @brief Size of the Tensor -- always equal to the product of the dimensions.
		size_t size = 1;
		
		/// @brief The current representation of the Tensor (i.e Dense or Sparse)
		Representation representation = Representation::Dense;
		
		/// @brief Single value representing a constant scaling factor.
		value_t factor = 1.0;
		
		/** 
		 * @brief Shared pointer to the dense data array, if representation is dense. 
		 * @details The data is stored such that indices increase from right to left (row-major order). 
		 * If the tensor is modified and not sole owner a deep copy is performed.
		 */
		std::shared_ptr<value_t> denseData;
		
		/** 
		 * @brief Shared pointer to the a map containing the entries of the SparseTensor. 
		 * @details The entries are stored in a map which uses the position of each entry assuming row-major ordering as key value.
		 * If the tensor is modified and not sole owner a deep copy is performed.
		 */
		std::shared_ptr<std::map<size_t, value_t>> sparseData;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/// @brief Constructs an order zero Tensor with the given inital representation
		explicit Tensor(const Representation _representation = Representation::Dense);
		
		/// @brief Tensors are copy constructable.
		implicit Tensor( const Tensor& _other ) = default;
		
		/// @brief Tensors are move constructable.
		implicit Tensor( Tensor&& _other ) = default;

		/** 
		 * @brief: Creates a new tensor with the given dimensions.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _representation (optional) the initial representation of the tensor.
		 * @param _init (optional) inital data treatment, i.e. whether the tensor is to be zero Initialized.
		 */
		explicit Tensor(const std::vector<size_t>& _dimensions, const Representation _representation = Representation::Dense, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief: Creates a new tensor with the given dimensions.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _representation (optional) the initial representation of the tensor.
		 * @param _init (optional) inital data treatment, i.e. whether the tensor is to be zero Initialized.
		 */
		explicit Tensor(      std::vector<size_t>&& _dimensions, const Representation _representation = Representation::Dense, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief: Creates a new tensor with the given dimensions.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _representation (optional) the initial representation of the tensor.
		 * @param _init (optional) inital data treatment, i.e. whether the tensor is to be zero Initialized.
		 */
		explicit Tensor(std::initializer_list<size_t>&& _dimensions, const Representation _representation = Representation::Dense, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief: Creates a new (dense) tensor with the given dimensions, using a provided data.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _data inital dense data in row-major order.
		 */
		template<ADD_MOVE(Vec, std::vector<size_t>), ADD_MOVE(SPtr, std::shared_ptr<value_t>)>
		explicit Tensor(Vec&& _dimensions, SPtr&& _data)
		: dimensions(std::forward<Vec>(_dimensions)), size(misc::product(dimensions)), representation(Representation::Dense), denseData(std::forward<SPtr>(_data)) { }
		
		
		
		/// @brief Returns a pointer containing a copy of the tensor with same type (i.e. FullTensor or SparseTensor).
		virtual Tensor* get_copy() const = 0;
		
		/// @brief Returns a pointer containing a moved copy of the object with same type (i.e. FullTensor or SparseTensor).
		virtual Tensor* get_moved_copy() = 0;
		
		/// @brief Returns a pointer to a newly constructed order zero tensor of same type (i.e. FullTensor or SparseTensor) with entry equals zero.
		virtual Tensor* construct_new() const = 0;
		
		/** 
		 * @brief: Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with all entries set to zero and global factor one.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions) const = 0;
		
		/** 
		 * @brief: Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with all entries set to zero and global factor one.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		virtual Tensor* construct_new(	  std::vector<size_t>&& _dimensions) const = 0;
		
		/** 
		 * @brief: Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with all entries set to zero and global factor one.
		 * @details The second parameter is a DONT_SET_ZERO helper object that is only used to provide the function overload.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		virtual Tensor* construct_new(const std::vector<size_t>&  _dimensions, _unused_ DONT_SET_ZERO) const = 0;
		
		/** 
		 * @brief: Returns a pointer to a newly constructed tensor of same type (i.e. FullTensor or SparseTensor) with all entries set to zero and global factor one.
		 * @details The second parameter is a DONT_SET_ZERO helper object that is only used to provide the function overload.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		virtual Tensor* construct_new(	  std::vector<size_t>&& _dimensions, _unused_ DONT_SET_ZERO) const = 0;
		
		/** 
		 * @brief: Returns a FullTensor with all entries equal to one.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		static FullTensor ones(const std::vector<size_t>& _dimensions);
		
		/** 
		 * @brief: Returns a SparseTensor representation of the identity operator with the given dimensions.
		 * @details That is combining the first half of the dimensions and the second half of the dimensions results in an identity matrix.
		 * @param _dimensions the dimensions of the new tensor. It is required that _dimensions[i] = _dimensions[d/2+i], otherwise this cannot be the identity operator.
		 */
		static SparseTensor identity(const std::vector<size_t>& _dimensions);
		
		/** 
		 * @brief: Returns a SparseTensor representation of the kronecker delta.
		 * @details That is each entry is one if all indices are equal and zero otherwise. Note iff d=2 this coincides with identity.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		static SparseTensor kronecker(const std::vector<size_t>& _dimensions);
		
		/** 
		 * @brief: Returns a SparseTensor with a single entry equals oen and all other zero.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _position The position of the one
		 */
		static SparseTensor dirac(const std::vector<size_t>& _dimensions, const std::vector<size_t>& _position);
		
		/** 
		 * @brief: Returns a SparseTensor with a single entry equals oen and all other zero.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _position The position of the one
		 */
		static SparseTensor dirac(const std::vector<size_t>& _dimensions, const size_t _position);
		
		
		
		
		/// @brief Destructor
		virtual ~Tensor();
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/// @brief Ensures that this tensor is the sole owner of its data. If needed new space is allocated and all entries are copied.
		virtual void ensure_own_data() = 0;
		
		/// @brief Ensures that this tensor is the sole owner of its data space. If needed new space is allocated with entries left undefined.
		virtual void ensure_own_data_no_copy() = 0;
		
		/// @brief Checks whether there is a non-trivial scaling factor and applies it if nessecary.
		virtual void apply_factor() = 0;
		
		/// @brief Checks whether there is a non-trivial factor and applies it. Even if no factor is applied ensure_own_data() is called.
		virtual void ensure_own_data_and_apply_factor() = 0;
		
		/** 
		 * @brief Resets the tensor to the given dimensions, factor equals one and undefined entries.
		 * @details The second parameter is a DONT_SET_ZERO helper object that is only used to provide the function overload.
		 * @param _newDim the new dimensions the tensor shall have
		 */
		virtual void reset(const std::vector<size_t>&  _newDim, _unused_ DONT_SET_ZERO) = 0;
		
		/** 
		 * @brief Resets the tensor to the given dimensions, factor equals one and undefined entries.
		 * @details The second parameter is a DONT_SET_ZERO helper object that is only used to provide the function overload.
		 * @param _newDim the new dimensions the tensor shall have
		 */
		virtual void reset(	  std::vector<size_t>&& _newDim, _unused_ DONT_SET_ZERO) = 0;
		
		/** 
		 * @brief Resets the tensor to the given dimensions, factor equals one and all entries equals zero.
		 * @param _newDim the new dimensions the tensor shall have
		 */
		virtual void reset(const std::vector<size_t>&  _newDim) = 0;
		
		/** 
		 * @brief Resets the tensor to the given dimensions, factor equals one and all entries equals zero.
		 * @param _newDim the new dimensions the tensor shall have
		 */
		virtual void reset(	  std::vector<size_t>&& _newDim) = 0;
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
		/** 
		 * @brief Performs the entrywise multiplication with a constant @a _factor.
		 * @details Internally this only results in a change in the global factor.
		 * @param _factor the factor,
		 * @return a reference to this Tensor.
		 */
		Tensor& operator*=(const value_t _factor);
		
		
		/** 
		 * @brief Performs the entrywise divison by a constant @a _divisor.
		 * @details Internally this only results in a change in the global factor.
		 * @param _divisor the divisor,
		 * @return a reference to this Tensor.
		 */ 
		Tensor& operator/=(const value_t _divisor);
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/** 
		 * @brief Read/Write access a single entry.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return a reference to the selected entry.
		 */
		virtual value_t& operator[](const size_t _position) = 0;
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return the selected entry.
		 */
		virtual value_t operator[](const size_t _position) const = 0;
		
		/** 
		 * @brief Read/Write access a single entry.
		 * @param _position the position of the desired entry.
		 * @return a reference to the selected entry.
		 */
		virtual value_t& operator[](const std::vector<size_t>& _positions) = 0;
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry.
		 * @return the selected entry.
		 */
		virtual value_t operator[](const std::vector<size_t>& _positions) const = 0;
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/** 
		 * @brief Indexes the Tensor for read/write use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor.
		 */
		template<typename... args>
		IndexedTensor<Tensor> operator()(args... _args) {
			return IndexedTensor<Tensor>(this, std::vector<Index>({_args...}), false);
		}
		
		
		/** 
		 * @brief Indexes the Tensor for read only use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor.
		 */
		template<typename... args>
		IndexedTensorReadOnly<Tensor> operator()(args... _args) const {
			return IndexedTensorReadOnly<Tensor>(this, std::vector<Index>({_args...}));
		}
		
		
		/** 
		 * @brief Indexes the tensor for read/write use.
		 * @param _indices several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor.
		 */
		IndexedTensor<Tensor> operator()(const std::vector<Index>&  _indices);
		
		
		/** 
		 * @brief Indexes the tensor for read/write use.
		 * @param _indices several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor.
		 */
		IndexedTensor<Tensor> operator()(	  std::vector<Index>&& _indices);

		
		/** 
		 * @brief Indexes the tensor for read only use.
		 * @param _indices several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor.
		 */
		IndexedTensorReadOnly<Tensor> operator()(const std::vector<Index>&  _indices) const;
		
		
		/** 
		 * @brief Indexes the tensor for read only use.
		 * @param _indices Several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor.
		 */
		IndexedTensorReadOnly<Tensor> operator()(	  std::vector<Index>&& _indices) const;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/** 
		 * @brief Checks whether the object is of type SparseTensor.
		 * @return true is sparse, false if not.
		 */
		virtual bool is_sparse() const = 0;
		
		
		/** 
		 * @brief Returns the degree of the tensor.
		 * @details The degree is always equals to dimensions.size()
		 * @return the degree of the tensor
		 */
		size_t degree() const;
		
		/** 
		 * @brief Checks whether the tensor has a non-trivial global scaling factor.
		 * @return true if there is a non-trivial factor, false if not.
		 */
		bool has_factor() const;
		
		/** 
		 * @brief Determines the number of non-zero entries.
		 * @param _eps (optional) epsilon detrmining the maximal value, that is still assumed to be zero.
		 * @return the number of non-zero entries found.
		 */
		virtual size_t count_non_zero_entries(const value_t _eps = std::numeric_limits<value_t>::epsilon()) const = 0;
		
		/** 
		 * @brief Checks the tensor for illegal entries, e.g. nan, inf,...
		 * @return TRUE there are no invalid entries, FALSE otherwise.
		 */
		virtual bool all_entries_valid() const = 0;
		
		/** 
		 * @brief Approximates the cost to reorder the tensor.
		 * @return the approximated costs.
		 */
		size_t reorder_costs() const;

		
		/** 
		 * @brief Calculates the frobenious norm of the tensor.
		 * @return the frobenious norm.
		 */
		virtual value_t frob_norm() const = 0;
		
		/** 
		 * @brief Reinterprets the dimensions of the tensor.
		 * @details For this simple reinterpretation it is nessecary that the size implied by the new dimensions is the same as to old size 
		 * (a vector with 16 entries cannot be interpreted as a 10x10 matrix, but it can be interpreted as a 4x4 matrix). If a real change in dimensions is 
		 * required use change_dimensions() instead.
		 * @param _newDimensions the dimensions the tensor shall be interpreted to have. 
		 */
		void reinterpret_dimensions(const std::vector<size_t>& _newDimensions);
		

		/** 
		 * @brief Reinterprets the dimensions of the tensor.
		 * @details For this simple reinterpretation it is nessecary that the size implied by the new dimensions is the same as to old size 
		 * (a vector with 16 entries cannot be interpreted as a 10x10 matrix, but it can be interpreted as a 4x4 matrix). If a real change in dimensions is 
		 * required use change_dimensions() instead.
		 * @param _newDimensions the dimensions the tensor shall be interpreted to have. 
		 */
		void reinterpret_dimensions(	  std::vector<size_t>&& _newDimensions);
		
		
		/** 
		 * @brief Reinterprets the dimensions of the tensor.
		 * @details For this simple reinterpretation it is nessecary that the size implied by the new dimensions is the same as to old size 
		 * (a vector with 16 entries cannot be interpreted as a 10x10 matrix, but it can be interpreted as a 4x4 matrix). If a real change in dimensions is 
		 * required use change_dimensions() instead.
		 * @param _newDimensions the dimensions the tensor shall be interpreted to have. 
		 */
		void reinterpret_dimensions( std::initializer_list<size_t>&& _newDimensions);
		
		/** 
		 * @brief Fixes a specific slate in one of the dimensions, effectively reducing the order by one.
		 * @param _dimension the dimension in which the slate shall be fixed, e.g. 0 to fix the first dimensions.
		 * @param _slatePosition the position in the corresponding dimensions that shall be used.
		 */
		virtual void fix_slate(const size_t _dimension, const size_t _slatePosition) = 0;
		
		/** 
		 * @brief Creates a string representation of the Tensor.
		 * @param _newDimensions the string representation. 
		 */
		virtual std::string to_string() const = 0;
		
		
		/** 
		 * @brief Compares the Tensor entriewise to the given data.
		 * @param _values data to compare to, must be of the same size as the Tensor has entries (i.e. size).
		 * @param _eps (optional) epsilon used to determine whether an entry is equal to a data point. 
		 * @return TRUE if all entries deviate by less than _eps from the correspondign data points, FALSE otherwise.
		 */
		virtual bool compare_to_data(const std::vector<value_t>& _values, const double _eps = std::numeric_limits<value_t>::epsilon()) const = 0;
		
		
		/** 
		 * @brief Compares the Tensor entriewise to the given data.
		 * @param _values data to compare to, must be of the same size as the Tensor has entries (i.e. size).
		 * @param _eps (optional) epsilon used to determine whether an entry is equal to a data point. 
		 * @return TRUE if all entries deviate by less than _eps from the correspondign data points, FALSE otherwise.
		 */
		virtual bool compare_to_data(const value_t* _values, const double _eps = std::numeric_limits<value_t>::epsilon()) const = 0;
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	protected:
		/// Internal: Assigns all member variables of tensor.
		void assign(const Tensor& _other);
		
		/// Internal: Move assigns all member variables of tensor.
		void assign(Tensor&& _other);
		
		/// Internal: Changes the dimensions of the tensor and recalculates the size of the tensor.
		void change_dimensions(const std::vector<size_t>& _newDimensions);
		
		/// Internal: Changes the dimensions of the tensor and recalculates the size of the tensor.
		void change_dimensions(	  std::vector<size_t>&& _newDimensions);
	};
	
	
	/** 
	* @brief Calculates the frobenius norm of the given tensor
	* @param _tensor the Tensor of which the frobenious norm shall be calculated.
	* @return the frobenius norm .
	*/
	static _inline_ value_t frob_norm(const Tensor& _tensor) { return _tensor.frob_norm(); }
	
	/** 
	* @brief Checks whether two Tensors are approximately equal.
	* @details Check whether ||@a _a - @a _b ||/(||@a a ||/2 + ||@a _b ||/2) < _eps, i.e. whether the relative difference in the frobenius norm is sufficently small.
	* @param _a the first test candidate.
	* @param _b the second test candidate
	* @param _eps the maximal relative difference between @a _a and @a _b.
	* @return TRUE if @a _a and @a _b are determined to be approximately equal, FALSE otherwise.
	*/
	bool approx_equal(const xerus::Tensor& _a, const xerus::Tensor& _b, const xerus::value_t _eps = EPSILON);
	
	/** 
	* @brief Checks whether two Tensors are approximately entrywise equal.
	* @details Check whether |@a _a[i] - @a _b[i] |/(|@a _a[i] | + |@a _b[i] | < _eps is fullfilled for all i.
	* @param _a the first test candidate.
	* @param _b the second test candidate
	* @param _eps the maximal relative difference between the entries of @a _a and @a _b.
	* @return TRUE if @a _a and @a _b are determined to be approximately equal, FALSE otherwise.
	*/
	bool approx_entrywise_equal(const xerus::Tensor& _a, const xerus::Tensor& _b, const xerus::value_t _eps = EPSILON);
}
