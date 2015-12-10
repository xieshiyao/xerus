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
	class Tensor;
	
	/// @brief Class that handles simple (non-decomposed) tensors in a dense or sparse representation.
	class Tensor {
	public:
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Auxiliary types- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/** 
		 * @brief Flags determining the initialisation of the data of Tensor objects. 
		 * @details None means that no initialisation is performed, i.e. the data can be random.
		 * Zero means that the data is zero initialized.
		 */
		enum class Initialisation : bool { None, Zero};
		
		/** 
		 * @brief Flags indicating the internal representation of the data of Tensor objects. 
		 * @details Dense means that an value_t array of 'size' is used to store each entry individually, 
		 * using row-major order. Sparse means that only the non-zero entries are stored explicitly in a set containing
		 * their value and position.
		 */
		enum class Representation : bool { Dense, Sparse};
		
		///@brief: Represention of the dimensions of a Tensor.
		typedef std::vector<size_t> DimensionTuple;
		
		///@brief: Represention of a MultiIndex, i.e. the tuple of positions for each dimension determining a single position in a Tensor.
		typedef std::vector<size_t> MultiIndex;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/// @brief Vector containing the individual dimensions of the tensor.
		DimensionTuple dimensions;
		
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
		 * @brief Shared pointer to the a map containing the non-zero entries, if representation is Sparse. 
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
// 		implicit Tensor( Tensor&& _other );

		/** 
		 * @brief: Creates a new tensor with the given dimensions.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _representation (optional) the initial representation of the tensor.
		 * @param _init (optional) inital data treatment, i.e. whether the tensor is to be zero Initialized.
		 */
		explicit Tensor(const DimensionTuple& _dimensions, const Representation _representation = Representation::Dense, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief: Creates a new tensor with the given dimensions.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _representation (optional) the initial representation of the tensor.
		 * @param _init (optional) inital data treatment, i.e. whether the tensor is to be zero Initialized.
		 */
		explicit Tensor(      DimensionTuple&& _dimensions, const Representation _representation = Representation::Dense, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief: Creates a new (dense) tensor with the given dimensions, using a provided data.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _data inital dense data in row-major order.
		 */
		template<ADD_MOVE(Vec, DimensionTuple), ADD_MOVE(SPtr, std::shared_ptr<value_t>)>
		explicit Tensor(Vec&& _dimensions, SPtr&& _data)
		: dimensions(std::forward<Vec>(_dimensions)), size(misc::product(dimensions)), representation(Representation::Dense), denseData(std::forward<SPtr>(_data)) { }
		
		/** 
		 * @brief: Creates a new (dense) tensor with the given dimensions, using a provided data.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _data inital dense data in row-major order.
		 */
		explicit Tensor(const DimensionTuple& _dimensions, std::unique_ptr<value_t[]>&& _data);
		
		/** 
		 * @brief Constructs a Tensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details In this overload no value is passed to _f, i.e. _f must determine the values of the entries independend of their position,
		 * or keep track of the position itself. _f may assume that it is called for the entries in the order they are stored (i.e. row-major order)
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to use to set the entries of the Tensor. 
		 */
		explicit Tensor(const DimensionTuple& _dimensions, const std::function<value_t()>& _f);
		
		/** 
		 * @brief Constructs a Tensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details In this overload the position of each entry assuming row-major order is passed to  _f.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to use to set the entries of the Tensor. 
		 */
		explicit Tensor(const DimensionTuple& _dimensions, const std::function<value_t(const size_t)>& _f);
		
		/** 
		 * @brief Constructs a Tensor with the given dimensions and uses the given function to assign the values to the entries.
		 * @details In this overload the complete position of each entry is passed to  _f.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to use to set the entries of the Tensor. 
		 */
		explicit Tensor(const DimensionTuple& _dimensions, const std::function<value_t(const MultiIndex&)>& _f);
		
		/** 
		 * @brief Constructs a Tensor with the given dimensions and uses the given function @a _f to create @a _N non zero entries.
		 * @details @a _f is called with the current number of entries present and the number of possible entries (i.e. size). @a _f shall return a pair containg the position
		 * and value of the next entry. @a _f is required not to return a position twice.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _f the function to be used to create each non zero entry. 
		 * @param _N the number of non-zero entries to be created.
		 */
		Tensor(const DimensionTuple& _dimensions, std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f, const size_t _N);
		
		
		
		/** 
		 * @brief Constructs a dense Tensor with the given dimensions and uses the given random generator and distribution to assign the values to the entries.
		 * @details The entries are assigned in the order they are stored (i.e. row-major order). Each assigned is a seperate call to the random distribution.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _rnd the random generator to be used.
		 * @param _dist the random distribution to be used.
		 */
		template<ADD_MOVE(Dim_T, DimensionTuple), class generator, class distribution>
		static Tensor random(Dim_T&& _dimensions, generator& _rnd, distribution& _dist) {
			Tensor result(std::forward<Dim_T>(_dimensions), Representation::Dense, Initialisation::None);
			value_t* const dataPtr = result.denseData.get();
			for(size_t i = 0; i < result.size; ++i) {
				dataPtr[i] = _dist(_rnd);
			}
			return result;
		}
		
		/** 
		 * @brief Constructs a dense Tensor with the given dimensions and uses the given random generator and distribution to assign the values to the entries.
		 * @details See the std::vector variant for details.
		 */
		template<class generator, class distribution>
		_inline_ static Tensor random(std::initializer_list<size_t>&& _dimensions, generator& _rnd, distribution& _dist) {
			return Tensor::random(std::vector<size_t>(std::move(_dimensions)), _rnd, _dist);
		}
		
		/** 
		 * @brief Constructs a random sparse Tensor with the given dimensions.
		 * @details The given random generator @a _rnd and distribution @a _dist are used to assign the values to @a _n randomly choosen entries.
		 * @param _dimensions the future dimensions of the Tensor.
		 * @param _N the number of non-zero entries to be created.
		 * @param _rnd the random generator to be used.
		 * @param _dist the random distribution to be used.
		 */
		template<ADD_MOVE(Vec, std::vector<size_t>), class generator, class distribution>
		static Tensor random(Vec&& _dimensions, const size_t _N, generator& _rnd, distribution& _dist) {
			Tensor result(std::forward<Vec>(_dimensions), Representation::Sparse, Initialisation::Zero);
			REQUIRE(_N <= result.size, " Cannot create " << _N << " non zero entries in a tensor with only " << result.size << " total entries!");
			
			std::uniform_int_distribution<size_t> entryDist(0, result.size-1);
			while(result.sparseData->size() < _N) {
				result.sparseData->emplace(entryDist(_rnd), _dist(_rnd));
			}
			return result;
		}
		
		/** 
		 * @brief Constructs a random sparse Tensor with the given dimensions.
		 * @details See the std::vector variant for details.
		 */
		template<class generator, class distribution>
		_inline_ static Tensor random(std::initializer_list<size_t>&& _dimensions, const size_t _N, generator& _rnd, distribution& _dist) {
			return Tensor::random(std::vector<size_t>(_dimensions), _N, _rnd, _dist);
		}
		
		
		
		/// @brief Returns a copy of this Tensor that uses a dense representation.
		Tensor dense_copy() const;
		
		/// @brief Returns a copy of this Tensor that uses a sparse representation.
		Tensor sparse_copy() const;
		
		
		/// @brief Returns a pointer containing a copy of the tensor with same type (i.e. Tensor or SparseTensor).
		Tensor* get_copy() const;
		
		
		/** 
		 * @brief: Returns a Tensor with all entries equal to one.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		static Tensor ones(const std::vector<size_t>& _dimensions);
		
		/** 
		 * @brief: Returns a Tensor representation of the identity operator with the given dimensions.
		 * @details That is combining the first half of the dimensions and the second half of the dimensions results in an identity matrix.
		 * @param _dimensions the dimensions of the new tensor. It is required that _dimensions[i] = _dimensions[d/2+i], otherwise this cannot be the identity operator.
		 */
		static Tensor identity(const std::vector<size_t>& _dimensions);
		
		/** 
		 * @brief: Returns a Tensor representation of the kronecker delta.
		 * @details That is each entry is one if all indices are equal and zero otherwise. Note iff d=2 this coincides with identity.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		static Tensor kronecker(const std::vector<size_t>& _dimensions);
		
		/** 
		 * @brief: Returns a Tensor with a single entry equals oen and all other zero.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _position The position of the one
		 */
		static Tensor dirac(const std::vector<size_t>& _dimensions, const std::vector<size_t>& _position);
		
		/** 
		 * @brief: Returns a Tensor with a single entry equals oen and all other zero.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _position The position of the one
		 */
		static Tensor dirac(const std::vector<size_t>& _dimensions, const size_t _position);
		
		
		
		
		/// @brief Destructor
		virtual ~Tensor();
		
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Standard assignment operator.
		 * @param _other the Tensor to be assinged to this one.
		 * @return a reference to this Tensor.
		 */
		Tensor& operator=(const Tensor&  _other) = default;
		
		/** 
		 * @brief Standard move-assignment operator.
		 * @param _other the Tensor to be move-assinged to this one.
		 * @return a reference to this Tensor.
		 */
// 		Tensor& operator=(      Tensor&& _other);
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Information - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/// @brief Returns whether the current representation is dense.
		bool is_dense() const;
		
		/// @brief Returns whether the current representation is sparse.
		bool is_sparse() const;
		
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		static size_t multiIndex_to_position(const std::vector<size_t>& _multiIndex, const std::vector<size_t>& _dimensions);
		
		template<int sign>
		static void plus_minus_equal(Tensor& _me, const Tensor& _other);
		
		/// @brief Adds the given sparse data to the given full data
		static void add_sparse_to_full(const std::shared_ptr<value_t>& _denseData, const value_t _factor, const std::shared_ptr<const std::map<size_t, value_t>>& _sparseData);
		
		/// @brief Adds the given sparse data to the given sparse data
		static void add_sparse_to_sparse(const std::shared_ptr<std::map<size_t, value_t>>& _sum, const value_t _factor, const std::shared_ptr<const std::map<size_t, value_t>>& _summand);
		
		/// @brief Ensures that this tensor is the sole owner of its data. If needed new space is allocated and all entries are copied.
		void ensure_own_data();
		
		/// @brief Ensures that this tensor is the sole owner of its data space. If needed new space is allocated with entries left undefined.
		void ensure_own_data_no_copy();
		
		/// @brief Checks whether there is a non-trivial scaling factor and applies it if nessecary.
		void apply_factor();
		
		/// @brief Checks whether there is a non-trivial factor and applies it. Even if no factor is applied ensure_own_data() is called.
		void ensure_own_data_and_apply_factor();
		
		/** 
		 * @brief Resets the tensor to the given dimensions and representation.
		 * @details Leaves the Tensor in the same state as if newly constructed with the the same arguments.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _representation the new representation of the tensor.
		 * @param _init (optional) data treatment, i.e. whether the tensor shall be zero initialized.
		 */
		void reset(const std::vector<size_t>&  _newDim, const Representation _representation, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief Resets the tensor to the given dimensions, preserving the current representation.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _init (optional) data treatment, i.e. whether the tensor shall be zero initialized.
		 */
		void reset(const std::vector<size_t>&  _newDim, const Initialisation _init = Initialisation::Zero);
		
		/** 
		 * @brief Resets the tensor to the given dimensions and uses the given data.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _newData new dense data in row-major order.
		 */
		void reset(const std::vector<size_t>&  _newDim, const std::shared_ptr<value_t>& _newData);
		
		/** 
		 * @brief Resets the tensor to the given dimensions, preserving the current representation.
		 * @param _dimensions the dimensions of the new tensor.
		 * @param _newData new dense data in row-major order.
		 */
		void reset(const std::vector<size_t>&  _newDim, std::unique_ptr<value_t[]>&& _newData);
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Adds the @a _other Tensor entrywise to this one.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the Tensor to be added to this one.
		 * @return a reference to this Tensor.
		 */
		Tensor& operator+=(const Tensor& _other);
		
		/** 
		 * @brief Subtracts the @a _other Tensor entrywise from this one.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the Tensor to be subtracted to this one.
		 * @return a reference to this Tensor.
		 */
		Tensor& operator-=(const Tensor& _other);
		
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
		
		/** 
		 * @brief Calculates the entrywise sum of this Tensor and @a _other.
		 * @details To be well-defined it is required that the dimensions of this and @a _other coincide.
		 * @param _other the second summand.
		 * @return the sum.
		 */
		Tensor  operator+( const Tensor& _other) const;

		/** 
		 * @brief Calculates the entrywise difference between this Tensor and @a _other.
		 * @details To be well-defined it is required that the dimensions of this and _other coincide.
		 * @param _other the subtrahend,
		 * @return the difference.
		 */
		Tensor  operator-( const Tensor& _other) const;
		
		/** 
		 * @brief Calculates the entrywise multiplication of this Tensor with a constant @a _factor.
		 * @details Internally this only results in a change in the global factor.
		 * @param _factor the factor,
		 * @return the resulting scaled Tensor.
		 */
		Tensor  operator*( const value_t _factor) const;
		
		/** 
		 * @brief Calculates the entrywise divison of this Tensor by a constant @a _divisor.
		 * @details Internally this only results in a change in the global factor.
		 * @param _divisor the divisor,
		 * @return the resulting scaled Tensor.
		 */
		Tensor  operator/( const value_t _divisor) const;
	
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/** 
		 * @brief Read/Write access a single entry.
		 * @details Note that the request for write access to this entry requires the a sparse representation to create it, i.e. to lower the sparsity. If only read access
		 * is required use at() instead.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return a reference to the selected entry.
		 */
		value_t& operator[](const size_t _position);
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return the value of the selected entry.
		 */
		value_t operator[](const size_t _position) const;
		
		/** 
		 * @brief Read/Write access a single entry.
		 * @param _positions the positions of the desired entry.
		 * @return a reference to the selected entry.
		 */
		value_t& operator[](const MultiIndex& _positions);
		
		/** 
		 * @brief Read access a single entry.
		 * @param _positions the positions of the desired entry.
		 * @return the value of the selected entry.
		 */
		value_t operator[](const MultiIndex& _positions) const;
		
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry, assuming row-major ordering.
		 * @return the value of the selected entry.
		 */
		value_t at(const size_t _position) const;
		
		/** 
		 * @brief Read access a single entry.
		 * @param _position the position of the desired entry.
		 * @return the value of the selected entry.
		 */
		value_t at(const MultiIndex& _positions) const;
		
		
		/** 
		 * @brief Returns a pointer for direct access to the dense data array in row major order. 
		 * @details Also takes care that this direct access is safe, i.e. that this tensor is using a dense representation, is the sole owner of the data and that no non trivial factor exists.
		 * @return pointer to the dense data array.
		 */
		value_t* get_dense_data();
		
		/** 
		 * @brief Gives access to the internal data pointer, without any checks.
		 * @details Note that the dense data array might not exist because a sparse representation is used, may shared with other tensors 
		 * or has to be interpreted considering a gloal factor. Both can be avoid if using data_pointer(). The tensor data itself is stored in row-major ordering.
		 * @return pointer to the internal dense data array.
		 */
		value_t* get_unsanitized_dense_data();
		
		/** 
		 * @brief Gives access to the internal data pointer, without any checks.
		 * @details Note that the dense data array might not exist because a sparse representation is used, may shared with other tensors 
		 * or has to be interpreted considering a gloal factor. Both can be avoid if using data_pointer(). The tensor data itself is stored in row-major ordering.
		 * @return pointer to the internal dense data array.
		 */
		const value_t* get_unsanitized_dense_data() const;
		
		/** 
		 * @brief Returns a pointer to the internal dense data array for complete rewrite purpose ONLY.
		 * @details This is equivalent to calling reset() with the current dimensions, dense representation and no initialisation and then
		 * calling get_unsanitized_dense_data().
		 * @return pointer to the internal dense data array.
		 */
		value_t* override_dense_data();
		
		/** 
		 * @brief Gives access to the internal shared data pointer, without any checks.
		 * @details Note that the data array might be shared with other tensors or has to be interpreted considering a global
		 * factor. Both can be avoid if using data_pointer(). The tensor data itself is stored in row-major ordering.
		 * @return The internal shared pointer to the data array.
		 */
		const std::shared_ptr<value_t>& get_internal_dense_data();
		
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
		
		
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Modifiers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
		/** 
		 * @brief Resizes a specific dimension of the Tensor.
		 * @param _n the dimension to resize.
		 * @param _newDim the new value that resized dimension shall have.
		 * @param _cutPos the index within the selected dimension after which new slates are inserted or removed, by default the last index.
		 */
		void resize_dimension(const size_t _n, const size_t _newDim, size_t _cutPos=~0ul);
		
		/** 
		 * @brief Removes a single slate from the Tensor.
		 * @param _indexNb the dimension to in defining the slate.
		 * @param _pos the index within the selected dimension for which the slate shall be removed.
		 */
		void remove_slate(const size_t _indexNb, const size_t _pos);
		
		/** 
		 * @brief Modifies the diagonal entries according to the given function.
		 * @details In this overload only the current diagonal entries are passed to @a _f, one at a time. At the moment this is only defined for matricies.
		 * @param _f the function to call to modify each entry.
		 */
		void modify_diag_elements(const std::function<void(value_t&)>& _f);
		
		/** 
		 * @brief Modifies the diagonal entries according to the given function.
		 * @details In this overload the current diagonal entries are passed to @a _f, one at a time, together with their position on the diagonal. At the moment this is only defined for matricies.
		 * @param _f the function to call to modify each entry.
		 */
		void modify_diag_elements(const std::function<void(value_t&, const size_t)>& _f);
		
		/** 
		 * @brief Modifies every entry according to the given function.
		 * @details In this overload only the current entry is passed to @a _f.
		 * @param _f the function to call to modify each entry.
		 */
		void modify_elements(const std::function<void(value_t&)>& _f);
		
		/** 
		 * @brief Modifies every entry according to the given function.
		 * @details In this overload the current entry together with its position, assuming row-major ordering is passed to @a _f.
		 * @param _f the function to call to modify each entry.
		 */
		void modify_elements(const std::function<void(value_t&, const size_t)>& _f);
		
		/** 
		 * @brief Modifies every entry according to the given function.
		 * @details In this overload the current entry together with its complete position is passed to @a _f.
		 * @param _f the function to call to modify each entry.
		 */
		void modify_elements(const std::function<void(value_t&, const std::vector<size_t>&)>& _f);
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		
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
		size_t count_non_zero_entries(const value_t _eps = std::numeric_limits<value_t>::epsilon()) const;
		
		/** 
		 * @brief Checks the tensor for illegal entries, e.g. nan, inf,...
		 * @return TRUE there are no invalid entries, FALSE otherwise.
		 */
		bool all_entries_valid() const;
		
		/** 
		 * @brief Approximates the cost to reorder the tensor.
		 * @return the approximated costs.
		 */
		size_t reorder_costs() const;

		
		/** 
		 * @brief Calculates the frobenious norm of the tensor.
		 * @return the frobenious norm.
		 */
		value_t frob_norm() const;
		
		
		/** 
		 * @brief Makes the Tensor use a dense representation.
		 */
		void use_dense_representation();
		
		/** 
		 * @brief Makes the Tensor use a sparse representation.
		 */
		void use_sparse_representation(const value_t _eps = std::numeric_limits<value_t>::epsilon());
		
		
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
		 * @brief Fixes a specific slate in one of the dimensions, effectively reducing the order by one.
		 * @param _dimension the dimension in which the slate shall be fixed, e.g. 0 to fix the first dimensions.
		 * @param _slatePosition the position in the corresponding dimensions that shall be used.
		 */
		void fix_slate(const size_t _dimension, const size_t _slatePosition);
		
		/** 
		 * @brief Creates a string representation of the Tensor.
		 * @return the string representation. 
		 */
		std::string to_string() const;
		
		
		/** 
		 * @brief Compares the Tensor entriewise to the given data.
		 * @param _values data to compare to, must be of the same size as the Tensor has entries (i.e. size).
		 * @param _eps (optional) epsilon used to determine whether an entry is equal to a data point. 
		 * @return TRUE if all entries deviate by less than _eps from the correspondign data points, FALSE otherwise.
		 */
		bool compare_to_data(const std::vector<value_t>& _values, const double _eps = std::numeric_limits<value_t>::epsilon()) const;
		
		
		/** 
		 * @brief Compares the Tensor entriewise to the given data.
		 * @param _values data to compare to, must be of the same size as the Tensor has entries (i.e. size).
		 * @param _eps (optional) epsilon used to determine whether an entry is equal to a data point. 
		 * @return TRUE if all entries deviate by less than _eps from the correspondign data points, FALSE otherwise.
		 */
		bool compare_to_data(const value_t* _values, const double _eps = std::numeric_limits<value_t>::epsilon()) const;
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	protected:
		/// Internal: Changes the dimensions of the tensor and recalculates the size of the tensor.
		void change_dimensions(const std::vector<size_t>& _newDimensions);
		
		/// Internal: Changes the dimensions of the tensor and recalculates the size of the tensor.
		void change_dimensions(	  std::vector<size_t>&& _newDimensions);
	};
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - External functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	/** 
	 * @brief Calculates the entrywise multiplication of the Tensor @a _lhs with a constant @a _rhs.
	 * @details Internally this only results in a change in the global factor.
	 * @param _lhs the Tensor that shall be scaled.
	 * @param _rhs the factor to be used.
	 * @return the resulting scaled Tensor.
	 */
	static _inline_ Tensor operator*(const value_t _lhs, const Tensor& _rhs) { return _rhs*_lhs; }
	
	
	/** 
	* @brief Calculates the frobenius norm of the given tensor
	* @param _tensor the Tensor of which the frobenious norm shall be calculated.
	* @return the frobenius norm .
	*/
	static _inline_ value_t frob_norm(const Tensor& _tensor) { return _tensor.frob_norm(); }
	
	
	
	/** 
	 * @brief Low-level contraction between Tensors.
	 * @details To be well-defined it is required that the dimensions of @a _lhs and @a _rhs coincide.
	 * @param _result Output for the result of the contraction. Must allready have the right dimensions!
	 * @param _lhs left hand side of the contraction.
	 * @param _lhsTrans Flags whether the LHS should be transposed (in the matrifications sense).
	 * @param _rhs right hand side of the contraction.
	 * @param _rhsTrans Flags whether the RHS should be transposed (in the matrifications sense).
	 * @param _numIndices number of indices that shall be contracted.
	 */
	void contract(Tensor& _result, const Tensor& _lhs, const bool _lhsTrans, const Tensor& _rhs, const bool _rhsTrans, const size_t _numIndices);
	
	
	/**
	 * @brief calculates the entrywise product of two Tensors
	 */
	Tensor entrywise_product(const Tensor &_A, const Tensor &_B);
	
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
