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
* @brief Header file for the TTNetwork class (and thus TTTensor and TTOperator).
*/

#pragma once

 
#include "misc/containerSupport.h"
#include "misc/check.h"

#include "index.h"
#include "indexedTensor_tensor_factorisations.h"
#include "tensor.h"
#include "tensorNetwork.h"
#include "indexedTensor.h"
#include "indexedTensorMoveable.h"
#include "indexedTensorList.h"

namespace xerus {
	/**
	* @brief Specialized TensorNetwork class used to represent TTTensor and TToperators.
	* @details TTTensors correspond to isOperator=FALSE and TTOperators correspond to isOperator=FALSE.
	*/
	template<bool isOperator>
	class TTNetwork final : public TensorNetwork {
	public:
		///@brief The number of external links in each node, i.e. one for TTTensors and two for TTOperators.
		static constexpr const size_t N = isOperator?2:1;
		
		/// @brief Flag indicating whether the TTNetwork is cannonicalized.
		bool cannonicalized;
		
		/**
		* @brief The position of the core.
		* @details If cannonicalized is TRUE, corePosition gives the position of the core tensor. All components
		* with smaller index are then left-orthogonalized and all components with larger index right-orthogonalized.
		*/
		size_t corePosition;		
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		* @brief Constructs an order zero TTNetwork.
		* @details This is an empty TensorNetwork. Internally the network contains one order zero node with entry zero.
		*/
		explicit TTNetwork();
		
		
		///@brief TTNetworks are default copy constructable.
		implicit TTNetwork(const TTNetwork & _cpy) = default;
		
		
		///@brief TTNetworks are default move constructable.
		implicit TTNetwork(      TTNetwork&& _mov) = default;
		
		
		/** 
		* @brief Constructs an zero initialized TTNetwork with the given degree and ranks all equal to one.
		* @details Naturally for TTOperators the degree must be even.
		*/
		explicit TTNetwork(const size_t _degree);
		
		
		/** 
		* @brief Constructs a TTNetwork from the given Tensor.
		* @details  The higher order SVD algorithm is used to decompose the given Tensor into the TT format.
		* @param _tensor The Tensor to decompose.
		* @param _eps the accuracy to be used in the decomposition.
		* @param _maxRank the maximal allowed rank (applies to all positions).
		*/
		explicit TTNetwork(const Tensor& _tensor, const double _eps=EPSILON, const size_t _maxRank=std::numeric_limits<size_t>::max());
		
		
		/** 
		* @brief Constructs a TTNetwork from the given Tensor.
		* @details  The higher order SVD algorithm is used to decompose the given Tensor into the TT format.
		* @param _tensor The Tensor to decompose.
		* @param _eps the accuracy to be used in the decomposition.
		* @param _maxRanks maximal ranks to be used
		*/
		explicit TTNetwork(const Tensor& _tensor, const double _eps, const RankTuple& _maxRanks);
		
		
		/** 
		* @brief Transforms a given TensorNetwork to a TTNetwork.
		* @details This is not yet implemented different from casting to Tensor and then using a HOSVD.
		* @param _network The network to transform.
		* @param _eps the accuracy to be used in the transformation.
		*/
		explicit TTNetwork(const TensorNetwork &_network, double _eps=EPSILON);
		
		
		/** 
		 * @brief Random constructs a TTNetwork with the given dimensions and ranks. 
		 * @details The entries of the componend tensors are sampled independendly using the provided random generator and distribution.
		 * @param _dimensions the dimensions of the to be created TTNetwork.
		 * @param _ranks the ranks of the to be created TTNetwork.
		 * @param _rnd the random engine to be passed to the constructor of the component tensors.
		 * @param _dist the random distribution to be passed to the constructor of the component tensors.
		 */
		template<class generator, class distribution>
		static TTNetwork random(const std::vector<size_t>& _dimensions, const std::vector<size_t> &_ranks, generator& _rnd, distribution& _dist) {
			REQUIRE(_dimensions.size()%N==0, "Illegal number of dimensions for ttOperator");
			REQUIRE(_ranks.size() == _dimensions.size()/N-1,"Non-matching amount of ranks given to TTNetwork::random");
			REQUIRE(!misc::contains(_dimensions, 0ul), "Trying to construct a TTTensor with dimension 0 is not possible.");
			REQUIRE(!misc::contains(_ranks, 0ul), "Trying to construct random TTTensor with rank 0 is illegal.");
			
			TTNetwork result(_dimensions.size());
			const size_t numComponents = _dimensions.size()/N;
			const std::vector<size_t> targetRank = reduce_to_maximal_ranks(_ranks, _dimensions);
			
			for(size_t i = 0; i < numComponents; ++i) {
				size_t leftRank = i==0 ? 1 : targetRank[i-1];
				size_t rightRank = i==numComponents-1 ? 1 : targetRank[i];

				if(isOperator) {
					result.set_component(i, Tensor::random({leftRank, _dimensions[i], _dimensions[numComponents+i], rightRank}, _rnd, _dist));
				} else {
					result.set_component(i, Tensor::random({leftRank, _dimensions[i], rightRank}, _rnd, _dist));
				}
			}
			result.cannonicalize_left();
			return result;
		}
		
		
		/** 
		 * @brief Random constructs a TTNetwork with the given dimensions and ranks limited by the given rank. 
		 * @details The entries of the componend tensors are sampled independendly using the provided random generator and distribution.
		 * @param _dimensions the dimensions of the to be created TTNetwork.
		 * @param _ranks the maximal allowed rank. 
		 * @param _rnd the random engine to be passed to the constructor of the component tensors.
		 * @param _dist the random distribution to be passed to the constructor of the component tensors.
		 */
		template<class generator, class distribution>
		static TTNetwork random(const std::vector<size_t>& _dimensions, const size_t _rank, generator& _rnd, distribution& _dist) {
			return TTNetwork::random(_dimensions, std::vector<size_t>(_dimensions.size()/N-1, _rank), _rnd, _dist);
		}
		
		
		/**
		 * @brief Random constructs a TTNetwork with the given dimensions and ranks. 
		 * @details The entries of the componend tensors are sampled independendly using the provided random generator and distribution.
		 *  the singular values of all matricisations M(1..n,n+1..N) are fixed according to the given function a posteriori
		 *  The callback function is assumed to take a reference to a diagonal tensor and modify it to represent the desired singular values.
		 */
		template<class generator, class distribution, class aloc = std::allocator<size_t>>
		static TTNetwork random(const std::vector<size_t, aloc>& _dimensions, const std::vector<size_t> &_ranks, 
								generator& _rnd, distribution& _dist,
								const std::function<void(Tensor&)> &_modifySingularValues) 
		{
			TTNetwork result = random(_dimensions, _ranks, _rnd, _dist);
			
			const Index i,j,k,l,m;
			
			for (size_t pos = 0; pos+1 < result.degree(); ++pos) {
				Tensor A;
				A(i,j^N,k^N,l) = result.component(pos)(i,j^N,m) * result.component(pos+1)(m,k^N,l);
				Tensor U,S,Vt;
				(U(i&1,j),S(j,k),Vt(k,l&1)) = SVD(A(i/2,l/2));
				_modifySingularValues(S);
				Vt(j,l&1) = S(j,k) * Vt(k,l&1);
				result.set_component(pos, U);
				result.set_component(pos+1, Vt);
				result.assume_core_position(pos+1);
			}
			
			REQUIRE(result.is_valid_tt(), "Internal Error.");
			REQUIRE(!result.exceeds_maximal_ranks(), "Internal Error");
			result.cannonicalize_left();
			return result;
		}
		
		
		/** 
		 * @brief: Returns a the (rank one) TT-Tensor with all entries equal to one.
		 * @param _dimensions the dimensions of the new tensor.
		 */
		static TTNetwork ones(const std::vector<size_t>& _dimensions);
		
		
		/** 
		 * @brief: Construct a TTOperator with the given dimensions representing the identity.
		 * @details Only applicable for TTOperators, i.e. not for TTtensors
		 * @param _dimensions the dimensions of the new TTOperator.
		 */
		template<bool B = isOperator, typename std::enable_if<B, int>::type = 0>
		static TTNetwork identity(const std::vector<size_t>& _dimensions);
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard Operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		///@brief TTNetworks are default assignable.
		TTNetwork& operator=(const TTNetwork&  _other) = default;
		
		
		///@brief TTNetworks are default move-assignable.
		TTNetwork& operator=(      TTNetwork&& _other) = default;
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	protected:
		///@brief Constructs a TTNetwork in _out by decomposing the given Tensor _A.
		static void construct_train_from_full(TensorNetwork& _out, const Tensor& _A, const double _eps);
		
		
		/** 
		 * @brief Tests whether any rank exceeds the theoretic maximal value it should have.
		 * @details Does not check for the actual minimal rank for this tensor. But if any rank exceeds the theoretic maximum it is guaranteed not to be the minimal rank.
		 * @return TRUE if any rank exceeds its theoretic maximum.
		 */
		bool exceeds_maximal_ranks() const;
		
		
		///@brief Return the number of ranks, i.e. 0 for degree zero and degree()/N-1 otherwise.
		size_t num_ranks() const;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	public:
		/** 
		* @brief Reduces the given ranks to the maximal possible.
		* @details If a given rank is allready smaller or equal it is left unchanged.
		* @param _ranks the inital ranks to be reduced.
		* @param _dimensions the dimensions used to calculate the maximal ranks.
		* @return the reduced ranks.
		*/
		static std::vector<size_t> reduce_to_maximal_ranks(std::vector<size_t> _ranks, const std::vector<size_t>& _dimensions);
		
		
		/** 
		* @brief Computes the dyadic product of @a _lhs and @a _rhs. 
		* @details This function is currently needed to keep the resulting network in the TTNetwork class.
		* Apart from that the result is the same as result(i^d1, j^d2) = _lhs(i^d1)*_rhs(j^d2).
		* @returns the dyadic product as a TTNetwork.
		*/
		static TTNetwork dyadic_product(const TTNetwork &_lhs, const TTNetwork &_rhs);
		
		
		/** 
		* @brief Computes the dyadic product of all given TTNetworks. 
		* @details This is nothing but the repeated application of dyadic_product() for the given TTNetworks.
		* @returns the dyadic product as a TTNetwork.
		*/
		static TTNetwork dyadic_product(const std::vector<std::reference_wrapper<TTNetwork>> &_tensors);
		
		
		/**
		* @brief Calculates the componentwise product of two tensors given in the TT format.
		* @details In general the resulting rank = rank(A)*rank(B). Retains the core position of @a _A
		*/
		static TTNetwork entrywise_product(const TTNetwork &_A, const TTNetwork &_B);
		
		
		/**
		* @brief Computes the entrywise square of the tensor.
		* @details In general the resulting equals rank*(rank+1)/2. Retains the core position.
		*/
		void entrywise_square();
		
		
		/** 
		* @brief Complete access to a specific component of the TT decomposition.
		* @note This function will not update rank and external dimension informations if it is used to set a component.
		* @details This function gives complete access to the components, only intended for internal use.
		* @param _idx index of the component to access.
		* @returns a reference to the requested component.
		*/
		Tensor& component(const size_t _idx);
		
		
		/** 
		* @brief Read access to a specific component of the TT decomposition.
		* @details This function should be used to access the components, instead of direct access via
		* nodes[...], because the implementation does not store the first component in nodes[0] but rather as
		* nodes[1] etc. nodes[0] is an order one node with dimension one only used to allow the first component
		* to be an order three tensor.
		* @param _idx index of the component to access.
		* @returns a const reference to the requested component.
		*/
		const Tensor& get_component(const size_t _idx) const;
		
		
		/** 
		* @brief Sets a specific component of the TT decomposition.
		* @details This function also takes care of adjusting the corresponding link dimensions and external dimensions
		* if needed. However this might still leave the TTNetwork in an invalid if the rank is changed. In this case it
		* is the callers responsibility to also update the other component tensors consistently to account for that rank
		* change.
		* @param _idx index of the component to set.
		* @param _T Tensor to use as the new component tensor.
		*/
		void set_component(const size_t _idx, Tensor _T);
		
		
		/** 
		* @brief Splits the TTNetwork into two parts by removing the node.
		* @param _position index of the component to be removed, thereby also defining the position 
		*of the split.
		* @return a std::pair containing the two remaining parts as TensorNetworks.
		*/
		std::pair<TensorNetwork, TensorNetwork> chop(const size_t _position) const;
		
		
		/** 
		* @brief Reduce all ranks up to a given accuracy and maximal number.
		* @param _maxRanks maximal allowed ranks. All current ranks that are larger than the given ones are reduced by truncation.
		* @param _eps the accuracy to use for truncation in the individual SVDs.
		*/
		void round(const std::vector<size_t>& _maxRanks, const double _eps = EPSILON);
		
		
		/** 
		* @brief Reduce all ranks to the given number.
		* @param _maxRank maximal allowed rank. All current ranks that are larger than this are reduced by truncation.
		*/
		void round(const size_t _maxRank);
		
		
		/** 
		* @brief Reduce all ranks to the given number.
		* @param _maxRank maximal allowed rank. All current ranks that are larger than this are reduced by truncation.
		*/
		void round(const int _maxRank);
		
		
		/** 
		* @brief Reduce all ranks up to a given accuracy.
		* @param _eps the accuracy to use for truncation in the individual SVDs.
		*/
		void round(const value_t _eps);
		
		
		/** 
		* @brief Applies the soft threshholding operation to all ranks.
		* @param _tau the soft threshholding parameter to be applied. I.e. all singular values are reduced to max(0, Lambda_ui - _tau).
		*/
		void soft_threshold(const double _tau, const bool _preventZero = false);
		
		
		/** 
		* @brief Applies soft threshholding operations to all ranks.
		* @param _taus the soft threshholding parameters to be applied. I.e. all singular values of the j-th matrification are reduced to max(0, Lambda_ui - _tau[j]).
		*/
		void soft_threshold(const std::vector<double>& _taus, const bool _preventZero = false);
		
		
		/** 
		* @brief Gets the ranks of the TTNetwork.
		* @return A vector containing the current ranks.
		*/
		std::vector<size_t> ranks() const;
		
		
		/** 
		* @brief Gets the rank of a specific egde of the TTNetwork.
		* @param _i Position of the edge in question.
		* @return The current rank of edge _i.
		*/
		size_t rank(const size_t _i) const;
		
		
		/** 
		* @brief Calculates the storage requirement of the current representation.
		* @return The datasize in sizeof(value_t).
		*/
		size_t datasize() const;
		
		
		/** 
		* @brief Move the core to a new position.
		* @details The core is moved to @a _position and the nodes between the old and the new position are orthogonalized
		* accordingly. If the TTNetwork is not yet cannonicalized it will be with @a _position as new corePosition.
		* @param _position the new core position.
		* @param _keepRank by default a rank revealing QR decomposition is used to move the core and the ranks are reduced
		* accordingly. If @a _keepRank is set the rank is not reduced, this is need e.g. in the ALS.
		*/
		void move_core(const size_t _position, const bool _keepRank=false);
		
		
		/**
		* @brief stores @a _pos as the current core position without verifying of ensuring that this is the case
		* @details this is particularly useful after constructing an own TT tensor with set_component calls
		* as these will assume that all orthogonalities are destroyed
		*/
		void assume_core_position(const size_t _pos);
		
		
		/** 
		* @brief Move the core to the left.
		* @details Basically calls move_core() with _position = 0
		*/
		void cannonicalize_left();
		
		
		/** 
		* @brief Move the core to the left.
		* @details Basically calls move_core() with _position = degree()-1
		*/
		void cannonicalize_right();
		
		
		/** 
		* @brief Transpose the TTOperator
		* @details Swaps all external indices to create the transposed operator.
		*/
		template<bool B = isOperator, typename std::enable_if<B, int>::type = 0>
		void transpose() {
			const std::vector<size_t> shuffle({0,2,1,3});
			for (size_t n = 0; n < degree()/N; ++n) {
				xerus::reshuffle(component(n), component(n), shuffle);
			}
		}
		
		
		virtual TensorNetwork* get_copy() const override;
		
		
		virtual value_t frob_norm() const override;
		
		
		/** 
		* @brief Finds the position of the approximately largest entry.
		* @details Uses an algorithms to find an entry that is at least of size @a _accuracy * X_max in absolute value,
		* where X_max is the largest entry of the tensor. The smaller @a _accuracy, the faster the algorithm will work.
		* @param _accuracy factor that determains the maximal deviation of the returned entry from the true largest entry.
		* @param _lowerBound a lower bound for the largest entry, i.e. there must be an entry in the tensor which is at least of
		* this size (in absolute value). The algorithm may fail completely if this is not fullfilled, but will work using its own 
		* approximation if no value (i.e. 0.0) is given.
		* @return the position of the entry found.
		*/
		size_t find_largest_entry(const double _accuracy, const value_t _lowerBound = 0.0) const;
		
		
		/** 
		 * @brief Tests whether the network resembles that of a TTTensor and checks consistency with the underlying tensor objects.
		 * @details Note that this will NOT check for orthogonality of cannonicalized TTNetworks.
		 */
		virtual void require_correct_format() const override;
		
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - -  Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - */
		/** 
		* @brief Adds a given TTNetwork to this one.
		* @details To be well-defined it is required that the dimensions of this and @a _other coincide. 
		* The rank of the result are in general the entrywise sum of the ranks of this and @a _other.
		* @param _other the TTNetwork to add.
		* @return reference to this TTNetwork.
		*/
		TTNetwork& operator+=(const TTNetwork& _other);
		
		
		/** 
		* @brief Subtracts the @a _other TTNetwork entrywise from this one.
		* @details To be well-defined it is required that the dimensions of this and @a _other coincide. 
		* The rank of the result are in general the entrywise sum of the ranks of this and @a _other.
		* @param _other the Tensor to be subtracted to this one.
		* @return a reference to this TTNetwork.
		*/
		TTNetwork& operator-=(const TTNetwork& _other);
		
		/** 
		* @brief Calculates the entrywise multiplication of this TensorNetwork with a constant @a _factor.
		* @details Internally this only results in a change in the global factor.
		* @param _factor the factor.
		*/
		virtual void operator*=(const value_t _factor) override;
		
		/** 
		* @brief Calculates the entrywise divison of this TensorNetwork by a constant @a _divisor.
		* @details Internally this only results in a change in the global factor.
		* @param _divisor the divisor.
		*/
		virtual void operator/=(const value_t _divisor) override;
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
		static bool specialized_contraction_f(std::unique_ptr<internal::IndexedTensorMoveable<TensorNetwork>>& _out, internal::IndexedTensorReadOnly<TensorNetwork>&& _me, internal::IndexedTensorReadOnly<TensorNetwork>&& _other);
		
		static bool specialized_sum_f(std::unique_ptr<internal::IndexedTensorMoveable<TensorNetwork>>& _out, internal::IndexedTensorReadOnly<TensorNetwork>&& _me, internal::IndexedTensorReadOnly<TensorNetwork>&& _other);
		
		virtual bool specialized_contraction(std::unique_ptr<internal::IndexedTensorMoveable<TensorNetwork>>& _out, internal::IndexedTensorReadOnly<TensorNetwork>&& _me, internal::IndexedTensorReadOnly<TensorNetwork>&& _other) const override {
			return specialized_contraction_f(_out, std::move(_me), std::move(_other));
		}
		
		virtual bool specialized_sum(std::unique_ptr<internal::IndexedTensorMoveable<TensorNetwork>>& _out, internal::IndexedTensorReadOnly<TensorNetwork>&& _me, internal::IndexedTensorReadOnly<TensorNetwork>&& _other) const override {
			return specialized_sum_f(_out, std::move(_me), std::move(_other));
		}
		
		virtual void specialized_evaluation(internal::IndexedTensorWritable<TensorNetwork>&& _me, internal::IndexedTensorReadOnly<TensorNetwork>&& _other) override;
		
	};

	typedef TTNetwork<false> TTTensor;
	typedef TTNetwork<true> TTOperator;
	

	/** 
	* @brief Calculates the entrywise sum of the given TTNetworks @a _lhs and @a _rhs.
	* @details To be well-defined it is required that the dimensions of @a _lhs and @a _rhs coincide. 
	* The rank of the result are in general the entrywise sum of the ranks of @a _lhs and @a _rhs.
	* @param _lhs the first summand.
	* @param _rhs the second summand.
	* @return the sum.
	*/
	template<bool isOperator>
	TTNetwork<isOperator> operator+(TTNetwork<isOperator> _lhs, const TTNetwork<isOperator>& _rhs);


	/** 
	* @brief Calculates the entrywise difference of the given TTNetworks @a _lhs and @a _rhs.
	* @details To be well-defined it is required that the dimensions of @a _lhs and @a _rhs coincide. 
	* The rank of the result are in general the entrywise sum of the ranks of @a _lhs and @a _rhs.
	* @param _lhs the minuend.
	* @param _rhs the subtrahend.
	* @return the difference.
	*/
	template<bool isOperator>
	TTNetwork<isOperator> operator-(TTNetwork<isOperator> _lhs, const TTNetwork<isOperator>& _rhs);
	
	
	/** 
	* @brief Calculates the entrywise multiplication of the given TTNetwork @a _network with a constant @a _factor.
	* @details Internally this only results in a change in the global factor.
	* @param _network the TTNetwork,
	* @param _factor the factor,
	* @return the resulting scaled TTNetwork.
	*/
	template<bool isOperator>
	TTNetwork<isOperator> operator*(TTNetwork<isOperator> _network, const value_t _factor);
	
	
	/** 
	* @brief Calculates the entrywise multiplication of the given TTNetwork @a _network with a constant @a _factor.
	* @details Internally this only results in a change in the global factor.
	* @param _factor the factor,
	* @param _network the TTNetwork,
	* @return the resulting scaled TTNetwork.
	*/
	template<bool isOperator>
	TTNetwork<isOperator> operator*(const value_t _factor, TTNetwork<isOperator> _network);


	/** 
	* @brief Calculates the entrywise divison of this TTNetwork by a constant @a _divisor.
	* @details Internally this only results in a change in the global factor.
	* @param _divisor the divisor,
	* @return the resulting scaled TTNetwork.
	*/
	template<bool isOperator>
	TTNetwork<isOperator> operator/(const TTNetwork<isOperator>& _network, const value_t _div);
}
