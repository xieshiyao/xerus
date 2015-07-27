// Xerus - A General Purpose Tensor Library
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
 * @brief Header file for the TensorNetwork class.
 */

#pragma once

#include "indexedTensor.h"
#include "tensor.h"
#include "misc/missingFunctions.h"

#include <map>
#include <set>
#include <memory>

namespace xerus {
    // Necessary forward declaritons
    class FullTensor;
    class SparseTensor;
    
	/** 
	* @brief Very general class used to represent arbitary tensor networks.
	* @details Used as a basis for tensor decompositions like the TTNetwork but also used for the lazy evaluation of Tensor contractions.
	*/
	class TensorNetwork {
    public:
		
		/**
		* @brief Class representing a link from a TensorNode to another node or an external index.
		*/
		class Link {
		public:
			///@brief The index of the otherNode this Link links to.
			size_t other; 
			
			///@brief IndexPosition on the other node or index of external index.
			size_t indexPosition;
			
			///@brief Dimension of the link, always equals to other->tensorObject->dimensions[indexPosition].
			size_t dimension;
			
			///@brief Flag indicating whether this link correspond to an external index.
			bool external;
			
			Link() = default;
			Link(const Link& ) = default;
			Link(      Link&&) = default;
			
			Link(const size_t _other, const size_t _indexPos, const size_t _dim, const bool _external);
			
			Link& operator=(const Link& ) = default;
			Link& operator=(      Link&&) = default;
			
			/**
			 * @brief Checks whether this link links to a particular node
			 * @param _other the other node for which the linkage shall be checked
			 * @return TRUE if _other is the target of this Link, FALSE otherwise.
			 */
			bool links(const size_t _other) const;
		};
			
		/**
		* @brief The TensorNode class is used by the class TensorNetwork to store the componentent tensors defining the network.
		*/
		class TensorNode {
		public:
			///@brief Save slot for the tensorObject associated with this node.
			std::unique_ptr<Tensor> tensorObject;
			
			///@brief Vector of links defining the connection of this node to the network.
			std::vector<Link> neighbors;
			
			///@brief Internal Flag
			bool erased;
			
			explicit TensorNode();
			
			implicit TensorNode(const TensorNode&  _other);
			implicit TensorNode(      TensorNode&& _other);
			
			explicit TensorNode(      std::unique_ptr<Tensor>&& _tensorObject);
			
			explicit TensorNode(std::unique_ptr<Tensor>&& _tensorObject, const std::vector<Link>& _neighbors);
			explicit TensorNode(std::unique_ptr<Tensor>&& _tensorObject,       std::vector<Link>&& _neighbors);
			
			TensorNode& operator=(const TensorNode&  _other);
			TensorNode& operator=(      TensorNode&& _other);

			TensorNode strippped_copy() const;
			
			void add_factor(const value_t _factor);
			
			// All getters are written without the use of tensorObject so that they also work for empty nodes
			
			size_t size() const;
			
			size_t degree() const;
			
			void erase();
		};
		
		

        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
            
        ///@brief Dimensions of the external indices, i.e. the dimensions of the tensor represented by the network.
        std::vector<size_t> dimensions;
        
        ///@brief The nodes constituting the network. The order determines the ids of the nodes.
		std::vector<TensorNode> nodes;
            
        ///@brief The open links of the network in order.
        std::vector<Link> externalLinks;
        
        ///@brief A single value representing a constant factor and/or the only entry of an order zero tensor
//         value_t factor;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
            
        /** 
		* @brief Constructs an order zero TensorNetwork.
		* @details The order of an empty TN is zero. The network will contain one node with the single
		* entry zero.
		*/
		explicit TensorNetwork(const misc::NoCast<bool> _addZeroNode = true);
		
	protected:
		
		///@brief Internal indicator to avoid magic false
		static const misc::NoCast<bool> NoZeroNode;
		
	public:
        
        ///@brief Copy Constructor
        implicit TensorNetwork(const TensorNetwork& _cpy);
        
        ///@brief Move Constructor
        implicit TensorNetwork(TensorNetwork&& _mv);
        
		/** 
		* @brief Constructs the trivial TensorNetwork containing the given Tensor as single node.
		* @details The global factor of the TN is set to 1.0 (the tensor may still have a factor on its own).
		*/
        implicit TensorNetwork(const Tensor& _other);
        
		/** 
		* @brief Constructs the trivial TensorNetwork containing the given Tensor as single node.
		* @details The global factor of the TN is set to 1.0 (the tensor may still have a factor on its own).
		*/
        implicit TensorNetwork(Tensor&& _other);
        
		/** 
		* @brief Constructs the trivial TensorNetwork containing the given Tensor as single node.
		* @details The global factor of the TN is set to 1.0 (the tensor may still have a factor on its own).
		* The TN takes the ownership of the pointer.
		*/
        implicit TensorNetwork(std::unique_ptr<Tensor>&&  _tensor);
        
		/** 
		* @brief Constructs the trivial TensorNetwork containing a FullTensor with the given degree.
		* @details All dimensions are set equals one, the global factor of the TN is one and the only entry 
		* of the tensor is zero.
		*/
        implicit TensorNetwork(size_t _degree);
        
        ///@brief Destructor
		virtual ~TensorNetwork() {}
        
        /** 
		* @brief Returns a new copy of the network.
		* @details All dimensions are set equals one, the global factor of the TN is one and the only entry 
		* of the tensor is zero.
		*/
		virtual TensorNetwork* get_copy() const;
            
    private:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		//TODO describtion
		std::vector<Link> init_from_dimension_array();
        
		/** 
		* @brief Contracts all nodes that are not connected to any external links.
		*/
        void contract_unconnected_subnetworks();

    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
            
        /** 
		* @brief Explicit cast to FullTensor
		* @details Contracts the complete network into a single FullTensor
		*/
        explicit operator FullTensor() const;
        
		/** 
		* @brief Explicit cast to SparseTensor
		* @details Contracts the complete network into a single SparseTensor
		*/
        explicit operator SparseTensor() const;
            
        /** 
		* @brief Fully contract the TensorNetwork
		* @details The complete TensorNetwork is contracted. The result can be both full or sparse.
		* @returns a pointer to the resulting single Tensor.
		*/
		std::unique_ptr<Tensor> fully_contracted_tensor() const;
        
        ///@brief TensorNetworks are copy assignable.
		TensorNetwork& operator=(const TensorNetwork &_cpy);
            
        ///@brief TensorNetworks are move assignable.
		TensorNetwork& operator=(TensorNetwork &&_mv);
            
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /** 
		* @brief Read the value at a specific position.
		* @details This allows the efficent calculation of a single entry of the TensorNetwork, by first fixing the external dimensions
		* and then completly contracting the network. Do NOT use this as a manual cast to FullTensor (there is an explicit cast for that).
		* @param _position the position of the entry to be read assuming row-major ordering and a single node.
		* @returns the calculated value (NO reference)
		*/
        value_t operator[](const size_t _position) const;
        
		/** 
		* @brief Read the value at a specific position.
		* @details This allows the efficent calculation of a single entry of the TensorNetwork, by first fixing the external dimensions
		* and then completly contracting the network. Do NOT use this as a manual cast to FullTensor (there is an explicit cast for that).
		* @param _position the position of the entry to be read assuming a single node.
		* @returns the calculated value (NO reference)
		*/
        value_t operator[](const std::vector<size_t>& _positions) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
        /** 
         * @brief Performs the entrywise multiplication with a constant @a _factor.
         * @details Internally this only results in a change in the global factor.
         * @param _factor the factor,
         * @return a reference to this TensorNetwork.
         */
        virtual void operator*=(const value_t _factor);
        
        
        /** 
         * @brief Performs the entrywise divison by a constant @a _divisor.
         * @details Internally this only results in a change in the global factor.
         * @param _divisor the divisor,
         * @return a reference to this TensorNetwork.
         */ 
        virtual void operator/=(const value_t _divisor);
        
		
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		/** 
		 * @brief Indexes the TensorNetwork for read/write use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor(Network).
		 */
		template<typename... args>
		IndexedTensor<TensorNetwork> operator()(args... _args) {
				return IndexedTensor<TensorNetwork>(this, std::vector<Index>({_args...}), false);
		}
		
		/** 
		 * @brief Indexes the TensorNetwork for read only use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor(Network).
		 */
		template<typename... args>
		IndexedTensorReadOnly<TensorNetwork> operator()(args... _args) const {
				return IndexedTensorReadOnly<TensorNetwork>(this, std::vector<Index>({_args...}));
		}
		
		/** 
		 * @brief Indexes the TensorNetwork for read/write use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor(Network).
		 */
		IndexedTensor<TensorNetwork> operator()(const std::vector<Index> & _indices);
        
		/** 
		 * @brief Indexes the TensorNetwork for read/write use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor(Network).
		 */
        IndexedTensor<TensorNetwork> operator()(      std::vector<Index>&& _indices);
            
		/** 
		 * @brief Indexes the TensorNetwork for read only use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor(Network).
		 */
		IndexedTensorReadOnly<TensorNetwork> operator()(const std::vector<Index> & _indices) const;
        
		/** 
		 * @brief Indexes the TensorNetwork for read only use.
		 * @param _args several [indices](@ref Index) determining the desired index order.
		 * @return an internal representation of an IndexedTensor(Network).
		 */
        IndexedTensorReadOnly<TensorNetwork> operator()(      std::vector<Index>&& _indices) const;
            
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
		/// Calculates the contraction between _me and _other and stores the result in _out. Requires that *this is the tensorObjectReadOnly of _me.
        virtual bool specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const;
        
        /// Calculates the sum between _me and _other and stores the result in _out. Requires that *this is the tensorObjectReadOnly of _me.
		virtual bool specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const;
        
        /// Evaluates _other into _me. Requires that *this is the tensorObjectReadOnly of _me.
		virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other);
            
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
        /** 
		 * @brief Gets the degree of the TensorNetwork.
		 * @details The degree is defined as the number of dimensions (i.e. dimensions.size()) 
		 * and is always equal to the number of externalLinks (i.e. externalLinks.size()).
		 * @return the degree.
		 */
        size_t degree() const;
           
		/// @brief reshuffled the nodes according to the given (from, to) map
		void reshuffle_nodes(const std::map<size_t, size_t> &_map);
		
		/// @brief reshuffled the nodes according to the given function
		void reshuffle_nodes(std::function<size_t(size_t)> _f);
		
		/** 
		 * @brief Sanity check for the network.
		 * @details Checks whether all links in the network are set consistently and matching the 
		 * underlying tensor objects. Note that this only checks whether the TensorNetwork is valid
		 * not whether the additional constrains of a specific format are fullfilled. For this purpose
		 * use is_in_expected_format().
		 * @return TRUE if the sanity check passes. If not an exception is throws and the function does not return.
		 */
		bool is_valid_network() const;
        
		/** 
		 * @brief Creates a dataless copy of a subnet.
		 * @details Creates a copy of this TensorNetwork containing the specified nodes,
		 * but do not propagate the data but use the nullptr as data for all nodes
		 * @param _ids the indices of the nodes to be copied. 
		 * @return the new TensorNetwork.
		 */
        TensorNetwork stripped_subnet(std::set<size_t> _ids) const;
        
        
        // TODO describtion
        void swap_external_links(const size_t _i, const size_t _j);
        
        /// shuffles the external links of _lhs according to the indices of the indexedTensors
        /// lhs contains a copy of rhs, thus we have to swap the rhs.indices to resemble those of the lhs
		static void shuffle_indices(std::vector<Index> &_currentIndices, const IndexedTensorWritable<TensorNetwork> &_lhs);
		
        // TODO describtion
		static void add_network_to_network(IndexedTensorWritable<TensorNetwork> & _base, const IndexedTensorReadOnly<TensorNetwork> & _toInsert);
        
        // TODO describtion
		static void trace_out_double_indices(std::vector<Index> &_modifiedIndices, const IndexedTensorWritable<TensorNetwork> & _base);
	protected:
		/** 
		 * @brief Finds the position of a single common edge between two nodes.
		 * @param _posA position of the common edge in the first node.
		 * @param _posB position of the common edge in the second node.
		 * @param _ba Index whose span is set equals the number of dimensions before the common edge in first node.
		 * @param _aa Index whose span is set equals the number of dimensions after the common edge in first node.
		 * @param _bb Index whose span is set equals the number of dimensions before the common edge in second node.
		 * @param _ab Index whose span is set equals the number of dimensions after the common edge in second node.
		 * @param _nodeA The first node.
		 * @param _nodeB The second node.
		 */
		void identify_common_edge(size_t& _posA, size_t& _posB, Index& _ba, Index& _aa, Index& _bb, Index& _ab, const size_t _nodeA, const size_t _nodeB) const;
		
	public:
		/**
		 * @brief Thresholds the rank between two given nodes.
		 * @details The given nodes must be joined by a single edge. Both nodes are contracted and an SVD is calculated to perform the thresholding.
		 * The obtained core is contracted to nodeB, i.e. nodeA remains orthogonal in the used matrification.
		 * @param _nodeA First node that takes part in the rank thresholding. This node remains orthogonalized.
		 * @param _nodeB Second node that takes part in the rank thresholding. This nodes carries the core.
		 * @param _maxRank Maximal allowed rank.
		 * @param _eps Epsilion to be used in the SVD to determine zero singular values.
		 * @param _softThreshold Softthreshold that is to be applied.
		 * @param _preventZero Flag set to prevent the result to become the zero tensor.
		 */
		void round_edge(const size_t _nodeA, const size_t _nodeB, const size_t _maxRank, const double _eps, const double _softThreshold, const bool _preventZero);
		
		/**
		 * @brief Transfers the core from one given node to another.
		 * @details The given nodes must be joined by a single edge. A QR decomposition of the first node is calculated and the core contracted to the second one.
		 * @param _nodeA First node, which remains orthogonalized.
		 * @param _nodeB Second node, which carries the core.
		 * @param _allowRankReduction Flag indicating whether a rank revealing decomposition is to be used which allows the reduction of the rank.
		 */
		void transfer_core(const size_t _nodeA, const size_t _nodeB, const bool _allowRankReduction = true);
		
		
		/** 
		 * @brief Fixes a specific slate in one of the dimensions, effectively reducing the order by one.
		 * @param _dimension the dimension in which the slate shall be fixed, e.g. 0 to fix the first dimensions.
		 * @param _slatePosition the position in the corresponding dimensions that shall be used.
		 */
		void fix_slate(const size_t _dimension, const size_t _slatePosition);
		
		/**
		* contracts the nodes with indices @a _node1 and @a _node2
		* replaces node1 with the contraction and node2 with an degree-0 tensor
		*/
		void contract(size_t _nodeId1, size_t _nodeId2);
		
		/** 
		 * @brief Approximates the cost of contraction two given nodes.
		 * @param _nodeId1 id of the first node.
		 * @param _nodeId2 id of the second node.
		 * @return The approxiamted contraction cost.
		 */
		double contraction_cost(size_t _nodeId1, size_t _nodeId2);
		
		
		/**
		* contracts the nodes with indices included in the set
		* replaces all but one node with degree-0 tensor
		* @returns the id of the contracted tensor
		*/
		size_t contract(std::set<size_t> _ids);
		
		/** 
		 * @brief Calculates the frobenious norm of the TensorNetwork.
		 * @return the frobenious norm of the TensorNetwork.
		 */
		virtual value_t frob_norm() const;
		
		/** 
		 * @brief Sanity check for the TensorNetwork and if applicable for the specific format.
		 * @details Checks whether all links in the network are set consistently and matching the 
		 * underlying tensor objects. This also checks whether the additional constrains of the specific 
		 * format (if any) are fullfilled.
		 * @return TRUE if the sanity check passes. If not an exception is throws and the function does not return.
		 */
		virtual bool is_in_expected_format() const;
		
		/**
		 * @brief Draws a graph representation of the TensorNetwork.
		 * @details The drawing is realized by a system call to "dot" which plots the graph structure.
		 * @param _filename path and name of the file where to save the image.
		 */
		void draw(const std::string& _filename) const;
    };
	
	/** 
	* @brief Calculates the frobenious norm of the given TensorNetwork.
	* @param _network the TensorNetwork of which the frobenious norm shall be calculated.
	* @return the frobenious norm.
	*/
    static _inline_ value_t frob_norm(const TensorNetwork& _network) { return _network.frob_norm(); }
    
    
    std::ostream &operator<<(std::ostream &_out, const TensorNetwork::Link &_rhs);
}
