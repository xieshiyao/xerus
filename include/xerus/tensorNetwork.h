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

#pragma once

#include "tensorNode.h"
#include "indexedTensor.h"
#include <map>
#include <set>


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
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Member variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
            
        ///@brief Dimensions of the external indices, i.e. the dimensions of the tensor represented by the network.
        std::vector<size_t> dimensions;
        
        ///@brief The nodes constituting the network. The order determines the ids of the nodes.
		std::vector<TensorNode> nodes;
            
        ///@brief The open links of the network in order.
        std::vector<TensorNode::Link> externalLinks;
        
        ///@brief A single value representing a constant factor and/or the only entry of an order zero tensor
        value_t factor;
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
            
        /** 
		* @brief Constructs an empty TensorNetwork.
		* @details The order of an empty TN is zero. The network will contain no nodes 
		* and the global factor, which also determines the only entry, is set to zero.
		*/
		explicit TensorNetwork();
        
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
		std::vector<TensorNode::Link> init_from_dimension_array();
        
		/** 
		* @brief Checks whether there is a non-trivial global scaling factor.
		* @details That is it checks whether factor != 1.0
		* @return True is there is a non trivial factor, FALSE otherwise.
		*/
        bool has_factor() const;
        
		/** 
		* @brief Applies the factor to the network.
		* @details If there is a non trivial global scaling factor, it is apllied to one node.
		* In case of special decompositions, e.g. TTTensor the canonicalisation is conservered.
		*/
        virtual void apply_factor();
        
		/** 
		* @brief Contracts all nodes that are not connected to any external links.
		*/
        void contract_unconnected_subnetworks();

    public:
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
            
        /// Allows explicit casts to FullTensor
        explicit operator FullTensor() const;
        
        /// Allows explicit casts to SparseTensor
        explicit operator SparseTensor() const;
            
        /// Fully contracts the network to a single tensor and returns it as a unique_ptr. Result can be both full or sparse.
        std::unique_ptr<Tensor> fully_contracted_tensor() const;
        
        /// TensorNetworks are copy assignable.
		TensorNetwork &operator=(const TensorNetwork &_cpy);
            
        /// TensorNetworks are move assignable.
		TensorNetwork &operator=(TensorNetwork &&_mv);
            
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        /// Allows read access to the entry at _position, assuming row-major ordering and a single node.
        value_t operator[](const size_t _position) const;
        
        /// Allows access to the entry at _position.
        value_t operator[](const std::vector<size_t>& _positions) const;
        
        
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
		template<typename... args>
		IndexedTensor<TensorNetwork> operator()(args... _args) {
				return IndexedTensor<TensorNetwork>(this, std::vector<Index>({_args...}), false);
		}
		
		template<typename... args>
		IndexedTensorReadOnly<TensorNetwork> operator()(args... _args) const {
				return IndexedTensorReadOnly<TensorNetwork>(this, std::vector<Index>({_args...}));
		}
		
		IndexedTensor<TensorNetwork> operator()(const std::vector<Index> & _indices);
        
        IndexedTensor<TensorNetwork> operator()(      std::vector<Index>&& _indices);
            
		IndexedTensorReadOnly<TensorNetwork> operator()(const std::vector<Index> & _indices) const;
        
        IndexedTensorReadOnly<TensorNetwork> operator()(      std::vector<Index>&& _indices) const;
            
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
        /// Calculates the contraction between _me and _other and stores the result in _out. Requires that *this is the tensorObjectReadOnly of _me.
        virtual bool specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const;
        
        /// Calculates the sum between _me and _other and stores the result in _out. Requires that *this is the tensorObjectReadOnly of _me.
		virtual bool specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const;
        
        /// Evaluates _other into _me. Requires that *this is the tensorObjectReadOnly of _me.
		virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other);
            
            
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
        /// Returns the degree of the tensor network , i.e. the number of externalLinks.
        size_t degree() const;
        
        /// Eleminates all erased Nodes
        void sanitize();
            
		/// reshuffled the nodes according to the given (from, to) map
		void reshuffle_nodes(const std::map<size_t, size_t> &_map);
		
		/// reshuffled the nodes according to the given function
		void reshuffle_nodes(std::function<size_t(size_t)> _f);
		
		/// check whether all links in the network are set consistently and matching the underlying tensor objects
		bool is_valid_network() const;
        
        /// Creates a copy of a subnet that only contains nullptr as data pointers
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
	
		/**
		* contracts the nodes with indices @a _node1 and @a _node2
		* replaces node1 with the contraction and node2 with an degree-0 tensor
		*/
		void contract(size_t _nodeId1, size_t _nodeId2);
		
        // TODO describtion
		double contraction_cost(size_t _nodeId1, size_t _nodeId2);
		
		
		/**
		* contracts the nodes with indices included in the set
		* replaces all but one node with degree-0 tensor
		* @returns the id of the contracted tensor
		*/
		size_t contract(std::set<size_t> _ids);
		
        /// Calculates the frobenious norm of the tensor represented by the tensor network.
		virtual value_t frob_norm() const;
		
		/**
		 * checks whether the given TensorNetwork adheres to the format definition it expects according to its other virtual overloads
		 * eg. calls is_valid_tt for tt tensors and is_valid_network by default
		 */
		virtual bool is_in_expected_format() const;
    };
}
