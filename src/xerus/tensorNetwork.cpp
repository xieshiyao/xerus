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

#include <xerus/tensorNetwork.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/contractionHeuristic.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <algorithm>

namespace xerus {    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    TensorNetwork::TensorNetwork() : factor(0.0) {}
    
    TensorNetwork::TensorNetwork(const TensorNetwork& _cpy) : dimensions(_cpy.dimensions), nodes(_cpy.nodes), externalLinks(_cpy.externalLinks), factor(_cpy.factor) {}
    
    TensorNetwork::TensorNetwork(TensorNetwork&& _mv) : dimensions(std::move(_mv.dimensions)), nodes(std::move(_mv.nodes)), externalLinks(std::move(_mv.externalLinks)), factor(_mv.factor) {}
    
    TensorNetwork::TensorNetwork(const Tensor& _other) : dimensions(_other.dimensions), factor(1.0) {
        nodes.emplace_back(std::shared_ptr<Tensor>(_other.get_copy()), init_from_dimension_array());
    }
    
    TensorNetwork::TensorNetwork(Tensor&& _other) : dimensions(_other.dimensions), factor(1.0) { //NOTE don't use std::move here, because we need _other to be untouched to move it later
        nodes.emplace_back(std::shared_ptr<Tensor>(_other.get_moved_copy()), init_from_dimension_array());
    }
    
    TensorNetwork::TensorNetwork(const std::shared_ptr<Tensor>& _tensor) : dimensions(_tensor->dimensions), factor(1.0) {
        nodes.emplace_back(_tensor, init_from_dimension_array());
    }
    
    /// Constructs the trivial network containing non-specified size-1 FullTensor
    TensorNetwork::TensorNetwork(size_t _degree) : dimensions(std::vector<size_t>(_degree,1)), factor(1.0) {
        nodes.emplace_back(std::shared_ptr<Tensor>(new FullTensor(_degree)), init_from_dimension_array());
    }
    
    TensorNetwork* TensorNetwork::get_copy() const {
		return new TensorNetwork(*this);
	}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    std::vector<TensorNode::Link> TensorNetwork::init_from_dimension_array() {
        std::vector<TensorNode::Link> newLinks;
        for (size_t d=0; d<dimensions.size(); ++d) {
            externalLinks.emplace_back(0, d, dimensions[d], false);
            newLinks.emplace_back(-1, d, dimensions[d], true);
        }
        return newLinks;
    }
    
    bool TensorNetwork::has_factor() const {
        #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wfloat-equal"
            return (factor != 1.0);
        #pragma GCC diagnostic pop
    }
    
    void TensorNetwork::apply_factor() {
        REQUIRE(is_valid_network(), "Cannot apply factor to inconsistent network.");
        if(has_factor()) {
            // Find an unerased node to apply factor to. If there is none, then the factor stays
            for(TensorNode& node : nodes) {
                if(!node.erased) {
                    node.add_factor(factor);
                    factor = 1.0;
                    break;
                }
            }
        }
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    std::unique_ptr<Tensor> TensorNetwork::fully_contracted_tensor() const {
		REQUIRE(is_valid_network(), "Cannot fully contract inconsistent network.");
        std::unique_ptr<Tensor> result;
        
        if (degree() == 0) {
            #ifndef DISABLE_RUNTIME_CHECKS_
                for(const TensorNode& node : nodes) {
                    REQUIRE(node.erased, "Tensor Network with degree zero must not have unerased nodes.");
                }
            #endif
            result.reset(new FullTensor());
            static_cast<FullTensor*>(result.get())->data.get()[0] = factor;
        } else {
            std::set<size_t> all;
            TensorNetwork cpy(*this);
            
            cpy.apply_factor();
            
            for (size_t i=0; i < nodes.size(); ++i) {
                if (!nodes[i].erased) all.insert(i);
            }
            size_t res = cpy.contract(all);
            
            std::vector<Index> externalOrder;
            for(size_t i = 0; i < cpy.nodes[res].neighbors.size(); ++i) { externalOrder.emplace_back(); }
            
            std::vector<Index> internalOrder;
            for(const TensorNode::Link& link: cpy.nodes[res].neighbors) {
                REQUIRE(link.external, "Internal Error");
                internalOrder.emplace_back(externalOrder[link.indexPosition]);
            }
            
            result.reset(cpy.nodes[res].tensorObject->construct_new());
            
            (*result)(externalOrder) = (*cpy.nodes[res].tensorObject)(internalOrder);
        }
        
        return result;
    }
    
    TensorNetwork::operator FullTensor() const {
        std::unique_ptr<Tensor> contractedTensor = fully_contracted_tensor();
        
        if(contractedTensor->is_sparse()) {
            return FullTensor(std::move(*static_cast<SparseTensor*>(fully_contracted_tensor().get())));
        } else {
            return FullTensor(std::move(*static_cast<FullTensor*>(fully_contracted_tensor().get())));
        }
    }
    
    TensorNetwork::operator SparseTensor() const {
        std::unique_ptr<Tensor> contractedTensor = fully_contracted_tensor();
        
        if(contractedTensor->is_sparse()) {
            return SparseTensor(std::move(*static_cast<SparseTensor*>(fully_contracted_tensor().get())));
        } else {
            LOG(error, "Casting a TensorNetwork containing FullTensors to SparseTensor. This is most likely not usefull.");
            return SparseTensor(std::move(*static_cast<FullTensor*>(fully_contracted_tensor().get())));
        }
    }
    
    
    TensorNetwork& TensorNetwork::operator=(const TensorNetwork& _other) {
        REQUIRE(_other.is_valid_network(), "Cannot assign inconsistent network.");
        
        dimensions = _other.dimensions;
        nodes = _other.nodes;
        externalLinks = _other.externalLinks;
        factor = _other.factor;
        return *this;
    }
    
    TensorNetwork& TensorNetwork::operator=(TensorNetwork &&_mv) {
        REQUIRE(_mv.is_valid_network(), "Cannot move inconsistent network.");
        
        dimensions = std::move(_mv.dimensions);
        nodes = std::move(_mv.nodes);
        externalLinks = std::move(_mv.externalLinks);
        factor = _mv.factor;
        
        // Ensure that _mv is in a legal state (ugly because _mv is usually deleted directly after move assignement...)
        _mv.dimensions.clear();
        _mv.nodes.clear();
        _mv.externalLinks.clear();
        return *this;
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    value_t TensorNetwork::operator[](const size_t _position) const {
        std::vector<size_t> positions(degree());
        size_t remains = _position;
        for(size_t i = degree()-1; i > 0; --i) {
            positions[i] = remains%dimensions[i];
            remains /= dimensions[i];
        }
        positions[0] = remains;
        return operator[](positions);
    }
    
    value_t TensorNetwork::operator[](const std::vector<size_t>& _positions) const {
        TensorNetwork partialCopy;
        partialCopy.factor = factor;
        partialCopy.nodes = nodes;
        
//         for(TensorNode& node : partialCopy.nodes) {
//             for(size_t i = 0; i < node.neighbors.size(); ++i) {
//                 if(node.neighbors[i].external) {continue;}
//                 LOG(bla, "neighbor: " << i << " / " << node.neighbors.size());
//                 LOG(vlaues, "Other node " << node.neighbors[i].other << " / " << partialCopy.nodes.size());
//                 LOG(vlaues, "Other Index " << node.neighbors[i].indexPosition << " / " << partialCopy.nodes[node.neighbors[i].other].neighbors.size());
//             }
//         }
        
        // Set all external indices in copy to the fixed values and evaluate the tensorObject accordingly
        for(TensorNode& node : partialCopy.nodes) {
            std::vector<Index> shrinkIndices, baseIndices;
            for(const TensorNode::Link& link : node.neighbors) {
                if(link.external) {
                    // Add fixed index to base
                    baseIndices.emplace_back(_positions[link.indexPosition]);
                } else {
                    // Add real index to both
                    shrinkIndices.emplace_back();
                    baseIndices.push_back(shrinkIndices.back());
                }
            }
            Tensor* tmp = node.tensorObject->construct_new();
            (*tmp)(std::move(shrinkIndices)) = (*node.tensorObject)(std::move(baseIndices));
            node.tensorObject.reset(tmp);
            
            // Remove all external links, because they don't exist anymore
            node.neighbors.erase(std::remove_if(node.neighbors.begin(), node.neighbors.end(), [](const TensorNode::Link& _test){return _test.external;}), node.neighbors.end());
            
            // Adjust the Links
            for(size_t i = 0; i < node.neighbors.size(); ++i) {
                partialCopy.nodes[node.neighbors[i].other].neighbors[node.neighbors[i].indexPosition].indexPosition = i;
            }
        }
        
        // Contract the network (there are not external Links)
        partialCopy.contract_unconnected_subnetworks();
        
        REQUIRE(partialCopy.nodes.empty(), "Internal Error.");
        
        return partialCopy.factor;
    }
    
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    IndexedTensor<TensorNetwork> TensorNetwork::operator()(const std::vector<Index> & _indices) {
        return IndexedTensor<TensorNetwork>(this, _indices, false);
    }
    
    IndexedTensor<TensorNetwork> TensorNetwork::operator()(      std::vector<Index>&& _indices) {
        return IndexedTensor<TensorNetwork>(this, std::move(_indices), false);
    }
        
    IndexedTensorReadOnly<TensorNetwork> TensorNetwork::operator()(const std::vector<Index> & _indices) const {
        return IndexedTensorReadOnly<TensorNetwork>(this, _indices);
    }
    
    IndexedTensorReadOnly<TensorNetwork> TensorNetwork::operator()(      std::vector<Index>&& _indices) const {
        return IndexedTensorReadOnly<TensorNetwork>(this, std::move(_indices));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
	bool TensorNetwork::specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out , const IndexedTensorReadOnly<TensorNetwork> &_me , const IndexedTensorReadOnly<TensorNetwork> &_other ) const {
		return false; // A general tensor Network can't do anything specialized
	}
	
	bool TensorNetwork::specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const {
        return false; // A general tensor Network can't do anything specialized
	}
	
	void TensorNetwork::specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) {
		// If tensorObject is already identical, don't attempt to copy it
		std::vector<Index> currentIndices(_other.get_assigned_indices());
        if (_other.tensorObjectReadOnly != _me.tensorObject) {
            *_me.tensorObject = *_other.tensorObjectReadOnly;
        } 
        TensorNetwork::trace_out_double_indices(currentIndices, _me);
        TensorNetwork::shuffle_indices(currentIndices, _me);
	}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    size_t TensorNetwork::degree() const {
        REQUIRE(externalLinks.size() == dimensions.size(), "invalid network, " << externalLinks.size() << " vs " << dimensions.size());
        return dimensions.size();
    }
    
    /// Eliminates all erased Nodes
    void TensorNetwork::sanitize() {
        std::vector<size_t> idMap(nodes.size(), -1);
        
        // Move nodes in vector
        size_t newId=0, oldId=0;
        for (; oldId < nodes.size(); ++oldId) {
            if (nodes[oldId].erased) { continue; }
            idMap[oldId] = newId;
            if (newId != oldId) { std::swap(nodes[newId], nodes[oldId]); }
            newId++;
        }
        
        // Update links
        nodes.resize(newId);
        for (TensorNode &n : nodes) {
            for (TensorNode::Link &l : n.neighbors) {
                if (!l.external) l.other = idMap[l.other];
            }
        }
        
        // Update external links
        for (TensorNode::Link &l : externalLinks) {
            l.other = idMap[l.other];
        }
    }
    
    void TensorNetwork::reshuffle_nodes(const std::map<size_t, size_t> &_map) {
		reshuffle_nodes([&](size_t i){return _map.at(i);});
	}
	
	void TensorNetwork::reshuffle_nodes(std::function<size_t(size_t)> _f) {
		std::vector<TensorNode> newOrder(nodes.size());
		size_t newSize = 0;
		for (size_t i=0; i<nodes.size(); ++i) {
			if (nodes[i].erased) { continue; }
			size_t newIndex = _f(i);
			newSize = std::max(newSize, newIndex+1);
			TensorNode &newNode = newOrder[ newIndex];
			REQUIRE(newNode.erased, "tried to shuffle two nodes to the same new position " << newIndex << " i= " << i);
			newNode = nodes[i];
			for (TensorNode::Link &l : newNode.neighbors) {
				if (!l.external) { l.other = _f(l.other); }
			}
		}
		nodes = newOrder;
		nodes.resize(newSize);
	}
    
#ifndef DISABLE_RUNTIME_CHECKS_
    /// check whether all links in the network are set consistently and matching the underlying tensor objects
    bool TensorNetwork::is_valid_network() const {
		REQUIRE(std::isfinite(factor), "factor = " << factor);
        
		// Per external link
		for (size_t n=0; n<externalLinks.size(); ++n) {
			const TensorNode::Link &el = externalLinks[n];
			REQUIRE(el.other < nodes.size(), "n=" << n);
            REQUIRE(el.dimension > 0, "n=" << n);
            REQUIRE(el.dimension == dimensions[n], "n=" << n);
			REQUIRE(!el.external, "n=" << n);
			
			const TensorNode &other = nodes[el.other];
			REQUIRE(other.degree() > el.indexPosition, "n=" << n);
			REQUIRE(other.neighbors[el.indexPosition].external, "n=" << n);
			REQUIRE(other.neighbors[el.indexPosition].indexPosition == n, "n=" << n);
			REQUIRE(other.neighbors[el.indexPosition].dimension == el.dimension, "n=" << n << " " << other.neighbors[el.indexPosition].dimension << " vs " << el.dimension);
		}
		
		// per node
		for (size_t n=0; n<nodes.size(); ++n) {
			const TensorNode &currNode = nodes[n];
			REQUIRE(!currNode.erased, "n=" << n);
			if (currNode.tensorObject) {
				REQUIRE(currNode.degree() == currNode.tensorObject->degree(), "n=" << n << " " << currNode.degree() << " vs " << currNode.tensorObject->degree());
			}
			// per neighbor
			for (size_t i=0; i<currNode.neighbors.size(); ++i) {
				const TensorNode::Link &el = currNode.neighbors[i];
				REQUIRE(el.dimension > 0, "n=" << n << " i=" << i);
				if (currNode.tensorObject) {
					REQUIRE(el.dimension==currNode.tensorObject->dimensions[i],  "n=" << n << " i=" << i << " " << el.dimension << " vs " << currNode.tensorObject->dimensions[i]);
				}
				
				if (!el.external) { // externals were already checked
					REQUIRE(el.other < nodes.size(), "n=" << n << " i=" << i << " " << el.other << " vs " << nodes.size());
					const TensorNode &other = nodes[el.other];
					REQUIRE(other.degree() > el.indexPosition, "n=" << n << " i=" << i << " " << other.degree() << " vs " << el.indexPosition);
					REQUIRE(!other.neighbors[el.indexPosition].external, "n=" << n << " i=" << i);
					REQUIRE(other.neighbors[el.indexPosition].other == n, "n=" << n << " i=" << i);
					REQUIRE(other.neighbors[el.indexPosition].indexPosition == i, "n=" << n << " i=" << i);
					REQUIRE(other.neighbors[el.indexPosition].dimension == el.dimension, "n=" << n << " i=" << i << " " << other.neighbors[el.indexPosition].dimension << " vs " << el.dimension);
				}
			}
		}
		
		return true;
	}
#else
	/// no checks are performed with disabled checks... 
    bool TensorNetwork::is_valid_network() const {
		return true;
	}
#endif
	
	bool TensorNetwork::is_in_expected_format() const {
		return is_valid_network();
	}
    
    /// Creates a copy of a subnet that only contains nullptr as data pointers
    TensorNetwork TensorNetwork::stripped_subnet(std::set<size_t> _ids) const {
        TensorNetwork cpy;
        cpy.nodes.resize(nodes.size());
        cpy.dimensions = dimensions;
        cpy.externalLinks = externalLinks;
        for (size_t id : _ids) {
            cpy.nodes[id] = nodes[id].strippped_copy();
            for (size_t i=0; i<cpy.nodes[id].neighbors.size(); ++i) {
                TensorNode::Link &l = cpy.nodes[id].neighbors[i];
                if (!l.external) { // Link was not external before
                    if (!contains(_ids, l.other)) { // ...but is "external" to this subnet
                        l.external = true;
                        l.indexPosition = cpy.externalLinks.size();
                        cpy.dimensions.emplace_back(l.dimension);
                        cpy.externalLinks.emplace_back(id, i, l.dimension, false);
                    } 
                }
            }
        }
        return cpy;
    }
    
    void TensorNetwork::swap_external_links(const size_t _i, const size_t _j) {
        TensorNode::Link &li = externalLinks[_i];
        TensorNode::Link &lj = externalLinks[_j];
        nodes[li.other].neighbors[li.indexPosition].indexPosition = _j;
        nodes[lj.other].neighbors[lj.indexPosition].indexPosition = _i;
        std::swap(externalLinks[_i], externalLinks[_j]);
        std::swap(dimensions[_i], dimensions[_j]);
    }
    

    void TensorNetwork::trace_out_double_indices(std::vector<Index> &_modifiedIndices, const IndexedTensorWritable<TensorNetwork> & _base) {
		TensorNetwork &base = *_base.tensorObject;
        
		REQUIRE(base.is_valid_network(), "Network that is supposed to be traced out is inconsistent.");
        
        #ifndef DISABLE_RUNTIME_CHECKS_
            std::set<Index> contractedIndices;
        #endif
            
        size_t j=0;
        size_t spanSumJ=0;
        while (j<_modifiedIndices.size()) {
            Index &ij = _modifiedIndices[j];
            if (!_base.is_contained_and_open(ij)) {
                REQUIRE(contractedIndices.count(ij) == 0, "Indices must occur at most twice per contraction");
                // trace out this index
                // find the id of the second occurance
                size_t k=j+1;
                size_t spanSumK = spanSumJ+ij.span;
                for (; k<_modifiedIndices.size(); spanSumK+=_modifiedIndices[k].span, ++k) {
                    if (ij == _modifiedIndices[k]) {
                        break;
                    }
                }
                REQUIRE(k<_modifiedIndices.size(), "ie " << _modifiedIndices << " k " << k << " j " << j << " ij " << ij << " contained and open? " << _base.is_contained_and_open(ij));
                
                for (size_t n=0; n<ij.span; ++n) {
                    TensorNode::Link &link1 = base.externalLinks[spanSumJ];
                    TensorNode::Link &link2 = base.externalLinks[spanSumK-n];
                    base.nodes[link1.other].neighbors[link1.indexPosition] = link2;
                    base.nodes[link2.other].neighbors[link2.indexPosition] = link1;
                    
                    // remove external links and dimensions from network
                    base.externalLinks.erase(base.externalLinks.begin()+(long)(spanSumK-n)); // note that spanSumK > spanSumJ
                    base.externalLinks.erase(base.externalLinks.begin()+(long)(spanSumJ));
                    base.dimensions.erase(base.dimensions.begin()+(long)(spanSumK-n)); // note that spanSumK > spanSumJ
                    base.dimensions.erase(base.dimensions.begin()+(long)(spanSumJ));
                    
                    for(size_t i = spanSumJ; i< spanSumK-n-1; ++i) {
                        base.nodes[base.externalLinks[i].other].neighbors[base.externalLinks[i].indexPosition].indexPosition -= 1;
                    }
                    for(size_t i = spanSumK-n-1; i < base.externalLinks.size(); ++i) {
                        base.nodes[base.externalLinks[i].other].neighbors[base.externalLinks[i].indexPosition].indexPosition -= 2;
                    }
                }
                #ifndef DISABLE_RUNTIME_CHECKS_
                    contractedIndices.insert(ij);
                #endif
                // remove index from indexed tensor and mark as contracted
                _modifiedIndices.erase(_modifiedIndices.begin() + (long)k); // note that k>j
                _modifiedIndices.erase(_modifiedIndices.begin() + (long)j);
            } else {
                spanSumJ += ij.span;
                j+=1;
            }
        }
        
        base.contract_unconnected_subnetworks();
		
		REQUIRE(base.is_valid_network(), "Network was broken in the process of tracing out double indices.");
    }
    
    // Fully contracts TN parts that are not connected to external links
    void TensorNetwork::contract_unconnected_subnetworks() {
        std::vector<bool> seen(nodes.size(), false);
        std::vector<size_t> expansionStack;
        
        // Starting at every external link...
        for (TensorNode::Link &el : externalLinks) {
            if(!seen[el.other]) {
                seen[el.other] = true;
                expansionStack.push_back(el.other);
            }
        }
        
        
        // ...traverse the connected nodes in a depth-first manner.
        while (!expansionStack.empty()) {
            const size_t curr = expansionStack.back();
            expansionStack.pop_back();
            
            // Add unseen neighbors
            for (const TensorNode::Link &n : nodes[curr].neighbors) {
                if ( !n.external && !seen[n.other] ) {
                    seen[n.other] = true;
                    expansionStack.push_back(n.other);
                }
            }
        }
        
        // Construct set of all unseen nodes...
        std::set<size_t> toContract;
        for (size_t i=0; i < nodes.size(); ++i) {
            if (!seen[i] && !nodes[i].erased) {
                toContract.insert(i);
            }
        }
        
        // ...and contract them
        if (!toContract.empty()) {
            size_t remaining = contract(toContract);
            
            // remove contracted degree-0 tensor
            REQUIRE(nodes[remaining].degree() == 0, "Internal Error.");
            factor *= (*nodes[remaining].tensorObject.get())[0];
            REQUIRE(nodes[remaining].neighbors.empty(), "Internal Error.");
            nodes[remaining].erased = true;
        }
        
        // Remove all erased nodes
        sanitize();
    }


    //TODO testcase A(i,j)*B(k,k) of TTtensors
    /**
     * contracts the nodes with indices @a _node1 and @a _node2
     * replaces node1 with the contraction and node2 with an degree-0 tensor
     */
    void TensorNetwork::contract(size_t _nodeId1, size_t _nodeId2) {
        TensorNode &node1 = nodes[_nodeId1];
        TensorNode &node2 = nodes[_nodeId2];
        
        REQUIRE(!node1.erased, "It appears node1 was already contracted?");
        REQUIRE(!node2.erased, "It appears node2 was already contracted?");
        REQUIRE(externalLinks.size() == degree(), "Internal Error: " << externalLinks.size() << " != " << degree());
        
        LOG(TNContract, "contraction of " << _nodeId1 << " and " << _nodeId2 << " size " << nodes.size());
        
        std::vector<Index> i1, i2, ri;
        std::vector<TensorNode::Link> newLinks;
        for (size_t d=0; d<node1.degree(); ++d) {
            // self-link?
            if (node1.neighbors[d].links(_nodeId1) && d > node1.neighbors[d].indexPosition) {
                i1.push_back(i1[node1.neighbors[d].indexPosition]);
            } else {
                i1.emplace_back();
            }
            
            if (!node1.neighbors[d].links(_nodeId2) && !node1.neighbors[d].links(_nodeId1)) {
                ri.push_back(i1.back());
                newLinks.push_back(node1.neighbors[d]);
            }
        }
        for (size_t d=0; d<node2.degree(); ++d) {
            if (node2.neighbors[d].links(_nodeId1)) {
                REQUIRE(i1.size()>node2.neighbors[d].indexPosition, "Internal Error.");
                REQUIRE(node1.neighbors[node2.neighbors[d].indexPosition].links(_nodeId2), "Internal Error."); // consistency check
                i2.push_back(i1[node2.neighbors[d].indexPosition]); // common index
            } else if (node2.neighbors[d].links(_nodeId2) && d > node2.neighbors[d].indexPosition) {
                i2.push_back(i2[node2.neighbors[d].indexPosition]);
            } else {
                i2.emplace_back(); // new (open) index
                ri.push_back(i2.back());
                newLinks.push_back(node2.neighbors[d]);
            }
        }
        
        // Construct new node
        std::shared_ptr<Tensor> newTensor;
        if (node1.tensorObject) {
            REQUIRE(node2.tensorObject, "Internal Error.");
            
            IndexedTensorMoveable<Tensor> result(xerus::contract((*node1.tensorObject)(i1), (*node2.tensorObject)(i2)));
            newTensor.reset(result.tensorObjectReadOnly->construct_new());
            (*newTensor.get())(ri) = result;
        } else {
            REQUIRE(!node2.tensorObject, "Internal Error.");
        } 
        TensorNode newNode(newTensor, newLinks);
        
        // replace node1 and node2
        nodes[_nodeId1] = newNode;
        nodes[_nodeId2].erase();
        
        // fix indices of other nodes // note that the indices that were previously part of node1 might also have changed
        for (size_t d=0; d<newLinks.size(); ++d) {
            TensorNode::Link *nn;
            if (newLinks[d].external) {
                nn = &externalLinks[newLinks[d].indexPosition];
            } else {
                nn = &nodes[newLinks[d].other].neighbors[newLinks[d].indexPosition];
            }
            REQUIRE(nn->other == _nodeId2 || nn->other == _nodeId1, "Internal Error."); // consistency check
            nn->other = _nodeId1;
            nn->indexPosition = d;
        }
    }

    double TensorNetwork::contraction_cost(size_t _nodeId1, size_t _nodeId2) {
        TensorNode &node1 = nodes[_nodeId1];
        TensorNode &node2 = nodes[_nodeId2];
        
        REQUIRE(!node1.erased, "It appears node1 was already contracted?");
        REQUIRE(!node2.erased, "It appears node2 was already contracted?");
        
        if (_nodeId1 == _nodeId2) {
            return (double)node1.size(); // cost of a trace
        }
        
        // assume cost of strassen with mxr * rxn = m*n*r (which is incorrect)
        //TODO add correct calculation for sparse matrices
        size_t m=1,n=1,r=1;
        for (size_t d=0; d<node1.degree(); ++d) {
            if (node1.neighbors[d].links(_nodeId2)) {
                r *= node1.neighbors[d].dimension;
            } else {
                m *= node1.neighbors[d].dimension;
            }
        }
        for (size_t d=0; d<node2.degree(); ++d) {
            if (!node2.neighbors[d].links(_nodeId1)) {
                n *= node2.neighbors[d].dimension;
            }
        }
        return (double) (m*n*r);
    }

    void TensorNetwork::add_network_to_network(IndexedTensorWritable<TensorNetwork> & _base, const IndexedTensorReadOnly<TensorNetwork> & _toInsert) {
        const std::vector<Index> toInsertIndices = _toInsert.get_assigned_indices();
        
        // Any factor has to be propageded
        _base.tensorObject->factor *= _toInsert.tensorObjectReadOnly->factor;
        
        // TODO after trace ensure connectedness (to external indices)
        TensorNetwork &base = *_base.tensorObject;
        const TensorNetwork &toInsert = *_toInsert.tensorObjectReadOnly;
        
        size_t firstNew = base.nodes.size();
        size_t firstNewExternal = base.externalLinks.size();
        
        // enlarge externalLinks vector by _rhs links
        for (const TensorNode::Link &l : toInsert.externalLinks) {
            base.externalLinks.emplace_back(l);
            base.externalLinks.back().other += firstNew;
        }
        
        // Merge indices
        _base.indices.insert(_base.indices.end(), toInsertIndices.begin(), toInsertIndices.end());
        base.dimensions.insert(base.dimensions.end(), toInsert.dimensions.begin(), toInsert.dimensions.end());
        
        #ifndef DISABLE_RUNTIME_CHECKS_
        for (const Index &idx : _base.indices) {
            REQUIRE(count(_base.indices, idx) < 3, "Internal Error.");
        }
        #endif
        
        // Add network (treating all external links as new external links)
        for (size_t i=0; i<toInsert.nodes.size(); ++i) {
            base.nodes.emplace_back(toInsert.nodes[i]);
            for(TensorNode::Link &l : base.nodes.back().neighbors) {
                if (!l.external) { // Link inside the added network
                    l.other += firstNew;
                } else { // External link
                    l.indexPosition += firstNewExternal;
                }
            }
        }
        
        // Find traces (former contractions have become traces due to the joining)
        trace_out_double_indices(_base.indices, _base);
		
		REQUIRE(_base.tensorObject->is_valid_network(), "ie");
    }


    //TODO vorwaerts fressende heuristik (ein knoten + nachbarn und nachbars-nachbarn -> exakte dreier berechnung)
    //TODO contract with std::function
    size_t TensorNetwork::contract(std::set<size_t> _ids) {
        LOG(TNContract, "contraction of " << _ids.size() << " nodes called");
        
        // TODO this and all the heuristics still assume that all tensors are dense (as there are no sparse tensors at the time of this writing)
        if (_ids.size() == 0) { return -1; }
        
        // trace out all single-node traces
        for (size_t id : _ids) {
            REQUIRE(!nodes[id].erased, "tried to contract erased node");
            bool traceNeeded=false;
            for (const TensorNode::Link &l : nodes[id].neighbors) {
                if (l.links(id)) {
                    traceNeeded = true;
                    break;
                }
            }
            if (traceNeeded) {
                LOG(TNContract, "single-node trace of id " << id);
                std::vector<Index> idxIn, idxOut;
                std::vector<TensorNode::Link> newLinks;
                for (size_t i=0; i<nodes[id].neighbors.size(); ++i) {
                    const TensorNode::Link &l = nodes[id].neighbors[i];
                    if (!l.links(id)) {
                        idxIn.emplace_back();
                        idxOut.emplace_back(idxIn.back());
                        newLinks.emplace_back(l);
                    } else if (l.indexPosition > i) {
                        idxIn.emplace_back();
                    } else {
                        idxIn.emplace_back(idxIn[l.indexPosition]);
                    }
                }
                // Perform the trace
                if (nodes[id].tensorObject) {
                    std::shared_ptr<Tensor> newTensor(nodes[id].tensorObject->construct_new());
                    (*newTensor)(idxOut) = (*nodes[id].tensorObject)(idxIn);
                    nodes[id].tensorObject = std::move(newTensor);
                    nodes[id].neighbors = std::move(newLinks);
                }
            }
        }
        
        if (_ids.size() == 1) { return *_ids.begin(); }
        if (_ids.size() == 2) {
            auto secItr = _ids.begin(); ++secItr;
            contract(*_ids.begin(), *secItr);
            return *_ids.begin();
        }
        
        if (_ids.size() == 3) {
            auto idItr = _ids.begin();
            size_t a=*idItr; TensorNode &na = nodes[a]; ++idItr;
            size_t b=*idItr; TensorNode &nb = nodes[b]; ++idItr;
            size_t c=*idItr; TensorNode &nc = nodes[c];
            float sa=1, sb=1, sc=1; // sizes devided by the link dimensions between a,b,c
            float sab=1, sbc=1, sac=1; // link dimensions
            for (size_t d=0; d<na.degree(); ++d) {
                if (na.neighbors[d].links(b)) {
                    sab *= (float)na.neighbors[d].dimension;
                } else if (na.neighbors[d].links(c)) {
                    sac *= (float)na.neighbors[d].dimension;
                } else {
                    sa *= (float)na.neighbors[d].dimension;
                }
            }
            for (size_t d=0; d<nb.degree(); ++d) {
                if (nb.neighbors[d].links(c)) {
                    sbc *= (float) nb.neighbors[d].dimension;
                } else if (!nb.neighbors[d].links(a)) {
                    sb *= (float) nb.neighbors[d].dimension;
                }
            }
            for (size_t d=0; d<nc.degree(); ++d) {
//                 size_t other = nc.neighbors[d].other;
                if (!nc.neighbors[d].links(a) && !nc.neighbors[d].links(b)) {
                    sc *= (float)nc.neighbors[d].dimension;
                }
            }
            // cost of contraction a-b first etc.
            float costAB = sa*sb*sac*sbc*(sab+sc); // (sa*sac)*sab*(sb*sbc) + sa*sb*sac*sbc*sc;
            float costAC = sa*sc*sab*sbc*(sac+sb); 
            float costBC = sb*sc*sab*sac*(sbc+sa);
            if (costAB < costAC && costAB < costBC) {
                LOG(TNContract, "contraction of ab first " << sa << " " << sb << " " << sc << " " << sab << " " << sbc << " " << sac);
                contract(a,b); contract(a,c); return a;
            } else if (costAC < costBC) {
                LOG(TNContract, "contraction of ac first " << sa << " " << sb << " " << sc << " " << sab << " " << sbc << " " << sac);
                contract(a,c); contract(a,b); return a;
            } else {
                LOG(TNContract, "contraction of bc first " << sa << " " << sb << " " << sc << " " << sab << " " << sbc << " " << sac);
                contract(b,c); contract(a,b); return a;
            }
        }
        
        
        TensorNetwork stripped = stripped_subnet(_ids); 
        float bestScore=1e32f;
        std::vector<std::pair<size_t, size_t>> *bestOrder=nullptr;
        
        // ask heuristics
        for (ContractionHeuristic &c : *ContractionHeuristic::list) {
            if (c.rescore(stripped) < bestScore) {
                bestScore = c.score;
                bestOrder = &c.contractions; //TODO not threadsafe
            }
            LOG(TNContract, "heuristic '" << c.name << "' contracts with cost " << c.score);
            //if (bestScore < 256*256*256) break;
        }
        
        REQUIRE(bestScore < 1e32f && bestOrder, "Internal Error.");
        // perform contractions
        for (const std::pair<size_t,size_t> &c : *bestOrder) {
            contract(c.first, c.second);
        }
        
        // note: no sanitization as eg. TTStacks require the indices not to change after calling this function
        return bestOrder->back().first;
    }

    value_t TensorNetwork::frob_norm() const {
        Index i;
        FullTensor res(0);
        res() = (*this)(i&0) * (*this)(i&0);
        return res.data.get()[0];
    }

}
