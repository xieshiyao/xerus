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
 * @brief Implementation of the TensorNetwork class.
 */

#include <xerus/tensorNetwork.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/contractionHeuristic.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/indexedTensor_tensor_factorisations.h>
#include <xerus/indexedTensorList.h>
#include <xerus/misc/stringUtilities.h>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <thread> // TODO define guard?

namespace xerus {
	const misc::NoCast<bool> TensorNetwork::NoZeroNode(false);
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    TensorNetwork::TensorNetwork(const misc::NoCast<bool> _addZeroNode) {
		if(_addZeroNode) {
			nodes.emplace_back(TensorNode(std::unique_ptr<Tensor>(new FullTensor())));
		}
	}
    
    TensorNetwork::TensorNetwork(const TensorNetwork& _cpy) : dimensions(_cpy.dimensions), nodes(_cpy.nodes), externalLinks(_cpy.externalLinks) {}
    
    TensorNetwork::TensorNetwork(TensorNetwork&& _mv) : dimensions(std::move(_mv.dimensions)), nodes(std::move(_mv.nodes)), externalLinks(std::move(_mv.externalLinks)){
		_mv = TensorNetwork();
	}
    
    TensorNetwork::TensorNetwork(const Tensor& _other) : dimensions(_other.dimensions) {
        nodes.emplace_back(std::unique_ptr<Tensor>(_other.get_copy()), init_from_dimension_array());
    }
    
    TensorNetwork::TensorNetwork(Tensor&& _other) : dimensions(_other.dimensions) { //NOTE don't use std::move here, because we need _other to be untouched to move it later
        nodes.emplace_back(std::unique_ptr<Tensor>(_other.get_moved_copy()), init_from_dimension_array());
    }
    
    TensorNetwork::TensorNetwork( std::unique_ptr<Tensor>&& _tensor) : dimensions(_tensor->dimensions) {
        nodes.emplace_back(std::move(_tensor), init_from_dimension_array());
    }
    
    /// Constructs the trivial network containing non-specified size-1 FullTensor
    TensorNetwork::TensorNetwork(size_t _degree) : dimensions(std::vector<size_t>(_degree,1)) {
        nodes.emplace_back(std::unique_ptr<Tensor>(new FullTensor(_degree)), init_from_dimension_array());
    }
    
    TensorNetwork* TensorNetwork::get_copy() const {
		return new TensorNetwork(*this);
	}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    std::vector<TensorNetwork::Link> TensorNetwork::init_from_dimension_array() {
        std::vector<TensorNetwork::Link> newLinks;
        for (size_t d=0; d<dimensions.size(); ++d) {
            externalLinks.emplace_back(0, d, dimensions[d], false);
            newLinks.emplace_back(-1, d, dimensions[d], true);
        }
        return newLinks;
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    std::unique_ptr<Tensor> TensorNetwork::fully_contracted_tensor() const {
		REQUIRE(is_valid_network(), "Cannot fully contract inconsistent network.");
        std::unique_ptr<Tensor> result;
        
		std::set<size_t> all;
		TensorNetwork cpy(*this);
		
		for (size_t i=0; i < nodes.size(); ++i) {
			if (!nodes[i].erased) all.insert(i);
		}
		size_t res = cpy.contract(all);
		
		std::vector<Index> externalOrder;
		for(size_t i = 0; i < cpy.nodes[res].neighbors.size(); ++i) { externalOrder.emplace_back(); }
		
		std::vector<Index> internalOrder;
		for(const TensorNetwork::Link& link: cpy.nodes[res].neighbors) {
			REQUIRE(link.external, "Internal Error");
			internalOrder.emplace_back(externalOrder[link.indexPosition]);
		}
		
		result.reset(cpy.nodes[res].tensorObject->construct_new());
		
		(*result)(externalOrder) = (*cpy.nodes[res].tensorObject)(internalOrder);
        
        return result;
    }
    
    TensorNetwork::operator FullTensor() const {
        std::unique_ptr<Tensor> contractedTensor = fully_contracted_tensor();
        
        if(contractedTensor->is_sparse()) {
            return FullTensor(std::move(*static_cast<SparseTensor*>(fully_contracted_tensor().get())));
        } else {
            return std::move(*static_cast<FullTensor*>(fully_contracted_tensor().get()));
        }
    }
    
    TensorNetwork::operator SparseTensor() const {
        std::unique_ptr<Tensor> contractedTensor = fully_contracted_tensor();
        
        if(contractedTensor->is_sparse()) {
            return std::move(*static_cast<SparseTensor*>(fully_contracted_tensor().get()));
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
        return *this;
    }
    
    TensorNetwork& TensorNetwork::operator=(TensorNetwork &&_mv) {
        REQUIRE(_mv.is_valid_network(), "Cannot move inconsistent network.");
        
        dimensions = std::move(_mv.dimensions);
        nodes = std::move(_mv.nodes);
        externalLinks = std::move(_mv.externalLinks);
        
        // Ensure that _mv is in a legal state (ugly because _mv is usually deleted directly after move assignement...)
        _mv.dimensions.clear();
        _mv.nodes.clear();
        _mv.externalLinks.clear();
        return *this;
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    value_t TensorNetwork::operator[](const size_t _position) const {
		REQUIRE(is_valid_network(), "Invalid Network");
		REQUIRE(nodes.size() > 0, "There must always be at least one node");
		if (degree() == 0) {
			REQUIRE(_position == 0, "Tried to access non-existing entry of TN");
			REQUIRE(nodes.size() == 1, "Internal Error");
			return (*nodes[0].tensorObject)[0];
		}
        std::vector<size_t> positions(degree());
        size_t remains = _position;
        for(size_t i = degree(); i > 1; --i) {
            positions[i-1] = remains%dimensions[i-1];
            remains /= dimensions[i-1];
        }
        positions[0] = remains;
        return operator[](positions);
    }
    
    value_t TensorNetwork::operator[](const std::vector<size_t>& _positions) const {
		REQUIRE(is_valid_network(), "Invalid Network");
        TensorNetwork partialCopy;
        partialCopy.nodes = nodes;
        
        // Set all external indices in copy to the fixed values and evaluate the tensorObject accordingly
        for(TensorNode& node : partialCopy.nodes) {
			// Fix slates in external links
			size_t killedDimensions = 0;
            for(size_t i = 0; i < node.neighbors.size(); ++i) {
                if(node.neighbors[i].external) {
                    node.tensorObject->fix_slate(i-killedDimensions, _positions[node.neighbors[i].indexPosition]);
					killedDimensions++;
                }
            }
            
            // Remove all external links, because they don't exist anymore
            node.neighbors.erase(std::remove_if(node.neighbors.begin(), node.neighbors.end(), [](const TensorNetwork::Link& _test){return _test.external;}), node.neighbors.end());
            
            // Adjust the Links
            for(size_t i = 0; i < node.neighbors.size(); ++i) {
                partialCopy.nodes[node.neighbors[i].other].neighbors[node.neighbors[i].indexPosition].indexPosition = i;
            }
        }
        
        // Contract the complete network (there are not external Links)
        partialCopy.contract_unconnected_subnetworks();
        
        REQUIRE(partialCopy.nodes.size() == 1, "Internal Error.");
        
        return (*partialCopy.nodes[0].tensorObject)[0];
    }
    
    
	void TensorNetwork::measure(std::vector<SinglePointMeasurment>& _measurments) const {
		std::vector<TensorNetwork> stack;
		stack.emplace_back(*this);
		stack.back().reduce_representation();
		
// 		std::map<size_t, size_t> histo;
		
		// Sort measurements
		std::sort(_measurments.begin(), _measurments.end(), SinglePointMeasurment::Comparator(degree()-1));
		
		bool firstTime = true;
		std::vector<size_t> previousPosition;
		for(SinglePointMeasurment& measurment : _measurments) {
			size_t rebuildIndex = ~0ul;
			if(firstTime) {
				rebuildIndex = 0;
				firstTime = false;
			} else {
				// Find the maximal recyclable stack position
				for(size_t i = 0; i < degree(); ++i) {
					if(previousPosition[i] != measurment.positions[i]) {
						rebuildIndex = i;
						break;
					}
				}
			}
			REQUIRE(rebuildIndex != ~0ul, "There were two identical measurements? pos: " << previousPosition);
			previousPosition = measurment.positions;
			
// 			histo[rebuildIndex] += 1;
			
			// Trash stack that is not needed anymore
			while(stack.size() > rebuildIndex+1) { stack.pop_back(); }
			
			// Rebuild stack
			for(size_t i = rebuildIndex; i < degree(); ++i) {
				stack.emplace_back(stack.back());
				stack.back().fix_slate(0, measurment.positions[i]); 
				stack.back().reduce_representation();
			}
			REQUIRE(stack.back().degree() == 0, "Internal Error");
			REQUIRE(stack.size() == degree()+1, "Internal Error");
/*			
			size_t bla = 0;
			for(const TensorNetwork& tn : stack) {
				tn.draw(std::string("stack_lvl") + misc::to_string(bla++));
			}*/
			
			measurment.value = stack.back()[0];
		}
		
// 		size_t sum=0;
// 		size_t integrated = 0;
// 		for (auto &a: histo) {
// 			integrated += a.second;
// 			LOG(histo, a.first << "\t" << a.second << "\t" << integrated);
// 			sum += (degree()-a.first) * a.second;
// 		}
// 		LOG(histo, "total: " << sum << " mean: " << double(sum)/double(_measurments.size()));
	}
    

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	void TensorNetwork::operator*=(const value_t _factor) {
		REQUIRE(nodes.size() > 0, "There must not be a TTNetwork without any node");
		*nodes[0].tensorObject *= _factor;
	}
	
	void TensorNetwork::operator/=(const value_t _divisor) {
		REQUIRE(nodes.size() > 0, "There must not be a TTNetwork without any node");
		*nodes[0].tensorObject /= _divisor;
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
			for (TensorNetwork::Link &l : newNode.neighbors) {
				if (!l.external) { l.other = _f(l.other); }
			}
		}
		nodes = newOrder;
		nodes.resize(newSize);
		for (size_t i=0; i<externalLinks.size(); ++i) {
			externalLinks[i].other = _f(externalLinks[i].other);
		}
		is_valid_network();
	}
    
#ifndef DISABLE_RUNTIME_CHECKS_
    /// check whether all links in the network are set consistently and matching the underlying tensor objects
    bool TensorNetwork::is_valid_network() const {
        REQUIRE(externalLinks.size() == dimensions.size(), "externalLinks.size() != dimensions.size()");
		REQUIRE(nodes.size() > 0, "There must always be at least one node!");
		
		// Per external link
		for (size_t n = 0; n < externalLinks.size(); ++n) {
			const TensorNetwork::Link &el = externalLinks[n];
			REQUIRE(el.other < nodes.size(), "n=" << n);
            REQUIRE(el.dimension > 0, "n=" << n);
            REQUIRE(el.dimension == dimensions[n], "n=" << n);
			REQUIRE(!el.external, "n=" << n);
			
			const TensorNode &other = nodes[el.other];
			REQUIRE(other.degree() > el.indexPosition, "n=" << n);
			REQUIRE(other.neighbors[el.indexPosition].external, "n=" << n);
			REQUIRE(other.neighbors[el.indexPosition].indexPosition == n, "n=" << n << " We have " << other.neighbors[el.indexPosition].indexPosition << " != " << n << " and el.other =" << el.other << " and el.indexPosition = " << el.indexPosition );
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
				const TensorNetwork::Link &el = currNode.neighbors[i];
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
	TensorNetwork TensorNetwork::stripped_subnet(std::function<bool(size_t)> _idF) const {
		TensorNetwork cpy(NoZeroNode);
		cpy.nodes.resize(nodes.size());
		cpy.dimensions = dimensions;
		cpy.externalLinks = externalLinks;
		for (size_t id =0; id<nodes.size(); ++id) {
			if (!_idF(id)) continue;
			cpy.nodes[id] = nodes[id].strippped_copy();
			for (size_t i=0; i<cpy.nodes[id].neighbors.size(); ++i) {
				TensorNetwork::Link &l = cpy.nodes[id].neighbors[i];
				if (!l.external) { // Link was not external before
					if (!_idF(l.other)) { // ...but is "external" to this subnet
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
	
	TensorNetwork TensorNetwork::stripped_subnet(std::set<size_t> _ids) const {
		return stripped_subnet([&](size_t _id){ return _ids.count(_id)>0; });
	}
    
    
    void TensorNetwork::swap_external_links(const size_t _i, const size_t _j) {
        TensorNetwork::Link &li = externalLinks[_i];
        TensorNetwork::Link &lj = externalLinks[_j];
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
                    TensorNetwork::Link &link1 = base.externalLinks[spanSumJ];
                    TensorNetwork::Link &link2 = base.externalLinks[spanSumK-n];
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
    void TensorNetwork::identify_common_edge(size_t& _posA, size_t& _posB, Index& _ba, Index& _aa, Index& _bb, Index& _ab, const size_t _nodeA, const size_t _nodeB) const {
		// Find common Index in nodeA
		bool foundCommon = false;
		for(size_t i = 0; i < nodes[_nodeA].neighbors.size(); ++i) {
			if(nodes[_nodeA].neighbors[i].other == _nodeB) {
				REQUIRE(!foundCommon, "TN round does not work if the two nodes share more than one link.");
				foundCommon = true;
				_posA = i;
			}
		}
		REQUIRE(foundCommon, "TN round does not work if the two nodes share no link.");
		
		// Find common Index in nodeB
		foundCommon = false;
		for(size_t i = 0; i < nodes[_nodeB].neighbors.size(); ++i) {
			if(nodes[_nodeB].neighbors[i].other == _nodeA) {
				REQUIRE(!foundCommon, "TN round does not work if the two nodes share more than one link.");
				foundCommon = true;
				_posB = i;
			}
		}
		REQUIRE(foundCommon, "TN round does not work if the two nodes share no link.");
		
		// Set the spans of the indices
		_ba.span = _posA;
		_aa.span = nodes[_nodeA].degree() - _posA - 1;
		_bb.span = _posB;
		_ab.span = nodes[_nodeB].degree() - _posB - 1;
	}
    
    void TensorNetwork::round_edge(const size_t _nodeA, const size_t _nodeB, const size_t _maxRank, const double _eps, const double _softThreshold, const bool _preventZero) {
		size_t posA, posB;
		Index ba, aa, bb, ab, c1, c2, k, l;
		identify_common_edge(posA, posB, ba, aa, bb, ab, _nodeA, _nodeB);
		
		Tensor& tensorA = (*nodes[_nodeA].tensorObject);
		Tensor& tensorB = (*nodes[_nodeB].tensorObject);
		
		FullTensor X; // TODO Sparse
		SparseTensor S;

		// TODO eventually use only one QC if only one is usefull.
		
		// Check whether prior QR makes sense
		if (tensorA.size > misc::sqr(tensorA.dimensions[posA]) 
		 || tensorB.size > misc::sqr(tensorB.dimensions[posB])) {
			// Calculate the cores
			FullTensor coreA, coreB;
			(tensorA(ba, c1, aa), coreA(c1, k)) = QC(tensorA(ba, k, aa));
			(tensorB(bb, c2, ab), coreB(k, c2)) = QC(tensorB(bb, k, ab)); // TODO use CQ when available
			
			// Contract the cores
			X(c1, c2) = coreA(c1, k)*coreB(k, c2);
			
			// Calculate SVD
			(coreA(c1, k), S(k,l), coreB(l, c2)) = SVD(X(c1, c2), _maxRank, _eps, _softThreshold, _preventZero);
			
			// Contract diagnonal matrix to coreB
			coreB(l, c2) = S(l, k) * coreB(k, c2);
			
			// Contract the "cores" back to their sides
			tensorA(ba, k, aa) = tensorA(ba, c1, aa) * coreA(c1, k);
			tensorB(bb, k, ab) = coreB(k, c2) * tensorB(bb, c2, ab);
			
		} else {
			// Contract the two
			X(ba, aa, bb, ab) = (*nodes[_nodeA].tensorObject)(ba, c1, aa) * (*nodes[_nodeB].tensorObject)(bb, c1, ab);
			
			// Calculate SVD
			(((*nodes[_nodeA].tensorObject)(ba, c1, aa)), S(c1, c2), ((*nodes[_nodeB].tensorObject)(bb, c2, ab))) = SVD(X(ba, aa, bb, ab), _maxRank, _eps, _softThreshold, _preventZero);
			
			// Contract diagnonal matrix to NodeB
			(*nodes[_nodeB].tensorObject)(bb, c1, ab) = S(c1, c2) * ((*nodes[_nodeB].tensorObject)(bb, c2, ab));
		}
		
		// Set the new dimension in the nodes
		nodes[_nodeA].neighbors[posA].dimension = S.dimensions[0];
		nodes[_nodeB].neighbors[posB].dimension = S.dimensions[0];
	}
	
	
	void TensorNetwork::transfer_core(const size_t _nodeA, const size_t _nodeB, const bool _allowRankReduction) {
		size_t posA, posB;
		Index ba, aa, bb, ab, c1, c2;
		identify_common_edge(posA, posB, ba, aa, bb, ab, _nodeA, _nodeB);
		
		// Calculate QR
		FullTensor X;
		if(_allowRankReduction) {
			((*nodes[_nodeA].tensorObject)(ba, c2, aa), X(c2, c1)) = QC((*nodes[_nodeA].tensorObject)(ba, c1, aa));
		} else {
			((*nodes[_nodeA].tensorObject)(ba, c2, aa), X(c2, c1)) = QR((*nodes[_nodeA].tensorObject)(ba, c1, aa));
		}
		
		// Contract diagnonal matrix to NodeB
		(*nodes[_nodeB].tensorObject)(bb, c1, ab) = X(c1, c2) * ((*nodes[_nodeB].tensorObject)(bb, c2, ab));
		
		// Set the new dimension in the nodes
		nodes[_nodeA].neighbors[posA].dimension = X.dimensions[0];
		nodes[_nodeB].neighbors[posB].dimension = X.dimensions[0];
	}
	
	
	void TensorNetwork::fix_slate(const size_t _dimension, const size_t _slatePosition) {
		is_valid_network();
		const size_t extNode = externalLinks[_dimension].other;
		const size_t extNodeIndexPos = externalLinks[_dimension].indexPosition;
		
		for(size_t i = _dimension+1; i < dimensions.size(); ++i) {
			nodes[externalLinks[i].other].neighbors[externalLinks[i].indexPosition].indexPosition--;
		}
		
		externalLinks.erase(externalLinks.begin()+(long)_dimension);
		dimensions.erase(dimensions.begin()+(long)_dimension);
		
		for(size_t i = extNodeIndexPos+1; i < nodes[extNode].neighbors.size(); ++i) {
			const Link& link = nodes[extNode].neighbors[i];
			if(link.external) {
				externalLinks[link.indexPosition].indexPosition--; 
			} else {
				nodes[link.other].neighbors[link.indexPosition].indexPosition--;
			}
		}
		
		nodes[extNode].tensorObject->fix_slate(extNodeIndexPos, _slatePosition);
		nodes[extNode].neighbors.erase(nodes[extNode].neighbors.begin() + (long)extNodeIndexPos);
		
		contract_unconnected_subnetworks();
		is_valid_network();
	}
    
    
    void TensorNetwork::contract_unconnected_subnetworks() {
		REQUIRE(is_valid_network(), "Invalid TensorNetwork");
        std::vector<bool> seen(nodes.size(), false);
        std::vector<size_t> expansionStack;
		expansionStack.reserve(nodes.size());
        
        // Starting at every external link...
        for (TensorNetwork::Link &el : externalLinks) {
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
            for (const TensorNetwork::Link &n : nodes[curr].neighbors) {
                if ( !n.external && !seen[n.other] ) {
                    seen[n.other] = true;
                    expansionStack.push_back(n.other);
                }
            }
        }
        
        // Construct set of all unseen nodes...
        std::set<size_t> toContract;
        for (size_t i=0; i < nodes.size(); ++i) {
            if (!seen[i]) {
                toContract.insert(i);
            }
        }
        
        const bool keepFinalNode = degree() == 0;
        
        // ...and contract them
        if (!toContract.empty()) {
            size_t remaining = contract(toContract);
            
            // remove contracted degree-0 tensor
            REQUIRE(nodes[remaining].degree() == 0, "Internal Error.");
            REQUIRE(nodes[remaining].neighbors.empty(), "Internal Error.");
			if(!keepFinalNode) {
				for(size_t i =0; i < nodes.size(); ++i) {
					if(i != remaining && !nodes[i].erased) {
						*nodes[i].tensorObject *= (*nodes[remaining].tensorObject.get())[0];
						break;
					}
					REQUIRE(i < nodes.size()-1, "Internal Error.");
				}
				nodes[remaining].erased = true;
			}
        }
        
        // Remove all erased nodes
        std::vector<size_t> idMap(nodes.size(), ~0ul);
        
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
            for (TensorNetwork::Link &l : n.neighbors) {
                if (!l.external) l.other = idMap[l.other];
            }
        }
        
        // Update external links
        for (TensorNetwork::Link &l : externalLinks) {
            l.other = idMap[l.other];
        }
        
        REQUIRE(nodes.size() > 0, "Internal error");
		REQUIRE(!keepFinalNode || nodes.size() == 1, "internal error!");
    }

    
    void TensorNetwork::reduce_representation() {
		REQUIRE(is_valid_network(), "");
		TensorNetwork strippedNet = stripped_subnet();
		std::vector<std::set<size_t>> contractions(strippedNet.nodes.size());
		for (size_t id1=0; id1 < strippedNet.nodes.size(); ++id1) {
			TensorNode &currNode = strippedNet.nodes[id1];
			if (currNode.erased) {
				continue;
			}
			for (Link &l : currNode.neighbors) {
				if (l.external) continue;
				size_t r=1;
				for (Link &l2 : currNode.neighbors) {
					if (l2.other == l.other) {
						r *= l2.dimension;
					}
				}
				if (r*r >= currNode.size() || r*r >= strippedNet.nodes[l.other].size()) {
					if (contractions[id1].empty()) {
						contractions[id1].insert(id1);
					}
					if (contractions[l.other].empty()) {
						contractions[id1].insert(l.other);
					} else {
						contractions[id1].insert(contractions[l.other].begin(), contractions[l.other].end());
						contractions[l.other].clear();
					}
					strippedNet.contract(id1, l.other);
					id1 -= 1; // check the same node again in the next iteration
					break; // for-each iterator is broken, so we have to break
				}
			}
		}
		
		// perform the collected contractions from above
		for (std::set<size_t> &ids : contractions) {
			if (ids.size() > 1) {
				contract(ids);
			}
		}
		
		sanitize();
		is_valid_network();
	}
	
	
	void TensorNetwork::sanitize() {
		size_t idCount=0;
		std::map<size_t, size_t> idMap;
		for (size_t i=0; i<nodes.size(); ++i) {
			if (!nodes[i].erased) {
				idMap[i] = idCount++;
			}
		}
		reshuffle_nodes(idMap);
	}
    

    //TODO testcase A(i,j)*B(k,k) of TTtensors
    /**
     * contracts the nodes with indices @a _node1 and @a _node2
     * replaces node1 with the contraction and node2 with an degree-0 tensor
     */
    void TensorNetwork::contract(size_t _nodeId1, size_t _nodeId2) {
        TensorNode &node1 = nodes[_nodeId1];
        TensorNode &node2 = nodes[_nodeId2];
        
        REQUIRE(!node1.erased, "It appears node1 = " << _nodeId1 << "  was already contracted?");
        REQUIRE(!node2.erased, "It appears node2 = " << _nodeId2 << "  was already contracted?");
        REQUIRE(externalLinks.size() == degree(), "Internal Error: " << externalLinks.size() << " != " << degree());
        
        LOG(TNContract, "contraction of " << _nodeId1 << " and " << _nodeId2 << " size " << nodes.size());
        
        std::vector<Index> i1, i2, ri;
		i1.reserve(node1.degree());
		i2.reserve(node2.degree());
		ri.reserve(node1.degree() + node2.degree());
		
        std::vector<TensorNetwork::Link> newLinks;
		newLinks.reserve(node1.degree() + node2.degree());
		
        for (size_t d=0; d<node1.degree(); ++d) {
            // self-link?
            if (node1.neighbors[d].links(_nodeId1) && d > node1.neighbors[d].indexPosition) {
                i1.emplace_back(i1[node1.neighbors[d].indexPosition]);
            } else {
                i1.emplace_back();
            }
            
            if (!node1.neighbors[d].links(_nodeId2) && !node1.neighbors[d].links(_nodeId1)) {
                ri.emplace_back(i1.back());
                newLinks.emplace_back(node1.neighbors[d]);
            }
        }
        
        for (size_t d = 0; d < node2.degree(); ++d) {
            if (node2.neighbors[d].links(_nodeId1)) {
                REQUIRE(i1.size() > node2.neighbors[d].indexPosition, "Internal Error.");
                REQUIRE(node1.neighbors[node2.neighbors[d].indexPosition].links(_nodeId2), "Internal Error."); // consistency check
                i2.emplace_back(i1[node2.neighbors[d].indexPosition]); // common index
            } else if (node2.neighbors[d].links(_nodeId2) && d > node2.neighbors[d].indexPosition) {
                i2.emplace_back(i2[node2.neighbors[d].indexPosition]);
            } else {
                i2.emplace_back(); // new (open) index
                ri.emplace_back(i2.back());
                newLinks.emplace_back(node2.neighbors[d]);
            }
        }
        
        // Construct new node
        std::unique_ptr<Tensor> newTensor;
        if (node1.tensorObject) {
            REQUIRE(node2.tensorObject, "Internal Error.");
            
            IndexedTensorMoveable<Tensor> result(xerus::contract((*node1.tensorObject)(i1), (*node2.tensorObject)(i2)));
            newTensor.reset(result.tensorObjectReadOnly->construct_new());
            (*newTensor.get())(ri) = result;
        } else {
            REQUIRE(!node2.tensorObject, "Internal Error.");
        } 
        TensorNode newNode(std::move(newTensor), newLinks);
        
        // replace node1 and node2
        nodes[_nodeId1] = newNode;
        nodes[_nodeId2].erase();
        
        // fix indices of other nodes // note that the indices that were previously part of node1 might also have changed
        for (size_t d=0; d<newLinks.size(); ++d) {
            TensorNetwork::Link *nn;
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
        
        REQUIRE(!node1.erased, "It appears node1 = " << _nodeId1 << " was already contracted?");
        REQUIRE(!node2.erased, "It appears node2 = " << _nodeId2 << " was already contracted?");
        
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
        
        // TODO after trace ensure connectedness (to external indices)
        TensorNetwork &base = *_base.tensorObject;
        const TensorNetwork &toInsert = *_toInsert.tensorObjectReadOnly;
        
        size_t firstNew = base.nodes.size();
        size_t firstNewExternal = base.externalLinks.size();
        
        // enlarge externalLinks vector by _rhs links
        for (const TensorNetwork::Link &l : toInsert.externalLinks) {
            base.externalLinks.emplace_back(l);
            base.externalLinks.back().other += firstNew;
        }
        
        // Merge indices
        _base.indices.insert(_base.indices.end(), toInsertIndices.begin(), toInsertIndices.end());
        base.dimensions.insert(base.dimensions.end(), toInsert.dimensions.begin(), toInsert.dimensions.end());
        
        #ifndef DISABLE_RUNTIME_CHECKS_
        for (const Index &idx : _base.indices) {
            REQUIRE(misc::count(_base.indices, idx) < 3, "Internal Error.");
        }
        #endif
        
        // Add network (treating all external links as new external links)
        for (size_t i=0; i<toInsert.nodes.size(); ++i) {
            base.nodes.emplace_back(toInsert.nodes[i]);
            for(TensorNetwork::Link &l : base.nodes.back().neighbors) {
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


    //TODO contract with std::function
    size_t TensorNetwork::contract(std::set<size_t> _ids) {
        LOG(TNContract, "contraction of " << _ids.size() << " nodes called");
        
        // TODO this and all the heuristics still assume that all tensors are dense (as there are no sparse tensors at the time of this writing)
        if (_ids.size() == 0) { return ~0ul; }
        
        // trace out all single-node traces
        for (size_t id : _ids) {
            bool traceNeeded=false;
            for (const TensorNetwork::Link &l : nodes[id].neighbors) {
                if (l.links(id)) {
                    traceNeeded = true;
                    break;
                }
            }
            if (traceNeeded) {
                LOG(TNContract, "single-node trace of id " << id);
                std::vector<Index> idxIn, idxOut;
                std::vector<TensorNetwork::Link> newLinks;
                for (size_t i=0; i<nodes[id].neighbors.size(); ++i) {
                    const TensorNetwork::Link &l = nodes[id].neighbors[i];
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
                    std::unique_ptr<Tensor> newTensor(nodes[id].tensorObject->construct_new());
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
            double sa=1, sb=1, sc=1; // sizes devided by the link dimensions between a,b,c
            double sab=1, sbc=1, sac=1; // link dimensions
            for (size_t d=0; d<na.degree(); ++d) {
                if (na.neighbors[d].links(b)) {
                    sab *= (double)na.neighbors[d].dimension;
                } else if (na.neighbors[d].links(c)) {
                    sac *= (double)na.neighbors[d].dimension;
                } else {
                    sa *= (double)na.neighbors[d].dimension;
                }
            }
            for (size_t d=0; d<nb.degree(); ++d) {
                if (nb.neighbors[d].links(c)) {
                    sbc *= (double) nb.neighbors[d].dimension;
                } else if (!nb.neighbors[d].links(a)) {
                    sb *= (double) nb.neighbors[d].dimension;
                }
            }
            for (size_t d=0; d<nc.degree(); ++d) {
//                 size_t other = nc.neighbors[d].other;
                if (!nc.neighbors[d].links(a) && !nc.neighbors[d].links(b)) {
                    sc *= (double)nc.neighbors[d].dimension;
                }
            }
            // cost of contraction a-b first etc.
            double costAB = sa*sb*sac*sbc*(sab+sc); // (sa*sac)*sab*(sb*sbc) + sa*sb*sac*sbc*sc;
            double costAC = sa*sc*sab*sbc*(sac+sb); 
            double costBC = sb*sc*sab*sac*(sbc+sa);
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
        
        
        TensorNetwork strippedNetwork = stripped_subnet(_ids); 
        double bestCost=std::numeric_limits<double>::max();
        std::vector<std::pair<size_t, size_t>> bestOrder;
        
        // ask heuristics
        for (const internal::ContractionHeuristic &c : internal::contractionHeuristics) {
            c(bestCost, bestOrder, strippedNetwork);
        }
        
        REQUIRE(bestCost < std::numeric_limits<double>::max() && !bestOrder.empty(), "Internal Error.");
        for (const std::pair<size_t,size_t> &c : bestOrder) {
            contract(c.first, c.second);
        }
        
        // note: no sanitization as eg. TTStacks require the indices not to change after calling this function
        return bestOrder.back().first;
    }

    value_t TensorNetwork::frob_norm() const {
        Index i;
        FullTensor res;
        res() = (*this)(i&0) * (*this)(i&0);
        return std::sqrt(res.data.get()[0]);
    }
    
    
	void TensorNetwork::draw(const std::string& _filename) const {
		static std::mutex fileNameMutex; // TODO define guards?
		fileNameMutex.lock();
		char* const tmpFileName = tempnam(nullptr, "xerusToDot");
		std::fstream graphLayout(tmpFileName, std::fstream::out);
		fileNameMutex.unlock();
		
		graphLayout << "graph G {" << std::endl;
		graphLayout << "graph [mclimit=1000, maxiter=1000, overlap = false, splines = true]" << std::endl;
		 
		for(size_t i = 0; i < nodes.size(); ++i) {
			// Create the Nodes
			if(nodes[i].erased) {
				graphLayout << "\tN"<<i<<" [label=\"N"<<i<<"\", shape=circle, fixedsize=shape, height=0.45];" << std::endl;
			} else {
				graphLayout << "\tN"<<i<<" [label=\"";
				for(size_t k=0; k+1 < nodes[i].degree(); ++k) {
					if(nodes[i].degree()/2 == k) {
						if(nodes[i].degree()%2 == 0) {
						graphLayout << "<i"<<k<<"> "<<i<<"| ";
						} else {
							graphLayout << "<i"<<k<<"> N"<<i<<"| ";
						}
					} else if(nodes[i].degree()%2 == 0 && nodes[i].degree()/2 == k+1) {
						graphLayout << "<i"<<k<<"> N| "; 
					} else {
						graphLayout << "<i"<<k<<"> | ";
					}
				}
				if(nodes[i].degree() <= 2) {
					graphLayout << "<i"<<nodes[i].degree()-1<<"> N"<<i<<"\", shape=record, fixedsize=shape, height=0.45, style=\"rounded,filled\"];" << std::endl;
				} else {
					graphLayout << "<i"<<nodes[i].degree()-1<<">\", shape=record, fixedsize=shape, height=0.45, style=\"rounded,filled\"];" << std::endl;
				}
				
				// Add all links to nodes with smaller index and externals
				for(size_t j = 0; j < nodes[i].neighbors.size(); ++j) {
					if(nodes[i].neighbors[j].external) {
						graphLayout << "\t"<<nodes[i].neighbors[j].indexPosition<<" [shape=diamond, fixedsize=shape, height=0.38, width=0.38, style=filled];" << std::endl;
						graphLayout << "\tN"<<i<<":i"<<j<<" -- " << nodes[i].neighbors[j].indexPosition << " [len=1, label=\""<<nodes[i].neighbors[j].dimension<<"\"];" << std::endl;
					} else if(nodes[i].neighbors[j].other < i) {
						graphLayout << "\tN"<<i<<":i"<<j<<" -- " << "N"<<nodes[i].neighbors[j].other << ":i"<< nodes[i].neighbors[j].indexPosition<<" [label=\""<<nodes[i].neighbors[j].dimension<<"\"];" << std::endl;
					}
				}
			}
		}
		graphLayout << "}" << std::endl;
		graphLayout.close();
		misc::exec(std::string("dot -Tsvg ") + tmpFileName + " > " + _filename+".svg");
		
		// Delete tmp File
		remove(tmpFileName);
		
		// Free stupid C memory!
		free(tmpFileName);
	}
}
