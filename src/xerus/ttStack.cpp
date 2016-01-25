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
 * @brief Implementation of the TTStack classes.
 */

#include <xerus/ttStack.h>
#include <xerus/basic.h>
#include <xerus/misc/check.h>

#include <xerus/index.h>
#include <xerus/tensor.h>
#include <xerus/ttNetwork.h>
 

namespace xerus {
	namespace internal {
		template<bool isOperator>
		TTStack<isOperator>::TTStack(const bool _canno, const size_t _corePos) : cannonicalization_required(_canno), futureCorePosition(_corePos) {}
		
		
		template<bool isOperator>
		TensorNetwork* TTStack<isOperator>::get_copy() const {
			return new TTStack(*this);
		}
		
		
		template<bool isOperator>
		void TTStack<isOperator>::operator*=(const value_t _factor) {
			REQUIRE(nodes.size() > 0, "There must not be a TTNetwork without any node");
			
			if(cannonicalization_required) {
				*nodes[futureCorePosition+1].tensorObject *= _factor;
			} else if(degree() > 0) {
				*nodes[1].tensorObject *= _factor;
			} else {
				*nodes[0].tensorObject *= _factor;
			}
		}
		
		
		template<bool isOperator>
		void TTStack<isOperator>::operator/=(const value_t _divisor) {
			operator*=(1/_divisor);
		}
		
		
		template<bool isOperator>
		void TTStack<isOperator>::contract_stack(IndexedTensorWritable<TensorNetwork>&& _me) {
			_me.tensorObject->require_valid_network();
			
			const size_t numComponents = _me.degree()/N;
			const size_t numNodes = _me.degree()/N+2;
			
			std::set<size_t> toContract;
			for (size_t currentNode = 0; currentNode < numNodes; ++currentNode) {
				toContract.clear();
				for (size_t i = currentNode; i < _me.tensorObject->nodes.size(); i+=numNodes) {
					toContract.insert(i);
				}
				_me.tensorObject->contract(toContract);
			}
			
			// all are contracted, reshuffle them to be in the correct order
			// after contraction the nodes will have one of the ids: node, node+numNodes, node+2*numNodes,... (as those were part of the contraction)
			// so modulus gives the correct wanted id
			_me.tensorObject->reshuffle_nodes([numNodes](size_t i){return i%(numNodes);});
			
			REQUIRE(_me.tensorObject->nodes.size() == numNodes, "Internal Error.");
			_me.tensorObject->require_valid_network();
			
			// Reset to new external links
			_me.tensorObject->externalLinks.clear();
			for(size_t i = 0; i < numComponents; ++i) {
				_me.tensorObject->externalLinks.emplace_back(i+1, 1, _me.tensorObject->dimensions[i], false);
			}
			if(N == 2) {
				for(size_t i = 0; i < numComponents; ++i) {
					_me.tensorObject->externalLinks.emplace_back(i+1, 2, _me.tensorObject->dimensions[numComponents+i], false);
				}
			}
			
			// Ensure right amount and order of links
			Index ext[N];
			size_t lastRank, externalDim[N], newRank=1;
			std::vector<Index> lastIndices, lastRight;
			std::vector<Index> oldIndices, newRight; // newLeft == lastRight
			std::vector<Index> newIndices;
			std::vector<size_t> newDimensions;
			_me.tensorObject->nodes.front().neighbors = std::vector<TensorNetwork::Link>({TensorNetwork::Link(1,0,1,false)});
			_me.tensorObject->nodes.front().tensorObject->reinterpret_dimensions({1});
			_me.tensorObject->nodes.back().neighbors = std::vector<TensorNetwork::Link>({TensorNetwork::Link(numComponents,N+1,1,false)});
			_me.tensorObject->nodes.back().tensorObject->reinterpret_dimensions({1});
			for (size_t i=0; i<numComponents; ++i) {
				lastIndices = std::move(oldIndices); oldIndices.clear();
				lastRight = std::move(newRight); newRight.clear();
				lastRank = newRank; newRank=1;
				TensorNode &n = _me.tensorObject->nodes[i+1];
				for (TensorNetwork::Link &l : n.neighbors) {
					if (l.external) {
						size_t externalNumber = 0;
						if (isOperator) {
							externalNumber = l.indexPosition>=numComponents?1:0;
						}
						oldIndices.push_back(ext[externalNumber]);
						externalDim[externalNumber] = l.dimension;
					} else if (l.links(i)) {
						if (i==0) {
							lastRight.emplace_back();
							oldIndices.push_back(lastRight.back());
						} else {
							REQUIRE(lastIndices.size() > l.indexPosition, "ie " << i << " " << lastIndices.size() << " " << l.indexPosition);
							oldIndices.push_back(lastIndices[l.indexPosition]);
						}
					} else if (l.links(i+2)) {
						oldIndices.emplace_back();
						newRight.push_back(oldIndices.back());
						newRank *= l.dimension;
					} else  {
						LOG(fatal, "Internal Error.");
					}
				}
				newIndices = std::move(lastRight);
				newIndices.insert(newIndices.end(), ext, ext+N);
				newIndices.insert(newIndices.end(), newRight.begin(), newRight.end());
				
				(*n.tensorObject)(newIndices) = (*n.tensorObject)(oldIndices);
				
				newDimensions.clear();
				n.neighbors.clear();
				n.neighbors.emplace_back(i, i==0?0:N+1,lastRank, false);
				newDimensions.push_back(lastRank);
				for (size_t j=0; j<N; ++j) {
					REQUIRE(_me.tensorObject->dimensions[i+j*numComponents] == externalDim[j], "Internal Error.");
					n.neighbors.emplace_back(0,i+j*numComponents,externalDim[j], true);
					newDimensions.push_back(externalDim[j]);
				}
				n.neighbors.emplace_back(i+2,0,newRank, false);
				newDimensions.push_back(newRank);
				n.tensorObject->reinterpret_dimensions(newDimensions);
			}
			
			// NOTE core position according to information in TTStack is set in evaluation
		}
		
		/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
		template<bool isOperator>
		void TTStack<isOperator>::specialized_evaluation(IndexedTensorWritable<TensorNetwork>&& _me _unused_ , IndexedTensorReadOnly<TensorNetwork>&& _other _unused_) {
			LOG(fatal, "TTStack not supported as a storing type");
		}
		
		
		template<bool isOperator>
		bool TTStack<isOperator>::specialized_contraction(std::unique_ptr<IndexedTensorMoveable<TensorNetwork>>& _out, IndexedTensorReadOnly<TensorNetwork>&& _me, IndexedTensorReadOnly<TensorNetwork>&& _other) const {
			return TTNetwork<isOperator>::specialized_contraction_f(_out, std::move(_me), std::move(_other));
		}
		
		
		template<bool isOperator>
		bool TTStack<isOperator>::specialized_sum(std::unique_ptr<IndexedTensorMoveable<TensorNetwork>>& _out, IndexedTensorReadOnly<TensorNetwork>&& _me, IndexedTensorReadOnly<TensorNetwork>&& _other) const {
			return TTNetwork<isOperator>::specialized_sum_f(_out, std::move(_me), std::move(_other));
		}
		
		
		template<bool isOperator>
		value_t TTStack<isOperator>::frob_norm() const {
			const Index i;
			TTNetwork<isOperator> tmp;
			tmp(i&0) = IndexedTensorMoveable<TensorNetwork>(this->get_copy(), {i&0});
			return tmp.frob_norm();
		}
		
		
		// Explicit instantiation of the two template parameters that will be implemented in the xerus library
		template class TTStack<false>;
		template class TTStack<true>;
	}
}
