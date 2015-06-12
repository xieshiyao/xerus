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

#include <xerus/ttNetwork.h>
#include <xerus/tensorNetwork.h>
#include <xerus/basic.h>
#include <xerus/index.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/indexedTensorList.h>
#include <xerus/indexedTensor_TN_operators.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/indexedTensor_tensor_factorisations.h>
#include <xerus/misc/blasLapackWrapper.h>
#include <xerus/misc/selectedFunctions.h>
#include <xerus/misc/check.h>

namespace xerus {
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	template<bool isOperator>
	TTNetwork<isOperator>::TTNetwork() : TensorNetwork(), cannonicalized(false) {}
	
	
	template<bool isOperator>
	TTNetwork<isOperator>::TTNetwork(const size_t _degree) : TensorNetwork(), cannonicalized(true), corePosition(0) {
		REQUIRE(_degree%N==0, "illegal degree for TTOperator");
		const size_t numComponents = _degree/N;
		factor = 0.0;
		
		if (numComponents == 0) {
			return;
		}
		
		dimensions = std::vector<size_t>(_degree, 1);
		
		// ExternalLinks
		for (size_t i = 1; i <= numComponents; ++i) {
			externalLinks.emplace_back(i, 1, 1, false);
		}
		if (N == 2) {
			for (size_t i = 1; i <= numComponents; ++i) {
				externalLinks.emplace_back(i, 2, 1, false);
			}
		}
		
		REQUIRE(externalLinks.size() == _degree, "Internal Error.");
		
		std::vector<TensorNode::Link> neighbors;
		
		neighbors.emplace_back(1,0,1,false);
		
		nodes.emplace_back(
			std::unique_ptr<Tensor>(new FullTensor({1},[](){return 1.0;})), 
			std::move(neighbors)
		);
		
		for (size_t i=0; i<numComponents; ++i) {
			neighbors.clear();
			neighbors.emplace_back(i, i==0?0:N+1, 1, false);
			for (size_t j=0; j< N; ++j) { 
				neighbors.emplace_back(-1, i+j*numComponents, 1, true);
			}
			neighbors.emplace_back(i+2, 0, 1, false);
			
			nodes.emplace_back(
				std::unique_ptr<Tensor>(new FullTensor(neighbors.size())), 
				std::move(neighbors)
			);
		}
		
		neighbors.clear();
		neighbors.emplace_back(numComponents, N+1, 1, false);
		nodes.emplace_back(
			std::unique_ptr<Tensor>(new FullTensor({1},[](){return 1.0;})), 
			std::move(neighbors)
		);
		
		REQUIRE(is_valid_tt(), "ie");
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>::TTNetwork(const Tensor& _tensor, const double _eps): TTNetwork(_tensor.degree()) {
		REQUIRE(_eps < 1, "_eps must be smaller than one. " << _eps << " was given.");
		REQUIRE(_tensor.degree()%N==0, "Number of indicis must be even for TTOperator");
		
		dimensions = _tensor.dimensions;
		
		if (_tensor.degree() == 0) { 
			factor = _tensor[0];
			return; 
		}
		factor = 1.0;
		
		const size_t numComponents = degree()/N;
		
		// Needed variables
		std::unique_ptr<Tensor> nxtTensor;
		std::unique_ptr<value_t[]> currentU, currentS;
		std::shared_ptr<value_t> workingData, currentVt;
		size_t leftDim=1, remainingDim=_tensor.size, maxRank, newRank=1, oldRank=1;
		
		// If we want a TTOperator we need to reshuffle the indices first, otherwise we want to copy the data because Lapack wants to destroy it
		if (!isOperator) {
			if(_tensor.is_sparse()) {
				FullTensor tmpTensor(static_cast<const SparseTensor&>(_tensor));
				workingData = std::move(tmpTensor.data);
			} else {
				workingData.reset(new value_t[_tensor.size], internal::array_deleter_vt);
				misc::array_copy(workingData.get(), static_cast<const FullTensor&>(_tensor).data.get(), _tensor.size);
			}
		} else {
			FullTensor tmpTensor(degree());
			std::vector<Index> presentIndices, newIndices;
			for(size_t i = 0; i < degree(); ++i) { presentIndices.emplace_back(); }
			for(size_t i = 0; i < numComponents; ++i) {
				newIndices.emplace_back(presentIndices[i]);
				newIndices.emplace_back(presentIndices[i+numComponents]);
			}
			tmpTensor(newIndices) = _tensor(presentIndices);
			workingData = std::move(tmpTensor.data);
		}
		
		for(size_t position = 0; position < numComponents-1; ++position) {
			// Determine the dimensions of the next matrification
			leftDim = oldRank*dimensions[position];
			if (isOperator) { leftDim *= dimensions[position+numComponents]; }
			remainingDim /= dimensions[position];
			if (isOperator) { remainingDim /= dimensions[position+numComponents]; }
			maxRank = std::min(leftDim, remainingDim);
			
			// create temporary space for the results
			currentU.reset(new value_t[leftDim*maxRank]);
			currentS.reset(new value_t[maxRank]);
			currentVt.reset(new value_t[maxRank*remainingDim], &internal::array_deleter_vt);
			
			blasWrapper::svd_destructive(currentU.get(), currentS.get(), currentVt.get(), workingData.get(), leftDim, remainingDim);
			
			// Determine the rank, keeping all singular values that are large enough
			newRank = maxRank;
			while (currentS[newRank-1] < _eps*currentS[0]) {
				newRank-=1;
			}
			
			// Create a FullTensor for U
			std::vector<size_t> constructionDim;
			constructionDim.emplace_back(oldRank);
			constructionDim.emplace_back(dimensions[position]);
			if (isOperator) { constructionDim.emplace_back(dimensions[position+numComponents]); }
			constructionDim.emplace_back(newRank);
			if (newRank == maxRank) {
				nxtTensor.reset(new FullTensor(std::move(constructionDim), std::move(currentU)) );
			} else {
				nxtTensor.reset(new FullTensor(std::move(constructionDim), DONT_SET_ZERO()) );
				for (size_t i = 0; i < leftDim; ++i) {
					misc::array_copy(static_cast<FullTensor*>(nxtTensor.get())->data.get()+i*newRank, currentU.get()+i*maxRank, newRank);
				}
			}
			
			// Update component tensor to U
			set_component(position, std::move(nxtTensor));
			
			// Calclate S*Vt by scaling the rows of Vt
			for (size_t row = 0; row < newRank; ++row) {
				misc::array_scale(currentVt.get()+row*remainingDim, currentS[row], remainingDim);
			}
			
			workingData = std::move(currentVt);
			oldRank = newRank;
		}
		
		// Create FullTensor for Vt
		if (!isOperator) {
			nxtTensor.reset(new FullTensor({oldRank, dimensions[numComponents-1], 1}, DONT_SET_ZERO()) );
		} else {
			nxtTensor.reset(new FullTensor({oldRank, dimensions[numComponents-1], dimensions[degree()-1], 1}, DONT_SET_ZERO()) );
		}
		misc::array_copy(static_cast<FullTensor*>(nxtTensor.get())->data.get(), workingData.get(), oldRank*remainingDim);
		
		// set last component tensor to Vt
		set_component(numComponents-1, std::move(nxtTensor));
		
		move_core(0); // TODO create with correct cannonicalization in the first place
		
		REQUIRE((N==1 && remainingDim == dimensions.back()) || (N==2 && remainingDim == dimensions[degree()/2-1]*dimensions[degree()-1]), "Internal Error");
		REQUIRE(is_in_expected_format(), "ie");
	}
	
	
	template<bool isOperator>
	TTNetwork<isOperator>::TTNetwork(const TTNetwork & _cpy) 
		: TensorNetwork(_cpy), cannonicalized(_cpy.cannonicalized), corePosition(_cpy.corePosition)
	{ }
	
	template<bool isOperator>
	TTNetwork<isOperator>::TTNetwork(TTNetwork&& _mov) 
		: TensorNetwork(std::move(_mov)), cannonicalized(_mov.cannonicalized), corePosition(_mov.corePosition) 
	{ }

	template<bool isOperator>
	TTNetwork<isOperator>::TTNetwork(const TensorNetwork &_network, double _eps) : TensorNetwork(_network) {
		LOG(fatal, "Cast of arbitrary tensor network to TT not yet supported"); // TODO
	}

	
	template<bool isOperator>
	TTNetwork<isOperator> TTNetwork<isOperator>::construct_identity(const std::vector<size_t>& _dimensions) {
		REQUIRE(isOperator, "TTTensor identity ill-defined"); // TODO enable if
		REQUIRE(_dimensions.size()%2==0, "Illegal number of dimensions for ttOperator");
		#ifndef DISABLE_RUNTIME_CHECKS_
			for (size_t d : _dimensions) {
				REQUIRE(d > 0, "trying to construct TTOperator with dimension 0");
			}
		#endif
		
		TTNetwork result(_dimensions.size());
		result.factor = 1.0;
		
		const size_t numComponents = _dimensions.size()/N;
		
		for (size_t i=0; i<numComponents; ++i) {
			std::vector<size_t> constructionVector;
			constructionVector.push_back(1);
			for (size_t j=0; j < N; ++j) { 
				constructionVector.push_back(_dimensions[i+j*numComponents]);
			}
			constructionVector.push_back(1);
			result.set_component(i, std::unique_ptr<Tensor>(new FullTensor(constructionVector, [](const std::vector<size_t> &_idx){
				if (_idx[1] == _idx[2]) {
					return 1.0;
				} else {
					return 0.0;
				}
			})));
		}
		
		REQUIRE(result.is_valid_tt(), "Internal Error.");
		result.cannonicalize_left();
		return result;
	}
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	template<bool isOperator>
	void TTNetwork<isOperator>::round_train(const std::vector<size_t>& _maxRanks, const double _eps) {
		REQUIRE(degree()%N==0, "Number of indicis must be even for TTOperator");
		REQUIRE(_eps < 1, "_eps must be smaller than one. " << _eps << " was given.");
		REQUIRE(_maxRanks.size() == degree()/N-1, "There must be exactly degree/N-1 maxRanks. Here " << _maxRanks.size() << " instead of " << degree()/N-1 << " are given.");
		REQUIRE(is_in_expected_format(), "round called on illegal TT tensor");
		
		// If there is no or only one node in the train object we are already finished as there is no rank that can be rounded
		if(degree() <= N) { return; }
		
		// Needed variables
		std::unique_ptr<value_t[]> LR, RR, M, U, S, Vt, newLeft, newRight;
		size_t leftDim, midDim, rightDim, maxLeftRank, maxRightRank, maxRank;
		
		for(size_t position = 0; position < degree()/N-1; ++position) {
			REQUIRE(dynamic_cast<FullTensor*>(nodes[position+1].tensorObject.get()), "Tensor Train nodes are required to be FullTensors for round(...)");
			REQUIRE(dynamic_cast<FullTensor*>(nodes[position+2].tensorObject.get()), "Tensor Train nodes are required to be FullTensors for round(...)");
			
			FullTensor& leftTensor = *static_cast<FullTensor*>(nodes[position+1].tensorObject.get());
			FullTensor& rightTensor = *static_cast<FullTensor*>(nodes[position+2].tensorObject.get());
			
			// Determine the dimensions of the next matrifications
			midDim   = leftTensor.dimensions.back();
			leftDim  = leftTensor.size/midDim;
			rightDim = rightTensor.size/midDim;
			maxLeftRank = std::min(leftDim, midDim);
			// blasWrapper::rank_revealing_split(leftTensor.data, LR, leftTensor.data.get(), leftDim, midDim, maxLeftRank);
			maxRightRank = std::min(midDim, rightDim);
			maxRank = std::min(maxLeftRank, maxRightRank);
			
			// Calculate QR and RQ decompositoins
			leftTensor.ensure_own_data();
			rightTensor.ensure_own_data();
			LR.reset(new value_t[maxLeftRank*midDim]);
			RR.reset(new value_t[midDim*maxRightRank]);
			blasWrapper::inplace_qr(leftTensor.data.get(), LR.get(), leftDim, midDim);
			blasWrapper::inplace_rq(RR.get(), rightTensor.data.get(), midDim, rightDim);
			
			// Calculate Middle Matrix M = LR*RR
			M.reset(new value_t[maxLeftRank*maxRightRank]);
			blasWrapper::matrix_matrix_product(M.get(), maxLeftRank, maxRightRank, 1.0, LR.get(), false, midDim, RR.get(), false);
			
			// Calculate SVD of M = U S Vt -- Reuse the space allocted for LR and RR and allow destruction of M
			U = std::move(LR);
			S.reset(new value_t[maxRank]);
			Vt = std::move(RR);
			blasWrapper::svd_destructive(U.get(), S.get(), Vt.get(), M.get(), maxLeftRank, maxRightRank);
			
			// Determine the Rank NOTE: this terminates as _eps is required to be < 1
			size_t currRank = std::min(maxRank, _maxRanks[position]);
			for(; S[currRank-1] < _eps*S[0]; --currRank) { }
			
			// Calclate S*Vt by scaling the rows of Vt appropriatly
			for(size_t row = 0; row < currRank; ++row) {
				misc::array_scale(Vt.get()+row*maxRightRank, S[row], maxRightRank);
			}
			
			// Multiply U and (S*Vt) back to the left and to the right respectively
			newLeft.reset(new value_t[leftDim*currRank]);
			newRight.reset(new value_t[currRank*rightDim]);
			blasWrapper::matrix_matrix_product(newLeft.get(), leftDim, currRank, 1.0, leftTensor.data.get(), maxLeftRank, false, maxLeftRank, U.get(), maxRank, false);
			blasWrapper::matrix_matrix_product(newRight.get(), currRank, rightDim, 1.0, Vt.get(), false, maxRightRank, rightTensor.data.get(), false);
			
			// Put the new Tensors to their place
			leftTensor.data.reset(newLeft.release(), internal::array_deleter_vt);
			rightTensor.data.reset(newRight.release(), internal::array_deleter_vt);
			leftTensor.dimensions.back() = currRank;
			leftTensor.size = misc::product(leftTensor.dimensions);
			rightTensor.dimensions.front() = currRank;
			rightTensor.size = misc::product(rightTensor.dimensions);
			nodes[position+1].neighbors.back().dimension  = currRank;
			nodes[position+2].neighbors.front().dimension = currRank;
		}
		
		cannonicalized = true;
		corePosition = degree()/N-1;
		REQUIRE(is_in_expected_format(), "ie");
	}
	
	
	template<bool isOperator>
	void TTNetwork<isOperator>::contract_stack(const IndexedTensorWritable<TensorNetwork> &_me) {
		REQUIRE(_me.tensorObject->is_valid_network(), "cannot contract inconsistent ttStack");
		const size_t numComponents = _me.degree()/N;
		const size_t numNodes = _me.degree()/N+2;
		std::set<size_t> toContract;
		for (size_t currentNode=0; currentNode < numNodes; ++currentNode) {
			toContract.clear();
			for (size_t i=currentNode; i<_me.tensorObject->nodes.size(); i+=numNodes) {
				toContract.insert(i);
			}
			_me.tensorObject->contract(toContract);
		}
		// all are contracted, reshuffle them to be in the correct order
		// after contraction the nodes will have one of the ids: node, node+numNodes, node+2*numNodes,... (as those were part of the contraction)
		// so modulus gives the correct wanted id
		_me.tensorObject->reshuffle_nodes([numNodes](size_t i){return i%(numNodes);});
		REQUIRE(_me.tensorObject->nodes.size() == numNodes, "ie");
		REQUIRE(_me.tensorObject->is_valid_network(), "ie: something went wrong in contract_stack");
		
		// reset to new external links
		_me.tensorObject->externalLinks.clear();
		for(size_t i = 0; i < numComponents; ++i) {
			_me.tensorObject->externalLinks.emplace_back(i+1, 1, _me.tensorObject->dimensions[i], false);
		}
		if(N == 2) {
			for(size_t i = 0; i < numComponents; ++i) {
				_me.tensorObject->externalLinks.emplace_back(i+1, 2, _me.tensorObject->dimensions[numComponents+i], false);
			}
		}
		
		// ensure right amount and order of links
		Index ext[N];
		size_t lastRank, externalDim[N], newRank=1;
		std::vector<Index> lastIndices, lastRight;
		std::vector<Index> oldIndices, newRight; // newLeft == lastRight
		std::vector<Index> newIndices;
		std::vector<size_t> newDimensions;
		_me.tensorObject->nodes.front().neighbors = std::vector<TensorNode::Link>({TensorNode::Link(1,0,1,false)});
		_me.tensorObject->nodes.front().tensorObject->reinterpret_dimensions({1});
		_me.tensorObject->nodes.back().neighbors = std::vector<TensorNode::Link>({TensorNode::Link(numComponents,N+1,1,false)});
		_me.tensorObject->nodes.back().tensorObject->reinterpret_dimensions({1});
		for (size_t i=0; i<numComponents; ++i) {
			lastIndices = std::move(oldIndices); oldIndices.clear();
			lastRight = std::move(newRight); newRight.clear();
			lastRank = newRank; newRank=1;
			TensorNode &n = _me.tensorObject->nodes[i+1];
			for (TensorNode::Link &l : n.neighbors) {
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
					LOG(fatal, "ie");
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
				REQUIRE(_me.tensorObject->dimensions[i+j*numComponents] == externalDim[j], "ie");
				n.neighbors.emplace_back(0,i+j*numComponents,externalDim[j], true);
				newDimensions.push_back(externalDim[j]);
			}
			n.neighbors.emplace_back(i+2,0,newRank, false);
			newDimensions.push_back(newRank);
			n.tensorObject->reinterpret_dimensions(newDimensions);
		}
		
		// NOTE core position according to information in TTStack is set in evaluation
		
		REQUIRE(_me.tensorObject->is_in_expected_format(), "something went wrong in contract_stack");
	}
	
	#ifndef DISABLE_RUNTIME_CHECKS_
		template<bool isOperator>
		bool TTNetwork<isOperator>::is_valid_tt() const {
			const size_t numComponents = degree()/N;
			const size_t numNodes = degree()/N + 2;
			REQUIRE(nodes.size() == numNodes, nodes.size() << " vs " << numNodes);
			REQUIRE(externalLinks.size() == degree(), externalLinks.size() << " vs " << degree());
			REQUIRE(std::isfinite(factor), factor);
			REQUIRE(!cannonicalized || corePosition < numComponents, corePosition << " vs " << numComponents);
			
			// per external link
			for (size_t n=0; n<externalLinks.size(); ++n) {
				const TensorNode::Link &l = externalLinks[n];
				REQUIRE(l.dimension == dimensions[n], "n=" << n << " " << l.dimension << " vs " << dimensions[n]);
				REQUIRE(!l.external, "n=" << n);
				REQUIRE(l.other == (n%numComponents)+1, "n=" << n << " " << l.other << " vs " << ((n%numComponents)+1));
				REQUIRE(l.indexPosition < nodes[l.other].neighbors.size(), "n=" << n << " " << l.indexPosition << " vs " << nodes[l.other].neighbors.size());
				REQUIRE(nodes[l.other].neighbors[l.indexPosition].external, "n=" << n);
				REQUIRE(nodes[l.other].neighbors[l.indexPosition].indexPosition == n, "n=" << n << " " << nodes[l.other].neighbors[l.indexPosition].indexPosition);
				REQUIRE(nodes[l.other].neighbors[l.indexPosition].dimension == l.dimension, "n=" << n << " " << nodes[l.other].neighbors[l.indexPosition].dimension << " vs " << l.dimension);
			}
			
			// virtual nodes
			REQUIRE(nodes.front().neighbors.size() == 1, nodes.front().neighbors.size());
			REQUIRE(nodes.front().neighbors[0].dimension == 1, nodes.front().neighbors[0].dimension);
			REQUIRE(nodes.front().neighbors[0].other == 1, nodes.front().neighbors[0].other);
			REQUIRE(nodes.front().neighbors[0].indexPosition == 0, nodes.front().neighbors[0].indexPosition);
			if (nodes.front().tensorObject) {
				REQUIRE(nodes.front().tensorObject->dimensions == std::vector<size_t>({1}), nodes.front().tensorObject->dimensions);
				#pragma GCC diagnostic push
				#pragma GCC diagnostic ignored "-Wfloat-equal"
				REQUIRE((*nodes.front().tensorObject)[0] == 1, (*nodes.front().tensorObject)[0]);
				#pragma GCC diagnostic pop
				REQUIRE(!nodes.front().tensorObject->has_factor(), "");
			}
			
			REQUIRE(nodes.back().neighbors.size() == 1, nodes.back().neighbors.size());
			REQUIRE(nodes.back().neighbors[0].dimension == 1, nodes.back().neighbors[0].dimension);
			REQUIRE(nodes.back().neighbors[0].other == numNodes-2, nodes.back().neighbors[0].other);
			REQUIRE(nodes.back().neighbors[0].indexPosition == N+1, nodes.back().neighbors[0].indexPosition);
			if (nodes.back().tensorObject) {
				REQUIRE(nodes.back().tensorObject->dimensions == std::vector<size_t>({1}), nodes.back().tensorObject->dimensions);
				#pragma GCC diagnostic push
				#pragma GCC diagnostic ignored "-Wfloat-equal"
				REQUIRE((*nodes.back().tensorObject)[0] == 1, (*nodes.back().tensorObject)[0]);
				#pragma GCC diagnostic pop
				REQUIRE(!nodes.back().tensorObject->has_factor(), "");
			}
			
			// per component
			for (size_t n=0; n<numComponents; ++n) {
				const TensorNode &node = nodes[n+1];
				REQUIRE(!node.erased, "n=" << n);
				REQUIRE(node.degree() == N+2, "n=" << n << " " << node.degree());
				if (node.tensorObject) {
					REQUIRE(node.tensorObject->degree() == N+2, "n=" << n << " " << node.tensorObject->degree());
				}
				REQUIRE(!node.neighbors[0].external, "n=" << n);
				REQUIRE(node.neighbors[0].other == n, "n=" << n);
				REQUIRE(node.neighbors[0].indexPosition == (n==0?0:N+1), "n=" << n << " " << node.neighbors[0].indexPosition);
				REQUIRE(node.neighbors[1].external, "n=" << n);
				REQUIRE(node.neighbors[1].indexPosition == n, "n=" << n << " " << node.neighbors[1].indexPosition);
				if (isOperator) {
					REQUIRE(node.neighbors[2].external, "n=" << n);
					REQUIRE(node.neighbors[2].indexPosition == numComponents+n, "n=" << n << " " << node.neighbors[2].indexPosition << " vs " << numComponents+n);
				}
				if (node.tensorObject) {
					if (cannonicalized && n!=corePosition) {
						REQUIRE(!node.tensorObject->has_factor(), "n="<<n);
					}
					REQUIRE(node.tensorObject->dimensions[0]==node.neighbors[0].dimension, "n=" << n);
					REQUIRE(node.tensorObject->dimensions[1]==node.neighbors[1].dimension, "n=" << n);
					REQUIRE(node.tensorObject->dimensions[2]==node.neighbors[2].dimension, "n=" << n);
					if (isOperator) {
						REQUIRE(node.tensorObject->dimensions[3]==node.neighbors[3].dimension, "n=" << n);
					}
				}
				REQUIRE(!node.neighbors.back().external, "n=" << n);
				REQUIRE(node.neighbors.back().other == n+2, "n=" << n << " " << node.neighbors.back().other);
				REQUIRE(node.neighbors.back().indexPosition == 0, "n=" << n << " " << node.neighbors.back().indexPosition);
				REQUIRE(!nodes[n+2].neighbors.empty(), "n=" << n);
				REQUIRE(node.neighbors.back().dimension == nodes[n+2].neighbors[0].dimension, "n=" << n << " " << node.neighbors.back().dimension << " vs " << nodes[n+1].neighbors[0].dimension);
			}
			
			return true;
		}
	#else
		template<bool isOperator>
		bool TTNetwork<isOperator>::is_valid_tt() const {
			return true;
		}
	#endif
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	template<bool isOperator>
	const Tensor &TTNetwork<isOperator>::get_component(size_t _idx) const {
		REQUIRE(_idx < degree()/N, "illegal index in TTNetwork::get_component");
		return *nodes[_idx+1].tensorObject;
	}
	
	
	template<bool isOperator>
	void TTNetwork<isOperator>::set_component(size_t _idx, const Tensor &_T) {
		REQUIRE(_idx < degree()/N, "illegal index in TTNetwork::set_component");
		TensorNode &currNode = nodes[_idx+1];
		REQUIRE(_T.degree() == N+2, "Component must have degree 3 (TTTensor) or 4 (TTOperator). Given: " << _T.degree());
		REQUIRE(_T.degree() == currNode.degree(), "Degree of _T does not match component tensors degree");
		currNode.tensorObject.reset(_T.get_copy());
		for (size_t i=0; i<currNode.degree(); ++i) {
			currNode.neighbors[i].dimension = currNode.tensorObject->dimensions[i];
			if (currNode.neighbors[i].external) {
				externalLinks[currNode.neighbors[i].indexPosition].dimension = currNode.tensorObject->dimensions[i];
				dimensions[currNode.neighbors[i].indexPosition] = currNode.tensorObject->dimensions[i];
			}
		}
		if (corePosition != _idx) {
			cannonicalized = false;
		}
	}
	
	template<bool isOperator>
	void TTNetwork<isOperator>::set_component(size_t _idx, std::unique_ptr<Tensor> &&_T) {
		REQUIRE(_idx < degree()/N, "Illegal index in TTNetwork::set_component");
		TensorNode &currNode = nodes[_idx+1];
		REQUIRE(_T->degree() == N+2, "Component must have degree 3 (TTTensor) or 4 (TTOperator). Given: " << _T->degree());
		REQUIRE(_T->degree() == currNode.degree(), "Degree of _T does not match component tensors degree");
		currNode.tensorObject = std::move(_T);
		for (size_t i=0; i<currNode.degree(); ++i) {
			currNode.neighbors[i].dimension = currNode.tensorObject->dimensions[i];
			if (currNode.neighbors[i].external) {
				externalLinks[currNode.neighbors[i].indexPosition].dimension = currNode.tensorObject->dimensions[i];
				dimensions[currNode.neighbors[i].indexPosition] = currNode.tensorObject->dimensions[i];
			}
		}
		if (corePosition != _idx) {
			cannonicalized = false;
		}
	}
	
	template<bool isOperator>
	TTNetwork<isOperator> TTNetwork<isOperator>::dyadic_product(const TTNetwork<isOperator> &_lhs, const TTNetwork<isOperator> &_rhs) {
		REQUIRE(_lhs.is_in_expected_format(), "");
		REQUIRE(_rhs.is_in_expected_format(), "");
		
		if (_lhs.degree() == 0) {
			TTNetwork result(_rhs);
			result.factor *= _lhs.factor;
			return result;
		}
		TTNetwork result(_lhs);
		result.factor *= _rhs.factor;
		if (_rhs.degree() == 0) {
			return result;
		}
		
		// add all nodes of rhs and fix neighbor relations
		result.nodes.pop_back();
		const size_t lhsNumComponents = _lhs.degree()/N;
		const size_t rhsNodesSize = _rhs.nodes.size();
		const size_t rhsNumComponents = _rhs.degree()/N;
		for (size_t i=1; i<rhsNodesSize; ++i) {
			const TensorNode &n = _rhs.nodes[i];
			result.nodes.emplace_back(n);
			for (TensorNode::Link &l : result.nodes.back().neighbors) {
				if (l.external) {
					if (l.indexPosition < rhsNumComponents) {
						l.indexPosition += lhsNumComponents;
					} else {
						l.indexPosition += 2*lhsNumComponents;
					}
				} else {
					if (l.other==0) {
						l.indexPosition = N+1;
					}
					l.other += lhsNumComponents;
				}
			}
		}
		
		// add all external indices of rhs
		result.externalLinks.clear();
		result.dimensions.clear();
		for (size_t i=0; i<lhsNumComponents; ++i) {
			const size_t d=_lhs.dimensions[i];
			result.externalLinks.emplace_back(i+1, 1, d, false);
			result.dimensions.push_back(d);
		}
		for (size_t i=0; i<rhsNumComponents; ++i) {
			const size_t d=_rhs.dimensions[i];
			result.externalLinks.emplace_back(lhsNumComponents+i+1, 1, d, false);
			result.dimensions.push_back(d);
		}
		if (isOperator) {
			for (size_t i=0; i<lhsNumComponents; ++i) {
				const size_t d=_lhs.dimensions[i];
				result.externalLinks.emplace_back(i+1, 2, d, false);
				result.dimensions.push_back(d);
			}
			for (size_t i=0; i<rhsNumComponents; ++i) {
				const size_t d=_rhs.dimensions[i];
				result.externalLinks.emplace_back(lhsNumComponents+i+1, 2, d, false);
				result.dimensions.push_back(d);
			}
		}
		
		if (_lhs.cannonicalized && _rhs.cannonicalized) {
			if (_lhs.corePosition == 0 && _rhs.corePosition == 0) {
				result.cannonicalized = true;
				result.corePosition = lhsNumComponents;
				result.move_core(0);
			}
			if (_lhs.corePosition == lhsNumComponents-1 && _rhs.corePosition == rhsNumComponents-1) {
				result.cannonicalized = true;
				result.corePosition = lhsNumComponents-1;
				result.move_core(lhsNumComponents + rhsNumComponents -1);
			}
		} else {
			result.cannonicalized = false;
		}
		
		REQUIRE(result.is_valid_tt(), "ie");
		return result;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator> TTNetwork<isOperator>::dyadic_product(const std::vector<std::reference_wrapper<TTNetwork<isOperator>>> &_tensors) {
		if (_tensors.size() == 0) {
			return TTNetwork();
		} 
		TTNetwork result(_tensors.back());
		// construct dyadic products right to left as default cannonicalization is left
		for (size_t i=_tensors.size()-1; i>0; --i) {
			REQUIRE_TEST;
			result = dyadic_product(_tensors[i-1], result);
		}
		return result;
	}
	
	
	template<bool isOperator>
	std::pair<TensorNetwork, TensorNetwork> TTNetwork<isOperator>::chop(const size_t _position) const {
		REQUIRE(is_valid_tt(), "Invalid TT cannot be chopped.");
		
		const size_t numComponents = degree()/N;
		REQUIRE(_position < numComponents, "Can't split a " << numComponents << " component TTNetwork at position " << _position);
		
		// Create the resulting TNs
		TensorNetwork left, right;
		left.factor = 1.;
		right.factor = 1.;
		
		left.nodes.push_back(nodes[0]);
		for (size_t i = 0; i < _position; ++i) {
			left.dimensions.push_back(dimensions[i]);
			left.externalLinks.push_back(externalLinks[i]);
			left.nodes.push_back(nodes[i+1]);
		}
		if(isOperator) {
			for(size_t i = 0; i < _position; ++i) {
				left.dimensions.push_back(dimensions[i+numComponents]);
				left.externalLinks.push_back(externalLinks[i+numComponents]);
			}
		}
		left.dimensions.push_back(left.nodes.back().neighbors.back().dimension);
		left.externalLinks.emplace_back(_position, _position==0?0:N+1, left.nodes.back().neighbors.back().dimension , false);
		left.nodes.back().neighbors.back().external = true;
		left.nodes.back().neighbors.back().indexPosition = isOperator ? 2*_position-1 : _position;
		
		right.dimensions.push_back(nodes[_position+2].neighbors.front().dimension);
		right.externalLinks.emplace_back(_position+2, 0, nodes[_position+2].neighbors.front().dimension , false); // NOTE other will be corrected to 0 in the following steps
		for(size_t i = _position+1; i < numComponents+1; ++i) {
			right.dimensions.push_back(dimensions[i]);
			right.externalLinks.push_back(externalLinks[i]);
			right.nodes.push_back(nodes[i+1]);
		}
		if(isOperator) {
			for(size_t i = _position+1; i < numComponents+1; ++i) {
				right.dimensions.push_back(dimensions[i+numComponents]);
				right.externalLinks.push_back(externalLinks[i+numComponents]);
			}
		}
		right.nodes.front().neighbors.front().external = true;
		right.nodes.front().neighbors.front().indexPosition = _position; // NOTE indexPosition will be corrected to 0 in the following steps
		
		// Account for the fact that the first _position+2 nodes do not exist
		for(TensorNode::Link& link : right.externalLinks) {
			link.other -= _position+2;
		}
		
		for(TensorNode& node : right.nodes) {
			for(TensorNode::Link& link : node.neighbors) {
				if(link.external) {
					link.indexPosition -= _position+1;
				} else {
					link.other -= _position+2;
				}
			}
		}
		
		REQUIRE(left.is_valid_network(), "Internal Error");
		REQUIRE(right.is_valid_network(), "Internal Error");
		
		return std::pair<TensorNetwork, TensorNetwork>(std::move(left), std::move(right));
	}
	
	
	template<bool isOperator>
	void TTNetwork<isOperator>::round(value_t _eps) {
		const bool oldCannon = cannonicalized;
		const size_t oldCorePos = corePosition;
		move_core(0);
		round_train(std::vector<size_t>(degree()/N-1, size_t(-1)), _eps);
		if (oldCannon) {
			move_core(oldCorePos);
		}
	}

	template<bool isOperator>
	void TTNetwork<isOperator>::round(size_t _maxRank) {
		const bool oldCannon = cannonicalized;
		const size_t oldCorePos = corePosition;
		move_core(0);
		round_train(std::vector<size_t>(degree()/N-1, _maxRank), 1e-15);
		if (oldCannon) {
			move_core(oldCorePos);
		}
	}

	template<bool isOperator>
	void TTNetwork<isOperator>::round(const std::vector<size_t> &_maxRanks) {
		const bool oldCannon = cannonicalized;
		const size_t oldCorePos = corePosition;
		move_core(0);
		round_train(_maxRanks, 1e-15);
		if (oldCannon) {
			move_core(oldCorePos);
		}
	}
	
	template<bool isOperator>
	void TTNetwork<isOperator>::round(int _maxRank) {
		REQUIRE( _maxRank > 0, "MaxRank must be positive");
		round(size_t(_maxRank));
	}

	template<bool isOperator>
	std::vector<size_t> TTNetwork<isOperator>::ranks() const {
		std::vector<size_t> res;
		for (size_t n=1; n<nodes.size()-2; ++n) {
			res.push_back(nodes[n].neighbors.back().dimension);
		}
		return res;
	}
	
	template<bool isOperator>
	size_t TTNetwork<isOperator>::rank(size_t _i) const {
		REQUIRE(_i+2 < nodes.size(), "requested illegal rank");
		return nodes[_i+1].neighbors.back().dimension;
	}
	
	template<bool isOperator>
	size_t TTNetwork<isOperator>::datasize() const {
		size_t result = 0;
		for (const TensorNode &n : nodes) {
			result += n.tensorObject->size;
		}
		return result;
	}
	
	template<bool isOperator>
	void TTNetwork<isOperator>::move_core(size_t _position, bool _keepRank) {
		REQUIRE(is_valid_tt(), "Cannot cannonicalize invalid TT");
		const size_t numComponents = degree()/N;
		REQUIRE(_position < numComponents, "Illegal position for core chosen");
		Index i,r,j;
		FullTensor core(2);
		const size_t fromLeft = cannonicalized? corePosition : 0;
		const size_t fromRight = cannonicalized? corePosition : (numComponents - 1);
		// move right?
		for (size_t n=fromLeft; n<_position; ++n) {
			Tensor &currTensor = *nodes[n+1].tensorObject;
			if (_keepRank) {
				( currTensor(i&1,r), core(r,j) ) = QR(currTensor(i&1,j));
			} else {
				( currTensor(i&1,r), core(r,j) ) = OrthogonalSplit(currTensor(i&1,j));
			}
			
			Tensor &nextTensor = *nodes[n+2].tensorObject;
			nextTensor(i,j&1) = core(i,r) * nextTensor(r,j&1);
			if (nextTensor.dimensions[0] != nodes[n+1].neighbors.front().dimension) {
				nodes[n+2].neighbors.front().dimension = nodes[n+1].neighbors.back().dimension = nextTensor.dimensions[0];
			}
		}
		for (size_t n=fromRight; n > _position; --n) {
			Tensor &currTensor = *nodes[n+1].tensorObject;
			if (_keepRank) {
				( core(j,r), currTensor(r,i&1) ) = RQ(currTensor(j,i&1));
			} else {
				( currTensor(i&1,r), core(r,j) ) = OrthogonalSplit(currTensor(j,i&1));
				currTensor(r, i&1) = currTensor(i&1, r);
				core(j,r) = core(r,j);
			}
			
			Tensor &nextTensor = *nodes[n].tensorObject;
			nextTensor(j&1,i) = nextTensor(j&1,r) * core(r,i);
			if (currTensor.dimensions[0] != nodes[n+1].neighbors.front().dimension) {
				nodes[n+1].neighbors.front().dimension = nodes[n].neighbors.back().dimension = currTensor.dimensions[0];
			}
		}
		cannonicalized = true;
		corePosition = _position;
	}
	
	template<bool isOperator>
	void TTNetwork<isOperator>::cannonicalize_left() {
		move_core(0);
	}
	
	template<bool isOperator>
	void TTNetwork<isOperator>::cannonicalize_right() {
		move_core(degree()/N-1);
	}
	
	template<bool isOperator>
	TensorNetwork* TTNetwork<isOperator>::get_copy() const {
		return new TTNetwork(*this);
	}
	
	template<bool isOperator>
	value_t TTNetwork<isOperator>::frob_norm() const {
		REQUIRE(is_valid_tt(), "frob_norm of illegal TT");
		if (cannonicalized) {
			return get_component(corePosition).frob_norm();
		} else {
			Index i;
			return (*((*this)(i&0)*(*this)(i&0)).tensorObject)[0];
		}
	}
	
	template<bool isOperator>
	bool TTNetwork<isOperator>::is_in_expected_format() const {
		return is_valid_tt();
	}
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - -  Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	template<bool isOperator>
	TTNetwork<isOperator>& TTNetwork<isOperator>::operator+=(const TTNetwork<isOperator>& _other) {
		Index i;
		(*this)(i&0) = (*this)(i&0) + _other(i&0);
		return *this;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>  TTNetwork<isOperator>::operator+(const TTNetwork<isOperator>& _other) const {
		TTNetwork cpy(*this);
		cpy += _other;
		return cpy;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>& TTNetwork<isOperator>::operator-=(const TTNetwork<isOperator>& _other) {
		Index i;
		(*this)(i&0) = (*this)(i&0) - _other(i&0);
		return *this;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>  TTNetwork<isOperator>::operator-(const TTNetwork<isOperator>& _other) const {
		TTNetwork cpy(*this);
		cpy -= _other;
		return cpy;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>& TTNetwork<isOperator>::operator*=(const value_t _prod) {
		factor *= _prod;
		return *this;
		
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>  TTNetwork<isOperator>::operator*(const value_t _prod) const {
		TTNetwork result(*this);
		result *= _prod;
		return result;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>& TTNetwork<isOperator>::operator/=(const value_t _div) {
		factor /= _div;
		return *this;
	}
	
	template<bool isOperator>
	TTNetwork<isOperator>  TTNetwork<isOperator>::operator/(const value_t _div) const {
		TTNetwork result(*this);
		result /= _div;
		return result;
	}
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
	
	template<bool isOperator>
	bool TTNetwork<isOperator>::specialized_contraction_f(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) {
		REQUIRE(!_out.tensorObject, "Internal Error.");
		
		// Only TTOperators construct stacks, so no specialized contractions for TTTensors
		if(!isOperator) { return false; }
		
		const std::vector<Index> myIndices = _me.get_assigned_indices();
		const std::vector<Index> otherIndices = _other.get_assigned_indices();
		
		const TTNetwork* const meTT = dynamic_cast<const TTNetwork*>(_me.tensorObjectReadOnly);
		const internal::TTStack<true>* const meTTStack = dynamic_cast<const internal::TTStack<true>*>(_me.tensorObjectReadOnly);
		REQUIRE(meTT || meTTStack, "ie");
		
		const TTTensor* const otherTT = dynamic_cast<const TTTensor*>(_other.tensorObjectReadOnly);
		const internal::TTStack<false>* const otherTTStack = dynamic_cast<const internal::TTStack<false>*>(_other.tensorObjectReadOnly);
		const TTOperator* const otherTTO = dynamic_cast<const TTOperator*>(_other.tensorObjectReadOnly);
		const internal::TTStack<true>* const otherTTOStack = dynamic_cast<const internal::TTStack<true>*>(_other.tensorObjectReadOnly);
		
		if (!otherTT && !otherTTStack && !otherTTO && !otherTTOStack) {
			return false;
		}
		
		bool cannoAtTheEnd = false;
		size_t coreAtTheEnd = 0;
		if (meTT) {
			cannoAtTheEnd = meTT->cannonicalized;
			coreAtTheEnd = meTT->corePosition;
		} else {
			cannoAtTheEnd = meTTStack->cannonicalization_required;
			coreAtTheEnd = meTTStack->futureCorePosition;
		}
		
		
		// TODO profiler should warn if other->corePosition is not identical to coreAtTheEnd
		
		// Determine my first half and second half of indices
		auto midIndexItr = myIndices.begin();
		size_t spanSum = 0;
		while (spanSum < _me.degree() / 2) {
			REQUIRE(midIndexItr != myIndices.end(), "ie");
			spanSum += midIndexItr->span;
			++midIndexItr;
		}
		if (spanSum > _me.degree() / 2) {
			return false; // an index spanned some links of the left and some of the right side
		}
		
		if (otherTT) {
			// ensure fitting indices
			if (misc::equal(myIndices.begin(), midIndexItr, otherIndices.begin(), otherIndices.end()) || misc::equal(midIndexItr, myIndices.end(), otherIndices.begin(), otherIndices.end())) {
				TensorNetwork *res = new internal::TTStack<false>(cannoAtTheEnd, coreAtTheEnd);
				*res = *_me.tensorObjectReadOnly;
				res->factor *= _other.tensorObjectReadOnly->factor;
				_out.reset(res, myIndices, true);
				TensorNetwork::add_network_to_network(_out, _other);
				return true;
			} else {
				return false;
			}
		} else { // other is operator or operator stack
			// determine other middle index
			auto otherMidIndexItr = otherIndices.begin();
			spanSum = 0;
			while (spanSum < _other.degree() / 2) {
				REQUIRE(otherMidIndexItr != otherIndices.end(), "ie");
				spanSum += otherMidIndexItr->span;
				++otherMidIndexItr;
			}
			if (spanSum > _other.degree() / 2) {
				return false; // an index spanned some links of the left and some of the right side
			}
			// or indices in fitting order to contract the TTOs
			if (   misc::equal(myIndices.begin(), midIndexItr, otherIndices.begin(), otherMidIndexItr) 
				|| misc::equal(midIndexItr, myIndices.end(), otherIndices.begin(), otherMidIndexItr)
				|| misc::equal(myIndices.begin(), midIndexItr, otherMidIndexItr, otherIndices.end()) 
				|| misc::equal(midIndexItr, myIndices.end(), otherMidIndexItr, otherIndices.end())	) 
			{
				TensorNetwork *res = new internal::TTStack<true>(cannoAtTheEnd, coreAtTheEnd);
				*res = *_me.tensorObjectReadOnly;
				res->factor *= _other.tensorObjectReadOnly->factor;
				_out.reset(res, myIndices, true);
				TensorNetwork::add_network_to_network(_out, _other);
				return true;
			} else {
				return false;
			}
		}
	}
	
	template<bool isOperator>
	bool TTNetwork<isOperator>::specialized_sum_f(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) {
		const std::vector<Index> myIndices = _me.get_assigned_indices();
		const std::vector<Index> otherIndices = _other.get_assigned_indices();
		
		// If the indices are in different order, we are lost. TODO inverse order is also ok...
		if (myIndices != otherIndices) { return false; }
		REQUIRE(_me.tensorObjectReadOnly->dimensions == _other.tensorObjectReadOnly->dimensions, "TT sum requires both operants to share the same dimensions");
		
		// If the other is not a TT tensor (or stack) we are also lost
		const TTNetwork* otherTT = dynamic_cast<const TTNetwork*>( _other.tensorObjectReadOnly);
		const internal::TTStack<isOperator>* otherTTStack = dynamic_cast<const internal::TTStack<isOperator>*>( _other.tensorObjectReadOnly);
		if (!otherTT && !otherTTStack) { return false; }
		
		// TODO the order is not canonical, because if I am no Stack I don't have to know whether or not i am moveable
		// If I am in fact a TTTensorStack, we have to evaluate me to TTNetwork
		std::unique_ptr<IndexedTensor<TensorNetwork>> meStorage;
		const IndexedTensorReadOnly<TensorNetwork> *realMePtr = &_me;
		const IndexedTensorMoveable<TensorNetwork> *movMe = dynamic_cast<const IndexedTensorMoveable<TensorNetwork> *>(&_me);
		if (movMe) {
			internal::TTStack<isOperator> *stackMe = dynamic_cast<internal::TTStack<isOperator> *>(movMe->tensorObject);
			if (stackMe) {
				meStorage.reset(new IndexedTensor<TensorNetwork>(new TTNetwork(_me.degree()), myIndices, true));
				(*meStorage) = _me;
				realMePtr = meStorage.get();
			}
		} else {
			REQUIRE(!dynamic_cast<const internal::TTStack<isOperator> *>(_me.tensorObjectReadOnly),"ie - non-moveable TTStack detected");
		}
		const IndexedTensorReadOnly<TensorNetwork> &realMe = *realMePtr;
		
		// If other is in fact a TTTensorStack, we have to evaluate it to tttensor
		std::unique_ptr<IndexedTensor<TensorNetwork>> otherStorage;
		const IndexedTensorReadOnly<TensorNetwork> *realOtherPtr = &_other;
		const IndexedTensorMoveable<TensorNetwork> *movOther = dynamic_cast<const IndexedTensorMoveable<TensorNetwork> *>(&_other);
		if (movOther) {
			internal::TTStack<isOperator> *stackOther = dynamic_cast<internal::TTStack<isOperator> *>(movOther->tensorObject);
			if (stackOther) {
				otherStorage.reset(new IndexedTensor<TensorNetwork>(new TTNetwork(_other.degree()), otherIndices, true));
				(*otherStorage) = _other;
				realOtherPtr = otherStorage.get();
			}
		} else {
			REQUIRE(!dynamic_cast<const internal::TTStack<isOperator> *>(_other.tensorObjectReadOnly),"ie - non-moveable TTStack detected");
		}
		const IndexedTensorReadOnly<TensorNetwork> &realOther = *realOtherPtr;
		
		// Number of components to create
		const size_t numComponents = realMe.degree()/N;
		
		TTNetwork* tmpPtr = new TTNetwork(realMe.degree());
		tmpPtr->factor = 1.0;
		
		//The external dimensions are the same as the ones of the input
		tmpPtr->dimensions = realMe.tensorObjectReadOnly->dimensions;
		REQUIRE(realOther.tensorObjectReadOnly->dimensions == realMe.tensorObjectReadOnly->dimensions, "Internal Error");
		
		IndexedTensor<TensorNetwork> tmpOut(tmpPtr, myIndices, true);
		TTNetwork& outTensor = *static_cast<TTNetwork*>(tmpOut.tensorObject);
		
		if (numComponents == 0) {
			outTensor.factor = _me.tensorObjectReadOnly->factor + _other.tensorObjectReadOnly->factor;
			return true;
		}
		
		if(numComponents == 1) {
			// Create the one Node
			std::unique_ptr<Tensor> nextTensor;
			const Tensor &myComponent = *realMe.tensorObjectReadOnly->nodes[1].tensorObject.get();
			const Tensor &otherComponent = *realOther.tensorObjectReadOnly->nodes[1].tensorObject.get();
			if(myComponent.is_sparse() && otherComponent.is_sparse()) { // Both Sparse
				nextTensor.reset(myComponent.get_copy());
				nextTensor->factor *= realMe.tensorObjectReadOnly->factor;
				*static_cast<SparseTensor*>(nextTensor.get()) += realOther.tensorObjectReadOnly->factor*(*static_cast<const SparseTensor*>(&otherComponent));
			} else { // at most one sparse
				if(myComponent.is_sparse()){
					nextTensor.reset(new FullTensor(*static_cast<const SparseTensor*>(&myComponent)));
				} else {
					nextTensor.reset(new FullTensor(*static_cast<const FullTensor*>(&myComponent)));
				}
				nextTensor->factor *= realMe.tensorObjectReadOnly->factor;
				if(otherComponent.is_sparse()){
					*static_cast<FullTensor*>(nextTensor.get()) += realOther.tensorObjectReadOnly->factor*static_cast<const SparseTensor&>(otherComponent);
				} else {
					*static_cast<FullTensor*>(nextTensor.get()) += realOther.tensorObjectReadOnly->factor*static_cast<const FullTensor&>(otherComponent);
				}
			}
			
			outTensor.set_component(0, std::move(nextTensor));
			_out.assign(std::move(tmpOut));
			return true;
		}
		
		const TTNetwork * const ttMe = static_cast<const TTNetwork*>(realMe.tensorObjectReadOnly);
		const TTNetwork * const ttOther = static_cast<const TTNetwork*>(realOther.tensorObjectReadOnly);
		
		PA_START;
		for(size_t position = 0; position < numComponents; ++position) {
			// Get current input nodes
			// TODO sparse
			FullTensor &myComponent = *static_cast<FullTensor*>(realMe.tensorObjectReadOnly->nodes[position+1].tensorObject.get());
			FullTensor &otherComponent = *static_cast<FullTensor*>(realOther.tensorObjectReadOnly->nodes[position+1].tensorObject.get());
			
			// Structure has to be (for degree 4)
			// (L1 R1) * ( L2 0  ) * ( L3 0  ) * ( L4 )
			// 			 ( 0  R2 )   ( 0  R3 )   ( R4 )
			
			// Create a FullTensor for Node
			std::vector<size_t> nxtDimensions;
			if (position == 0) { 
				nxtDimensions.emplace_back(1);
			} else {
				nxtDimensions.emplace_back(myComponent.dimensions.front()+otherComponent.dimensions.front());
			}
			nxtDimensions.emplace_back(outTensor.dimensions[position]);
			if (isOperator) { nxtDimensions.emplace_back(outTensor.dimensions[position+numComponents]); }
			if (position == numComponents-1) {
				nxtDimensions.emplace_back(1);
			} else {
				nxtDimensions.emplace_back(myComponent.dimensions.back()+otherComponent.dimensions.back());
			}
			
			std::unique_ptr<Tensor> newComponent(new FullTensor(std::move(nxtDimensions)) );
			value_t * const componentData = static_cast<FullTensor*>(newComponent.get())->data.get();
			
			
			const size_t leftIdxOffset = newComponent->size/newComponent->dimensions.front();
			const size_t extIdxOffset = newComponent->dimensions.back();
			const size_t myLeftIdxOffset = myComponent.size/myComponent.dimensions.front();
			const size_t myExtIdxOffset = myComponent.dimensions.back();
			const size_t otherLeftIdxOffset = otherComponent.size/otherComponent.dimensions.front();
			const size_t otherExtIdxOffset = otherComponent.dimensions.back();
			const size_t otherGeneralOffset = (position == 0 ? 0 : myComponent.dimensions.front()*leftIdxOffset) + (position == numComponents-1 ? 0 : myComponent.dimensions.back());
			const size_t extDimSize = myComponent.dimensions[1] * (isOperator? myComponent.dimensions[2] : 1);
			
			// Copy own Tensor into place
			if (!ttMe->cannonicalized || position == ttMe->corePosition) {
				for(size_t leftIdx = 0; leftIdx < myComponent.dimensions.front(); ++leftIdx) {
					for(size_t extIdx = 0; extIdx < extDimSize; ++extIdx) {
						// RightIdx can be copied in one piece
						misc::array_scaled_copy(componentData + leftIdx*leftIdxOffset + extIdx*extIdxOffset, 
												myComponent.factor*realMe.tensorObjectReadOnly->factor, 
												myComponent.data.get() + leftIdx*myLeftIdxOffset + extIdx*myExtIdxOffset, 
												myComponent.dimensions.back());
					}
				}
			} else {
				REQUIRE(!myComponent.has_factor(), "Only Core node is allowed to have a factor");
				for(size_t leftIdx = 0; leftIdx < myComponent.dimensions.front(); ++leftIdx) {
					for(size_t extIdx = 0; extIdx < extDimSize; ++extIdx) {
						// RightIdx can be copied as one piece
						misc::array_copy(componentData + leftIdx*leftIdxOffset + extIdx*extIdxOffset, 
										 myComponent.data.get() + leftIdx*myLeftIdxOffset + extIdx*myExtIdxOffset, 
										 myComponent.dimensions.back());
					}
				}
			}
			
			
			// Copy other Tensor into place
			if(!ttOther->cannonicalized || position == ttOther->corePosition) {
				for(size_t leftIdx = 0; leftIdx < otherComponent.dimensions.front(); ++leftIdx) {
					for(size_t extIdx = 0; extIdx < extDimSize; ++extIdx) {
						// RightIdx can be copied as one piece
						misc::array_scaled_copy(componentData + leftIdx*leftIdxOffset + extIdx*extIdxOffset + otherGeneralOffset, 
												otherComponent.factor*realOther.tensorObjectReadOnly->factor, 
												otherComponent.data.get() + leftIdx*otherLeftIdxOffset + extIdx*otherExtIdxOffset, 
												otherComponent.dimensions.back());
					}
				}
			} else {
				REQUIRE(!otherComponent.has_factor(), "Only Core node is allowed to have a factor");
				for(size_t leftIdx = 0; leftIdx < otherComponent.dimensions.front(); ++leftIdx) {
					for(size_t extIdx = 0; extIdx < extDimSize; ++extIdx) {
						// RightIdx can be copied as one piece
						misc::array_copy(componentData + leftIdx*leftIdxOffset + extIdx*extIdxOffset + otherGeneralOffset, 
										 otherComponent.data.get() + leftIdx*otherLeftIdxOffset + extIdx*otherExtIdxOffset, 
										 otherComponent.dimensions.back());
					}
				}
			}
			
			outTensor.set_component(position, std::move(newComponent));
		}
		
		PA_END("ADD/SUB", "TTNetwork ADD/SUB", std::string("Dims:")+misc::to_string(outTensor.dimensions)+" Ranks: "+misc::to_string(outTensor.ranks()));
		
		// TODO profiler should warn if other->corePos differs
		
		if (ttMe->cannonicalized) {
			REQUIRE(!outTensor.cannonicalized, "ie");
			outTensor.move_core(ttMe->corePosition);
		}
		
		REQUIRE(outTensor.is_valid_tt(), "ie");
		_out.assign(std::move(tmpOut));
		return true;
	}
	
	
	template<bool isOperator>
	void TTNetwork<isOperator>::specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) {
		REQUIRE(_me.tensorObject == this, "Internal Error.");
		
		const std::vector<Index> myIndices = _me.get_assigned_indices(_other.degree()); // TODO this wont work if we have fixed indices in TT tensors.
		const std::vector<Index> otherIndices = _other.get_assigned_indices();
		const size_t numComponents = _other.degree()/N;
		
		// First check whether the other is a TTNetwork as well, otherwise we can skip to fallback
		const TTNetwork* const otherTTN = dynamic_cast<const TTNetwork*>(_other.tensorObjectReadOnly);
		TTNetwork* const meTTN = dynamic_cast<TTNetwork*>(_me.tensorObject);
		REQUIRE(meTTN, "ie");
		const internal::TTStack<isOperator>* const otherTTStack = dynamic_cast<const internal::TTStack<isOperator>*>(_other.tensorObjectReadOnly);
		const IndexedTensorMoveable<TensorNetwork> *movOther = dynamic_cast<const IndexedTensorMoveable<TensorNetwork> *>(&_other);
		if (otherTTStack) {
			REQUIRE(movOther, "not moveable TTStack encountered...");
			contract_stack(*movOther);
		}
		if(otherTTN || otherTTStack) {
			// Check whether the index order coincides
			if (myIndices == otherIndices) {
				if (otherTTN) {
					*meTTN = *otherTTN;
				} else {
					static_cast<TensorNetwork*>(_me.tensorObject)->operator=(*_other.tensorObjectReadOnly);
					if (otherTTStack->cannonicalization_required) {
						meTTN->move_core(otherTTStack->futureCorePosition);
					} else {
						meTTN->cannonicalized = false;
					}
				}
				return;
			}
			
			// For TTOperators also check whether the index order is transposed
			if (isOperator) {
				bool transposed = false;
				
				auto midIndexItr = myIndices.begin();
				size_t spanSum = 0;
				while (spanSum < numComponents) {
					REQUIRE(midIndexItr != myIndices.end(), "Internal Error.");
					spanSum += midIndexItr->span;
					++midIndexItr;
				}
				if (spanSum == numComponents) {
					// Transposition possible on my end
					auto otherMidIndexItr = otherIndices.begin();
					spanSum = 0;
					while (spanSum < numComponents) {
						REQUIRE(otherMidIndexItr != otherIndices.end(), "Internal Error.");
						spanSum += otherMidIndexItr->span;
						++otherMidIndexItr;
					}
					if (spanSum == numComponents) {
						// Other tensor also transposable
						transposed = (misc::equal(myIndices.begin(), midIndexItr, otherMidIndexItr, otherIndices.end())) 
									&& (misc::equal(midIndexItr, myIndices.end(), otherIndices.begin(), otherMidIndexItr));
					}
				}
				
				if (transposed) {
					if (otherTTN) {
						*meTTN = *otherTTN;
					} else {
						static_cast<TensorNetwork*>(_me.tensorObject)->operator=(*_other.tensorObjectReadOnly);
						if (otherTTStack->cannonicalization_required) {
							meTTN->move_core(otherTTStack->futureCorePosition);
						} else {
							meTTN->cannonicalized = false;
						}
					}
					REQUIRE(is_valid_tt(), "Internal Error.");
					dynamic_cast<TTOperator*>(_me.tensorObject)->transpose(); // NOTE will never be called if !isOperator
					return;
				}
			}
		}
		 
		// Use FullTensor fallback
		CHECK(_other.tensorObjectReadOnly->nodes.size() <= 1, warning, "Assigning a general tensor network to TTOperator not yet implemented. casting to fullTensor first");
		std::unique_ptr<Tensor> otherFull(_other.tensorObjectReadOnly->fully_contracted_tensor());
		std::unique_ptr<Tensor> otherReordered(otherFull->construct_new());
		(*otherReordered)(myIndices) = (*otherFull)(otherIndices);
		
		// Cast to TTNetwork
		if(otherReordered->is_sparse()) {
			LOG(fatal, "Not yet implemented." ); //TODO
		} else {
			*_me.tensorObject = TTNetwork(std::move(*static_cast<FullTensor*>(otherReordered.get())));
		}
	}
	
	// Explicit instantiation of the two template parameters that will be implemented in the xerus library
	template class TTNetwork<false>;
	template class TTNetwork<true>;
}
