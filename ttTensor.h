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

#include "contractionHeuristic.h"

namespace xerus {

template<bool isOperator>
class TTNetwork : public TensorNetwork {
protected:
	static void construct_train_from_full(TensorNetwork& _out, const FullTensor& _A, const double _eps) {
        const size_t N = isOperator?2:1;
        REQUIRE(_eps < 1, "_eps must be smaller than one. " << _eps << " was given.");
        REQUIRE(_A.degree()%N==0, "Number of indicis must be even for TTOperator");
        
        //The external dimensions are the same as for the FullTensor
        _out.dimensions = _A.dimensions;
        
        // Number of Nodes to create
        const size_t numNodes = _A.degree()/N;
        
        // If _A has degree zero we are finished
        if(_A.degree() == 0) { 
			_out.factor = _A.data.get()[0];
			return; 
		}
		_out.factor = 1.0;
        
        // If there is only one node in the resulting train object we are already finished
        if(_A.degree() == N) {
            _out.nodes.emplace_back(std::shared_ptr<Tensor>(new FullTensor(_A)));
            _out.externalLinks.emplace_back(0, 0, _A.dimensions[0], false);
            _out.nodes.back().neighbors.emplace_back(0, 0, _A.dimensions[0], true);
            if(N==2) {
                _out.externalLinks.emplace_back(0, 1, _A.dimensions[1], false);
                _out.nodes.back().neighbors.emplace_back(0, 1, _A.dimensions[1], true);
            }
            return;
        }
        
        // Create the externalLinks first, as we know their position in advance
        _out.externalLinks.emplace_back(0, 0, _A.dimensions[0], false);
        for(size_t i = 1; i < numNodes; ++i) {
            _out.externalLinks.emplace_back(i, 1, _A.dimensions[i], false);
        }
        if(N == 2) {
            _out.externalLinks.emplace_back(0, 1, _A.dimensions[numNodes], false);
            for(size_t i = 1; i < numNodes; ++i) {
                _out.externalLinks.emplace_back(i, 2, _A.dimensions[numNodes+i], false);
            }
        }
        
        // Needed variables
        std::shared_ptr<FullTensor> nxtTensor;
        std::unique_ptr<value_t[]> currentU, currentS, currentVt, workingData, oldS;
        size_t leftDim, remainingDim=_A.size, maxRank, newRank, oldRank=1;
        
        // If we want a TTOperator we need to reshuffle the indices first, otherwise we want to copy the data because Lapack wants to destroy it
        if(N==1) {
            workingData.reset(new value_t[_A.size]);
            array_copy(workingData.get(), _A.data.get(), _A.size);
        } else {
            FullTensor tmpTensor(_A.degree());
            std::vector<Index> presentIndices, newIndices;
            for(size_t i = 0; i < _A.degree(); ++i) { presentIndices.emplace_back(); }
            for(size_t i = 0; i < numNodes; ++i) {
                newIndices.emplace_back(presentIndices[i]);
                newIndices.emplace_back(presentIndices[i+numNodes]);
            }
            tmpTensor(newIndices) = _A(presentIndices);
            workingData.reset(new value_t[tmpTensor.size]); // TODO unnecessary copy
            array_copy(workingData.get(), tmpTensor.data.get(), tmpTensor.size);
        }
        
        // This was the first step. Now automated steps until we reach the end
        for(size_t position = 0; position < numNodes-1; ++position) {
            // Determine Dimensions of the next matrification
            leftDim = oldRank*_A.dimensions[position];
            if(N == 2) { leftDim *= _A.dimensions[position+numNodes]; }
            remainingDim /= _A.dimensions[position];
            if(N == 2) { remainingDim /= _A.dimensions[position+numNodes]; }
            maxRank = std::min(leftDim, remainingDim);
            
            // Create temporary space for the results
            currentU.reset(new value_t[leftDim*maxRank]);
            currentS.reset(new value_t[maxRank]);
            currentVt.reset(new value_t[maxRank*remainingDim]);
            
            // Actually calculate the SVD --  We may destroy the current workingData
            blasWrapper::svd_destructive(currentU.get(), currentS.get(), currentVt.get(), workingData.get(), leftDim, remainingDim);
            
            //Determine the Rank
            for(newRank = maxRank; currentS[newRank-1] < _eps*currentS[0]; --newRank) { }
            
            // Create a FullTensor for U
            std::vector<size_t> dimensions;
            if(position != 0) { dimensions.emplace_back(oldRank); }
            dimensions.emplace_back(_A.dimensions[position]);
            if(N == 2) { dimensions.emplace_back(_A.dimensions[position+numNodes]); }
            dimensions.emplace_back(newRank);
            if(newRank == maxRank) {
                nxtTensor.reset(new FullTensor(std::move(dimensions), std::move(currentU)) );
            } else {
                nxtTensor.reset(new FullTensor(std::move(dimensions), DONT_SET_ZERO()) );
                for(size_t i = 0; i < leftDim; ++i) {
                    array_copy(nxtTensor->data.get()+i*newRank, currentU.get()+i*maxRank, newRank);
                }
            }
            
            // Create a node for U
            _out.nodes.emplace_back(std::static_pointer_cast<Tensor>(nxtTensor));
            if(position > 0) { _out.nodes.back().neighbors.emplace_back(_out.nodes.size()-2, ((position == 1) ? 0:1)+N, oldRank, false); }
            _out.nodes.back().neighbors.emplace_back(-1, position, _A.dimensions[position], true);
            if(N == 2) { _out.nodes.back().neighbors.emplace_back(-1, position+numNodes, _A.dimensions[position+numNodes], true); }
            _out.nodes.back().neighbors.emplace_back(_out.nodes.size(), 0, newRank, false);
                
            // Calclate S*Vt by scaling the rows of Vt
            for(size_t row = 0; row < newRank; ++row) {
                array_scale(currentVt.get()+row*remainingDim, currentS[row], remainingDim);
            }
            
            // Save the Stuff still needed
            workingData = std::move(currentVt);
            oldS = std::move(currentS);
            oldRank = newRank;
        }
        
        // Create FullTensor for Vt
        if(N==1) { 
			nxtTensor.reset(new FullTensor({oldRank, _A.dimensions[_A.degree()-1],                             }, DONT_SET_ZERO()) ); 
		} else { 
			nxtTensor.reset(new FullTensor({oldRank, _A.dimensions[_A.degree()/N-1], _A.dimensions[_A.degree()-1]}, DONT_SET_ZERO()) ); 
		}
        array_copy(nxtTensor->data.get(), workingData.get(), oldRank*remainingDim);
        
        // Create final Node for Vt
        _out.nodes.emplace_back(std::static_pointer_cast<Tensor>(nxtTensor));
        _out.nodes.back().neighbors.emplace_back(_out.nodes.size()-2, (numNodes == 2 ? N : 1+N), oldRank, false);
        _out.nodes.back().neighbors.emplace_back(-1, _A.degree()/N-1, _A.dimensions[_A.degree()/N-1], true);
        if(N == 2) { _out.nodes.back().neighbors.emplace_back(-1, _A.degree()-1, _A.dimensions[_A.degree()-1], true); }
        
        REQUIRE((N==1 && remainingDim == _A.dimensions.back()) || (N==2 && remainingDim == _A.dimensions[_A.degree()/2-1]*_A.dimensions[_A.degree()-1]), "Internal Error");
		REQUIRE(_out.is_in_expected_format(), "ie");
    }
    
    static void round_train(TensorNetwork& _me, const std::vector<size_t>& _maxRanks, const double _eps) {
        const size_t N = isOperator?2:1;
        REQUIRE(_me.degree()%N==0, "Number of indicis must be even for TTOperator");
        REQUIRE(_eps < 1, "_eps must be smaller than one. " << _eps << " was given.");
        REQUIRE(_maxRanks.size() == _me.degree()/N-1, "There must be exactly degree/N-1 maxRanks. Here " << _maxRanks.size() << " instead of " << _me.degree()/N-1 << " are given.");
        
        // If there is no or only one node in the train object we are already finished as there is no rank that can be rounded
        if(_me.degree() <= N) { return; }
        
        // Needed variables
        std::unique_ptr<value_t[]> LR, RR, M, U, S, Vt, newLeft, newRight;
        size_t leftDim, midDim, rightDim, maxLeftRank, maxRightRank, maxRank, rank;
        
        for(size_t position = 0; position < _me.degree()/N-1; ++position) {
            REQUIRE(dynamic_cast<FullTensor*>(_me.nodes[position].tensorObject.get()), "Tensor Train nodes are required to be FullTensors");
            REQUIRE(dynamic_cast<FullTensor*>(_me.nodes[position+1].tensorObject.get()), "Tensor Train nodes are required to be FullTensors");
            
            FullTensor& leftTensor = *static_cast<FullTensor*>(_me.nodes[position].tensorObject.get());
            FullTensor& rightTensor = *static_cast<FullTensor*>(_me.nodes[position+1].tensorObject.get());
            
            // Determine the dimensions of the next matrifications
            REQUIRE(leftTensor.dimensions.back() == rightTensor.dimensions.front(), "Internal Error");
            midDim   = leftTensor.dimensions.back();
            leftDim  = leftTensor.size/midDim;
            rightDim = rightTensor.size/midDim;
            maxLeftRank = std::min(leftDim, midDim);
            maxRightRank = std::min(midDim, rightDim);
            maxRank = std::min(maxLeftRank, maxRightRank);
            
            // Calculate QR and RQ decompositoins
            LR.reset(new value_t[maxLeftRank*midDim]);
            RR.reset(new value_t[midDim*maxRightRank]);
            blasWrapper::inplace_qr(leftTensor.data.get(), LR.get(), leftDim, midDim);
            blasWrapper::inplace_rq(RR.get(), rightTensor.data.get(), midDim, rightDim);
            
            // Calculate Middle Matrix M = LR*RR
            M.reset(new value_t[maxLeftRank*maxRightRank]);
            blasWrapper::matrix_matrix_product(M.get(), maxLeftRank, maxRightRank, 1.0, LR.get(), false, midDim, RR.get(), false);
            
            // Calculate SVD of U S Vt = M -- Reuse the space allocted for LR and RR and allow destruction of M
            U = std::move(LR);
            S.reset(new value_t[maxRank]);
            Vt = std::move(RR);
            blasWrapper::svd_destructive(U.get(), S.get(), Vt.get(), M.get(), maxLeftRank, maxRightRank);
            
            //Determine the Rank
            for(rank = maxRank; S[rank-1] < _eps*S[0]; --rank) { }
            rank = std::min(rank, _maxRanks[position]);
            
            // Calclate S*Vt by scaling the rows of Vt appropriatly
            for(size_t row = 0; row < rank; ++row) {
                array_scale(Vt.get()+row*maxRightRank, S[row], maxRightRank);
            }
            
            //Multiply U and (S*Vt) back to the left and to the right
            newLeft.reset(new value_t[leftDim*rank]);
            newRight.reset(new value_t[rank*rightDim]);
            blasWrapper::matrix_matrix_product(newLeft.get(), leftDim, rank, 1.0, leftTensor.data.get(), maxLeftRank, false, maxLeftRank, U.get(), maxRank, false);
            blasWrapper::matrix_matrix_product(newRight.get(), rank, rightDim, 1.0, Vt.get(), false, maxRightRank, rightTensor.data.get(), false);
            
            // Put the new Tensors to their place
            leftTensor.data.reset(newLeft.release(), internal::array_deleter_vt);
            rightTensor.data.reset(newRight.release(), internal::array_deleter_vt);
            leftTensor.dimensions.back() = rank;
            leftTensor.size = product(leftTensor.dimensions);
            rightTensor.dimensions.front() = rank;
            rightTensor.size = product(rightTensor.dimensions);
            _me.nodes[position  ].neighbors.back().dimension  = rank;
            _me.nodes[position+1].neighbors.front().dimension = rank;
        }
        REQUIRE(_me.is_in_expected_format(), "ie");
    }
    
    static void contract_stack(const IndexedTensorWritable<TensorNetwork> &_me) {
		REQUIRE(_me.tensorObject->is_valid_network(), "cannot contract inconsistent ttStack");
		const size_t N = isOperator?2:1;
		const size_t numNodes = _me.degree()/N;
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
		REQUIRE(_me.tensorObject->is_valid_network(), "something went wrong in contract_stack");
		
		// reset to new external links
		_me.tensorObject->externalLinks.clear();
		_me.tensorObject->externalLinks.emplace_back(0, 0, _me.tensorObject->dimensions[0], false);
		for(size_t i = 1; i < numNodes; ++i) {
			_me.tensorObject->externalLinks.emplace_back(i, 1, _me.tensorObject->dimensions[i], false);
		}
		if(N == 2) {
			_me.tensorObject->externalLinks.emplace_back(0, 1, _me.tensorObject->dimensions[numNodes], false);
			for(size_t i = 1; i < numNodes; ++i) {
				_me.tensorObject->externalLinks.emplace_back(i, 2, _me.tensorObject->dimensions[numNodes+i], false);
			}
		}
		
		// ensure right amount and order of links
		Index ext[N];
		size_t lastRank, externalDim[N], newRank;
		std::vector<Index> lastIndices, lastRight;
		std::vector<Index> oldIndices, newRight; // newLeft == lastRight
		std::vector<Index> newIndices;
		std::vector<size_t> newDimensions;
		for (size_t i=0; i<numNodes; ++i) {
			lastIndices = std::move(oldIndices); oldIndices.clear();
			lastRight = std::move(newRight); newRight.clear();
			lastRank = newRank; newRank=1;
			TensorNode &n = _me.tensorObject->nodes[i];
			for (TensorNode::Link &l : n.neighbors) {
				if (l.external) {
					size_t externalNumber = 0;
					if (N==2) {
						externalNumber = l.indexPosition>=numNodes?1:0;
					}
					oldIndices.push_back(ext[externalNumber]);
					externalDim[externalNumber] = l.dimension;
				} else if (i >= 1 && l.links(i-1)) {
					REQUIRE(lastIndices.size() > l.indexPosition, "ie " << lastIndices.size() << " " << l.indexPosition);
					oldIndices.push_back(lastIndices[l.indexPosition]);
				} else if (l.links(i+1)) {
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
			if (i>0) {
				n.neighbors.emplace_back(i-1,(i>1? N+1 : N),lastRank, false);
				newDimensions.push_back(lastRank);
			}
			for (size_t j=0; j<N; ++j) {
				REQUIRE(_me.tensorObject->dimensions[i+j*numNodes] == externalDim[j], "ie");
				n.neighbors.emplace_back(0,i+j*numNodes,externalDim[j], true);
				newDimensions.push_back(externalDim[j]);
			}
			if (i<numNodes-1) {
				n.neighbors.emplace_back(i+1,0,newRank, false);
				newDimensions.push_back(newRank);
			}
			n.tensorObject->reinterpret_dimensions(newDimensions);
		}
		
		
		REQUIRE(_me.tensorObject->is_in_expected_format(), "something went wrong in contract_stack");
	}
    
public:
    explicit TTNetwork() : TensorNetwork() {}
    
    explicit TTNetwork(size_t _degree) : TensorNetwork() {
		const size_t N=isOperator?2:1;
		REQUIRE(_degree%N==0, "illegal degree for TTOperator");
		const size_t numNodes = _degree/N;
		dimensions = std::vector<size_t>(_degree, 1);
		factor = 0.0;
		
		if (numNodes == 0) {
			return;
		}
		
		// externalLinks
		externalLinks.emplace_back(0, 0, 1, false);
		for(size_t i = 1; i < numNodes; ++i) {
			externalLinks.emplace_back(i, 1, 1, false);
		}
		if(N == 2) {
			externalLinks.emplace_back(0, 1, 1, false);
			for(size_t i = 1; i < numNodes; ++i) {
				externalLinks.emplace_back(i, 2, 1, false);
			}
		}
		
		REQUIRE(externalLinks.size() == _degree, "ie");
		
		std::vector<TensorNode::Link> neighbors;
		for (size_t i=0; i<numNodes; ++i) {
			neighbors.clear();
			if (i!=0) {
				neighbors.emplace_back(i-1, N+(i>=2?1:0), 1, false);
			}
			for (size_t j=0; j< N; ++j) { 
				neighbors.emplace_back(-1, i+j*numNodes, 1, true);
			}
			if (i!=numNodes-1) {
				neighbors.emplace_back(i+1, 0, 1, false);
			}
			
			nodes.emplace_back(
				std::shared_ptr<Tensor>(new FullTensor(neighbors.size())), 
				std::move(neighbors)
			);
		}
		REQUIRE(nodes.size() == _degree/N,"ie");
		REQUIRE(is_valid_tt(), "ie");
	} 
    
    explicit TTNetwork(const FullTensor& _full, const double _eps=1e-15) { //TODO no magic numbers
		construct_train_from_full(*this, _full, _eps);
	}
    
	ALLOW_MOVE(TTNetwork, TTT)
	explicit TTNetwork(TTT &&_cpy) : TensorNetwork(std::forward<TTT>(_cpy)) {
		// nothing to be done here
	}
	
	explicit TTNetwork(const TensorNetwork &_cpy, double _eps=1e-14) : TensorNetwork(_cpy) {
		LOG(fatal, "cast of arbitrary tensor network to TT not yet supported");
	}
	
	explicit TTNetwork(TensorNetwork &&_mov, double _eps=1e-14) : TensorNetwork(std::move(_mov)) {
		LOG(fatal, "cast of arbitrary tensor network to TT not yet supported");
	}
	
	template<class generator, class distribution>
	static TTNetwork construct_random(const std::vector<size_t>& _dimensions, const std::vector<size_t> &_ranks, generator& _rnd, distribution& _dist) {
        const size_t N = isOperator?2:1;
		REQUIRE(_ranks.size() == _dimensions.size()/N-1,"non-matching amount of ranks given to TTNetwork::construct_random");
		#ifndef DISABLE_RUNTIME_CHECKS_
		for (size_t d : _dimensions) {
			REQUIRE(d > 0, "trying to construct random TTTensor with dimension 0");
		}
		for (size_t d : _ranks) {
			REQUIRE(d > 0, "trying to construct random TTTensor with rank 0");
		}
		#endif
        
        TTNetwork result;
        if(_dimensions.size() == 0) {
// 			std::shared_ptr<Tensor> newT(new FullTensor(0));
// 			newT->data.get()[0] = dis(rnd);
// 			result.nodes.emplace_back(std::move(newT));
			result.factor = _dist(_rnd); //TODO create testcase for degree-0 TTtensors
			result.dimensions.clear();
        } else {
            result.nodes.clear();
            result.dimensions = _dimensions;
			result.factor = 1.0;
			
			// externalLinks
			size_t numNodes = _dimensions.size()/N;
			result.externalLinks.emplace_back(0, 0, result.dimensions[0], false);
			for(size_t i = 1; i < numNodes; ++i) {
				result.externalLinks.emplace_back(i, 1, result.dimensions[i], false);
			}
			if(N == 2) {
				result.externalLinks.emplace_back(0, 1, result.dimensions[numNodes], false);
				for(size_t i = 1; i < numNodes; ++i) {
					result.externalLinks.emplace_back(i, 2, result.dimensions[numNodes+i], false);
				}
			}
			
            size_t maxDim1 = product(_dimensions);
            size_t maxDim2 = 1;
			std::vector<size_t> constructVector;
			std::vector<TensorNode::Link> neighbors;
            for (size_t i=0; i<_dimensions.size()/N; ++i) {
                constructVector.clear();
				neighbors.clear();
				if (i!=0) {
					size_t r = std::min(_ranks[i-1],std::min(maxDim1, maxDim2));
					constructVector.emplace_back(r);
					neighbors.emplace_back(i-1, N+(i>=2?1:0), r, false);
				}
                for (size_t j=0; j< N; ++j) { 
					size_t r = _dimensions[i+j*numNodes];
					constructVector.emplace_back(r); 
					neighbors.emplace_back(-1, i+j*numNodes, r, true);
				}
                if (i!=_dimensions.size()/N-1) {
					size_t r = std::min(_ranks[i],std::min(maxDim1/_dimensions[i], maxDim2*_dimensions[i]));
					constructVector.emplace_back(r);
					neighbors.emplace_back(i+1, 0, r, false);
				}
				
				result.nodes.emplace_back(
					std::shared_ptr<Tensor>(new FullTensor(FullTensor::construct_random(constructVector, _rnd, _dist))), 
					std::move(neighbors)
				);
				
                maxDim1 /= _dimensions[i];
                maxDim2 *= _dimensions[i];
            }
        }
        REQUIRE(result.nodes.size() == _dimensions.size()/N,"ie");
		result.cannonicalize_right();
		REQUIRE(result.is_valid_tt(), "ie");
        return result;
    }
    
    template<class generator, class distribution>
	_inline_ static TTNetwork construct_random(const std::vector<size_t>& _dimensions, size_t _rank, generator& _rnd, distribution& _dist) {
		const size_t N = isOperator?2:1;
		return construct_random(_dimensions, std::vector<size_t>(_dimensions.size()/N-1, _rank), _rnd, _dist);
	}
    
    static TTNetwork construct_identity(const std::vector<size_t>& _dimensions) {
        const size_t N = isOperator?2:1;
		REQUIRE(isOperator, "tttensor identity ill-defined");
		#ifndef DISABLE_RUNTIME_CHECKS_
		for (size_t d : _dimensions) {
			REQUIRE(d > 0, "trying to construct random TTTensor with dimension 0");
		}
		#endif
        
        TTNetwork result;
        if(_dimensions.size() == 0) {
// 			std::shared_ptr<Tensor> newT(new FullTensor(0));
// 			newT->data.get()[0] = dis(rnd);
// 			result.nodes.emplace_back(std::move(newT));
			result.factor = 1; 
			result.dimensions.clear();
        } else {
            result.nodes.clear();
            result.dimensions = _dimensions;
			result.factor = 1.0;
			
			// externalLinks
			size_t numNodes = _dimensions.size()/N;
			result.externalLinks.emplace_back(0, 0, result.dimensions[0], false);
			for(size_t i = 1; i < numNodes; ++i) {
				result.externalLinks.emplace_back(i, 1, result.dimensions[i], false);
			}
			if(N == 2) {
				result.externalLinks.emplace_back(0, 1, result.dimensions[numNodes], false);
				for(size_t i = 1; i < numNodes; ++i) {
					result.externalLinks.emplace_back(i, 2, result.dimensions[numNodes+i], false);
				}
			}
			
			std::vector<size_t> constructVector1;
			std::vector<size_t> constructVector2;
			std::vector<TensorNode::Link> neighbors;
            for (size_t i=0; i<_dimensions.size()/N; ++i) {
                constructVector1.clear();
				constructVector2.clear();
				neighbors.clear();
				if (i!=0) {
					size_t r = 1;
					constructVector1.emplace_back(r);
					neighbors.emplace_back(i-1, N+(i>=2?1:0), r, false);
				}
                for (size_t j=0; j< N; ++j) { 
					size_t r = _dimensions[i+j*numNodes];
					constructVector1.emplace_back(r);
					constructVector2.emplace_back(r);
					neighbors.emplace_back(-1, i+j*numNodes, r, true);
				}
                if (i!=_dimensions.size()/N-1) {
					size_t r = 1;
					constructVector1.emplace_back(r);
					neighbors.emplace_back(i+1, 0, r, false);
				}
				
				// construct identity
				std::shared_ptr<Tensor> tI(new FullTensor(constructVector2, [](const std::vector<size_t> &_idx){
						if (_idx[0] == _idx[1]) {
							return 1.0;
						} else {
							return 0.0;
						}
					}));
				tI->reinterpret_dimensions(constructVector1);
				result.nodes.emplace_back(
					std::move(tI), 
					std::move(neighbors)
				);
            }
        }
        REQUIRE(result.nodes.size() == _dimensions.size()/N,"ie");
		result.cannonicalize_right();
		REQUIRE(result.is_valid_tt(), "ie");
        return result;
    }
    
    static TTNetwork dyadic_product(const TTNetwork &_lhs, const TTNetwork &_rhs) {
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
		const size_t lhsNodesSize = _lhs.nodes.size();
		const size_t rhsNodesSize = _rhs.nodes.size();
		for (const TensorNode &n : _rhs.nodes) {
			result.nodes.emplace_back(n);
			for (TensorNode::Link &l : result.nodes.back().neighbors) {
				if (l.external) {
					if (l.indexPosition < rhsNodesSize) {
						l.indexPosition += lhsNodesSize;
					} else {
						l.indexPosition += 2*lhsNodesSize;
					}
				} else {
					l.other += lhsNodesSize;
				}
			}
		}
		// add rank-1 connection between the two
		std::vector<size_t> dimensions(result.nodes[lhsNodesSize-1].tensorObject->dimensions);
		dimensions.push_back(1);
		result.nodes[lhsNodesSize-1].tensorObject->reinterpret_dimensions(dimensions);
		result.nodes[lhsNodesSize-1].neighbors.emplace_back(lhsNodesSize, 0, 1, false);
		const size_t linkPos = dimensions.size()-1;
		
		dimensions.clear();
		dimensions.push_back(1);
		dimensions.insert(dimensions.begin(), result.nodes[lhsNodesSize].tensorObject->dimensions.begin(), result.nodes[lhsNodesSize].tensorObject->dimensions.end());
		result.nodes[lhsNodesSize].tensorObject->reinterpret_dimensions(dimensions);
		std::vector<TensorNode::Link> newLinks;
		newLinks.emplace_back(lhsNodesSize-1, linkPos, 1, false);
		newLinks.insert(newLinks.begin(), result.nodes[lhsNodesSize].neighbors.begin(), result.nodes[lhsNodesSize].neighbors.end());
		result.nodes[lhsNodesSize].neighbors = newLinks;
		
		// add all external indices of rhs
		newLinks.clear();
		for (size_t i=0; i<lhsNodesSize; ++i) {
			newLinks.push_back(_lhs.externalLinks[i]);
		}
		for (size_t i=0; i<rhsNodesSize; ++i) {
			newLinks.push_back(_rhs.externalLinks[i]);
		}
		if (isOperator) {
			for (size_t i=0; i<lhsNodesSize; ++i) {
				newLinks.push_back(_lhs.externalLinks[lhsNodesSize+i]);
			}
			for (size_t i=0; i<rhsNodesSize; ++i) {
				newLinks.push_back(_rhs.externalLinks[rhsNodesSize+i]);
			}
		}
		result.externalLinks = newLinks;
		
		REQUIRE(result.is_valid_tt(), "ie");
		return result;
	}
	
	static TTNetwork dyadic_product(const std::vector<std::reference_wrapper<TTNetwork>> &_tensors) {
		if (_tensors.size() == 0) {
			return TTNetwork();
		} 
		TTNetwork result(_tensors.front());
		for (size_t i=0; i<_tensors.size(); ++i) {
			dyadic_product(result, _tensors[i]);
		}
		return result;
	}
	
	void round(value_t _eps) {
		cannonicalize_left();
		const size_t N = isOperator?2:1;
		round_train(*this, std::vector<size_t>(degree()/N-1, size_t(-1)), _eps);
	}

	void round(size_t _maxRank) {
		cannonicalize_left();
		const size_t N = isOperator?2:1;
		round_train(*this, std::vector<size_t>(degree()/N-1 ,_maxRank), 1e-15);
	}

	void round(const std::vector<size_t> &_maxRanks) {
		cannonicalize_left();
		round_train(*this, _maxRanks, 1e-15);
	}
    
    _inline_ void round(int _maxRank) {
        REQUIRE( _maxRank > 0, "MaxRank must be positive");
        round(size_t(_maxRank));
    }
	
	
	std::vector<size_t> ranks() const {
		std::vector<size_t> res;
		for (size_t n=0; n<nodes.size()-1; ++n) {
			res.push_back(nodes[n].neighbors.back().dimension);
		}
		return res;
	}
	
	size_t rank(size_t _i) const {
		REQUIRE(_i < nodes.size()-1, "requested illegal rank");
		return nodes[_i].neighbors.back().dimension;
	}
	
	size_t datasize() const {
		size_t result = 0;
		for (const TensorNode &n : nodes) {
			result += n.tensorObject->size;
		}
		return result;
	}
	
#ifndef DISABLE_RUNTIME_CHECKS_
	/// tests whether the network resembles that of a TTTensor and checks consistency with the udnerlying tensor objects
	/// @note will not check for orthogonality
	bool is_valid_tt() const {
		const size_t N = isOperator?2:1;
		const size_t numNodes = degree()/N;
		REQUIRE(nodes.size() == numNodes, nodes.size() << " vs " << numNodes);
		REQUIRE(externalLinks.size() == degree(), externalLinks.size() << " vs " << degree());
		REQUIRE(std::isfinite(factor), factor);
		
		// per external link
		for (size_t n=0; n<externalLinks.size(); ++n) {
			const TensorNode::Link &l = externalLinks[n];
			REQUIRE(l.dimension == dimensions[n], "n=" << n << " " << l.dimension << " vs " << dimensions[n]);
			REQUIRE(!l.external, "n=" << n);
			REQUIRE(l.other < numNodes, "n=" << n << " " << l.other << " vs " << numNodes);
			REQUIRE(l.indexPosition < nodes[l.other].neighbors.size(), "n=" << n << " " << l.indexPosition << " vs " << nodes[l.other].neighbors.size());
			REQUIRE(nodes[l.other].neighbors[l.indexPosition].external, "n=" << n);
			REQUIRE(nodes[l.other].neighbors[l.indexPosition].indexPosition == n, "n=" << n << " " << nodes[l.other].neighbors[l.indexPosition].indexPosition);
			REQUIRE(nodes[l.other].neighbors[l.indexPosition].dimension == l.dimension, "n=" << n << " " << nodes[l.other].neighbors[l.indexPosition].dimension << " vs " << l.dimension);
		}
		// per node
		for (size_t n=0; n<numNodes; ++n) {
			const TensorNode &node = nodes[n];
			REQUIRE(!node.erased, "n=" << n);
			if (node.tensorObject) {
				REQUIRE(n == numNodes-1 || !node.tensorObject->has_factor(), "n="<<n);
			}
			if (n==0) { // first node (or only node)
				REQUIRE(node.degree() == N+(numNodes>1?1:0), "n=" << n << " " << node.degree());
				if (node.tensorObject) {
					REQUIRE(node.tensorObject->degree() == N+(numNodes>1?1:0), "n=" << n << " " << node.tensorObject->degree());
				}
				REQUIRE(node.neighbors[0].external, "n=" << n);
				REQUIRE(node.neighbors[0].indexPosition == n, "n=" << n << " " << node.neighbors[0].indexPosition);
				if (isOperator) {
					REQUIRE(node.neighbors[1].external, "n=" << n);
					REQUIRE(node.neighbors[1].indexPosition == numNodes+n, "n=" << n << " " << node.neighbors[1].indexPosition << " vs " << numNodes+n);
				}
			} else {
				REQUIRE(node.degree() == N+(n<numNodes-1?2:1), "n=" << n << " " << node.degree());
				if (node.tensorObject) {
					REQUIRE(node.tensorObject->degree() == N+(n<numNodes-1?2:1), "n=" << n << " " << node.tensorObject->degree());
				}
				REQUIRE(!node.neighbors[0].external, "n=" << n);
				REQUIRE(node.neighbors[0].other == n-1, "n=" << n);
				REQUIRE(node.neighbors[0].indexPosition == N+(n>1?1:0), "n=" << n << " " << node.neighbors[0].indexPosition);
				REQUIRE(node.neighbors[1].external, "n=" << n);
				REQUIRE(node.neighbors[1].indexPosition == n, "n=" << n << " " << node.neighbors[0].indexPosition);
				if (isOperator) {
					REQUIRE(node.neighbors[2].external, "n=" << n);
					REQUIRE(node.neighbors[2].indexPosition == numNodes+n, "n=" << n << " " << node.neighbors[1].indexPosition << " vs " << numNodes+n);
				}
			}
			if (n < numNodes-1) {
				REQUIRE(!node.neighbors.back().external, "n=" << n);
				REQUIRE(node.neighbors.back().other == n+1, "n=" << n << " " << node.neighbors.back().other);
				REQUIRE(node.neighbors.back().indexPosition == 0, "n=" << n << " " << node.neighbors.back().indexPosition);
				REQUIRE(!nodes[n+1].neighbors.empty(), "n=" << n);
				REQUIRE(node.neighbors.back().dimension == nodes[n+1].neighbors[0].dimension, "n=" << n << " " << node.neighbors.back().dimension << " vs " << nodes[n+1].neighbors[0].dimension);
				
			} 
		}
		
		return true;
	}
#else
	/// performs no tests if checks are disabled
	bool is_valid_tt() const {
		return true;
	}
#endif
	
	/// moves core to the left (eg. to perform a round operation that then moves it to the right)
	void cannonicalize_left() {
		Index i,r,j;
		FullTensor core(2);
		for (size_t n=nodes.size()-1; n>0; --n) {
            REQUIRE(!nodes[n].erased, "ie n="<<n);
			Tensor &currTensor = *nodes[n].tensorObject;
			( core(j,r), currTensor(r,i&1) ) = RQ(currTensor(j,i&1));  //TODO we want a rank-detecting QR at this point?
			Tensor &nextTensor = *nodes[n-1].tensorObject;
			nextTensor(j&1,i) = nextTensor(j&1,r) * core(r,i);
			if (currTensor.dimensions[0] != nodes[n].neighbors.front().dimension) {
				nodes[n].neighbors.front().dimension = nodes[n-1].neighbors.back().dimension = currTensor.dimensions[0];
			}
		}
	}
	
	/// moves core to the right
	void cannonicalize_right() {
		Index i,r,j;
		FullTensor core(2);
		for (size_t n=0; n<nodes.size()-1; ++n) {
			Tensor &currTensor = *nodes[n].tensorObject;
			( currTensor(i&1,r), core(r,j) ) = QR(currTensor(i&1,j)); //TODO we want a rank-detecting QR at this point?
			Tensor &nextTensor = *nodes[n+1].tensorObject;
			nextTensor(i,j&1) = core(i,r) * nextTensor(r,j&1);
			if (nextTensor.dimensions[0] != nodes[n+1].neighbors.front().dimension) {
				nodes[n+1].neighbors.front().dimension = nodes[n].neighbors.back().dimension = nextTensor.dimensions[0];
			}
		}
	}
	    
    /// swaps all external indices to create the transposed operator
    template<bool B = isOperator, typename std::enable_if<B, int>::type = 0>
	void transpose() {
		Index i,r,l,j;
		for (size_t n=0; n<nodes.size(); ++n) {
			size_t ld = n==0?0:1;
			size_t rd = n==nodes.size()-1?0:1;
			REQUIRE(nodes[n].neighbors.size() == ld + rd + 2, "ie " << ld << " " << rd << " " << nodes[n].neighbors.size() << " n= " << n);
			(*nodes[n].tensorObject)(l^ld,j,i,r^rd) = (*nodes[n].tensorObject)(l^ld,i,j,r^rd);
			TensorNode::Link &l1 = nodes[n].neighbors[ld];
			TensorNode::Link &l2 = nodes[n].neighbors[ld+1];
			REQUIRE(l1.external, "ie "<< l1.other << " " << ld << " " << n << " " << nodes[n].neighbors.size() << " " << nodes.size());
			REQUIRE(l2.external, "ie "<< l2.other);
			std::swap(externalLinks[l1.indexPosition].dimension, externalLinks[l2.indexPosition].dimension);
			std::swap(l1.dimension, l2.dimension);
		}
	}
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - -  Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - */
	TTNetwork& operator+=(const TTNetwork& _other) {
		Index i;
		(*this)(i&0) = (*this)(i&0) + _other(i&0);
		return *this;
	}
	TTNetwork  operator+(const TTNetwork& _other) const {
		TTNetwork cpy(*this);
		cpy += _other;
		return cpy;
	}
	
	TTNetwork& operator-=(const TTNetwork& _other) {
		Index i;
		(*this)(i&0) = (*this)(i&0) - _other(i&0);
		return *this;
	}
	TTNetwork  operator-(const TTNetwork& _other) const {
		TTNetwork cpy(*this);
		cpy -= _other;
		return cpy;
	}
	
	TTNetwork& operator*=(const value_t _prod) {
		factor *= _prod;
		return *this;
	}
	TTNetwork  operator*(const value_t _prod) const {
		TTNetwork result(*this);
		result *= _prod;
		return result;
	}
	
	TTNetwork& operator/=(const value_t _div) {
		factor /= _div;
		return *this;
	}
	TTNetwork  operator/(const value_t _div) const {
		TTNetwork result(*this);
		result /= _div;
		return result;
	}
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
	virtual bool specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const override;
	virtual bool specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const override;
	virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) override;
	
	virtual TensorNetwork* get_copy() const override {
		return new TTNetwork(*this);
	}
	
	virtual value_t frob_norm() const override {
		REQUIRE(is_valid_tt(), "frob_norm of illegal TT");
		return nodes.back().tensorObject->frob_norm();
	}
	
	virtual bool is_in_expected_format() const override {
		return is_valid_tt();
	}
};

template<bool isOperator>
_inline_ TTNetwork<isOperator> operator*(const value_t _lhs, const TTNetwork<isOperator>& _rhs) { return _rhs*_lhs; }

/// Returns the frobenius norm of the given tensor
template<bool isOperator>
_inline_ value_t frob_norm(const TTNetwork<isOperator>& _tensor) { return _tensor.frob_norm(); }


typedef TTNetwork<false> TTTensor;
typedef TTNetwork<true> TTOperator;

namespace internal {

template<bool isOperator>
class TTStack : public TTNetwork<isOperator> {
public:
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
	virtual void specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me _unused_ , const IndexedTensorReadOnly<TensorNetwork> &_other _unused_) override {
		LOG(fatal, "TTStack not supported as a storing type");
	}
	virtual TensorNetwork* get_copy() const override {
		LOG(fatal, "forbidden");
		return nullptr;
	}
	virtual value_t frob_norm() const override {
		Index i;
		TTNetwork<isOperator> tmp(this->degree());
		tmp(i&0) = (*this)(i&0);
		return tmp.frob_norm();
	}
};



}
}


