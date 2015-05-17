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

#include <xerus/ttTensor.h>
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

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork() : TensorNetwork() {}

    
    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork(const size_t _degree) : TensorNetwork() {
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
    
    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork(const FullTensor& _full, const double _eps) { 
        construct_train_from_full(*this, _full, _eps);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - TTTensor - - - - - - - - - - - - - - - - - - - - - - - - - - */

    template<>
    bool TTNetwork<false>::specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out _unused_, const IndexedTensorReadOnly<TensorNetwork> &_me _unused_, const IndexedTensorReadOnly<TensorNetwork> &_other _unused_) const {
        return false; // Operator builds stacks, not the tttensor
    }

    template<bool isOperator>
    bool TTNetwork<isOperator>::specialized_sum(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const {
        const size_t N = isOperator?2:1;
        
        const std::vector<Index> myIndices = _me.get_assigned_indices();
        const std::vector<Index> otherIndices = _other.get_assigned_indices();
        
        // If the indices are in different order, we are lost. TODO inverse order is also ok...
        if(myIndices != otherIndices) { return false; }
        REQUIRE(_me.tensorObjectReadOnly->dimensions == _other.tensorObjectReadOnly->dimensions, "TT sum requires both opearants to share the same dimensions");
        
        // If the other is not a TT tensor we are also lost
        const TTTensor* _otherPtr = dynamic_cast<const TTTensor*>( _other.tensorObjectReadOnly);
        if(!_otherPtr) { return false; }
        
        // TODO the order is not canonical, because if I am no Stack I don't have to know whether or not i am moveable
        // If I am in fact a TTTensorStack, we have to evaluate me to TTTensor
        std::unique_ptr<IndexedTensor<TensorNetwork>> meStorage;
        const IndexedTensorReadOnly<TensorNetwork> *realMePtr = &_me;
        const IndexedTensorMoveable<TensorNetwork> *movMe = dynamic_cast<const IndexedTensorMoveable<TensorNetwork> *>(&_me);
        if (movMe) {
            internal::TTStack<isOperator> *stackMe = dynamic_cast<internal::TTStack<isOperator> *>(movMe->tensorObject);
            if (stackMe) {
                meStorage.reset(new IndexedTensor<TensorNetwork>(new TTTensor(_me.degree()), myIndices, true));
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
                otherStorage.reset(new IndexedTensor<TensorNetwork>(new TTTensor(_other.degree()), otherIndices, true));
                (*otherStorage) = _other;
                realOtherPtr = otherStorage.get();
            }
        } else {
            REQUIRE(!dynamic_cast<const internal::TTStack<isOperator> *>(_other.tensorObjectReadOnly),"ie - non-moveable TTStack detected");
        }
        const IndexedTensorReadOnly<TensorNetwork> &realOther = *realOtherPtr;
        
        // Number of Nodes to create
        const size_t numNodes = realMe.degree()/N;
        
        TTTensor* tmpPtr = new TTTensor();
        tmpPtr->factor = 1.0;
        
        //The external dimensions are the same as the ones of the input
        tmpPtr->dimensions = realMe.tensorObjectReadOnly->dimensions;
        REQUIRE(realOther.tensorObjectReadOnly->dimensions == realMe.tensorObjectReadOnly->dimensions, "Internal Error");
        
        IndexedTensor<TensorNetwork> tmpOut(tmpPtr, myIndices, true);
        TTTensor& outTensor = *static_cast<TTTensor*>(tmpOut.tensorObject);
        
        
        // Create the externalLinks first, as we know their position in advance
        outTensor.externalLinks.emplace_back(0, 0, outTensor.dimensions[0], false);
        for(size_t i = 1; i < numNodes; ++i) {
            outTensor.externalLinks.emplace_back(i, 1, outTensor.dimensions[i], false);
        }
        if(N == 2) {
            outTensor.externalLinks.emplace_back(0, 1, outTensor.dimensions[numNodes], false);
            for(size_t i = 1; i < numNodes; ++i) {
                outTensor.externalLinks.emplace_back(i, 2, outTensor.dimensions[numNodes+i], false);
            }
        }
        
        
        if(realMe.degree() == N) {
            // Create the one Node
            std::shared_ptr<Tensor> nextTensor;
            if(realMe.tensorObjectReadOnly->nodes[0].tensorObject->is_sparse() && realOther.tensorObjectReadOnly->nodes[0].tensorObject->is_sparse()) { // Both Sparse
                nextTensor.reset(realMe.tensorObjectReadOnly->nodes[0].tensorObject->get_copy());
                nextTensor->factor *= realMe.tensorObjectReadOnly->factor;
                *std::static_pointer_cast<SparseTensor>(nextTensor) += realOther.tensorObjectReadOnly->factor*(*std::static_pointer_cast<SparseTensor>(realOther.tensorObjectReadOnly->nodes[0].tensorObject));
            } else { // Maximal one sparse
                if(realMe.tensorObjectReadOnly->nodes[0].tensorObject->is_sparse()){
                    nextTensor.reset(new FullTensor(*static_cast<SparseTensor*>(realMe.tensorObjectReadOnly->nodes[0].tensorObject.get())));
                } else {
                    nextTensor.reset(new FullTensor(*static_cast<FullTensor*>(realMe.tensorObjectReadOnly->nodes[0].tensorObject.get())));
                }
                nextTensor->factor *= realMe.tensorObjectReadOnly->factor;
                if(realOther.tensorObjectReadOnly->nodes[0].tensorObject->is_sparse()){
                    *static_cast<FullTensor*>(nextTensor.get()) += realOther.tensorObjectReadOnly->factor*static_cast<SparseTensor&>(*realOther.tensorObjectReadOnly->nodes[0].tensorObject.get());
                } else {
                    *static_cast<FullTensor*>(nextTensor.get()) += realOther.tensorObjectReadOnly->factor*static_cast<FullTensor&>(*realOther.tensorObjectReadOnly->nodes[0].tensorObject.get());
                }
            }
            
            outTensor.nodes.emplace_back(std::static_pointer_cast<Tensor>(nextTensor));
            outTensor.nodes.back().neighbors.emplace_back(-1, 0, outTensor.dimensions[0], true);
            if(N == 2) { outTensor.nodes.back().neighbors.emplace_back(-1, 1, outTensor.dimensions[1], true); }
            _out.assign(std::move(tmpOut));
            return true;
        }
        
        for(size_t position = 0; position < numNodes; ++position) {
            // Get current input nodes
            // TODO sparse
            FullTensor &myNode = *static_cast<FullTensor*>(realMe.tensorObjectReadOnly->nodes[position].tensorObject.get());
            FullTensor &otherNode = *static_cast<FullTensor*>(realOther.tensorObjectReadOnly->nodes[position].tensorObject.get());
            
            // Structure has to be (for degree 4)
            // (L1 R1) * ( L2 0  ) * ( L3 0  ) * ( L4 )
            //           ( 0  R2 )   ( 0  R3 )   ( R4 )
            
            // Create a FullTensor for Node
            std::vector<size_t> nxtDimensions;
            if(position != 0) { 
                nxtDimensions.emplace_back(myNode.dimensions.front()+otherNode.dimensions.front());
            }
            nxtDimensions.emplace_back(outTensor.dimensions[position]);
            if(N == 2) { nxtDimensions.emplace_back(outTensor.dimensions[position+numNodes]); }
            if(position != numNodes-1) {
                nxtDimensions.emplace_back(myNode.dimensions.back()+otherNode.dimensions.back());
            }
            std::shared_ptr<FullTensor> nxtTensor(new FullTensor(std::move(nxtDimensions)) );
            
            
            // Create the Node
            outTensor.nodes.emplace_back(std::static_pointer_cast<Tensor>(nxtTensor));
            if(position != 0) { outTensor.nodes.back().neighbors.emplace_back(position-1, ((position == 1) ? 0:1)+N, nxtTensor->dimensions.front(), false); }
            outTensor.nodes.back().neighbors.emplace_back(-1, position, outTensor.dimensions[position], true);
            if(N == 2) { outTensor.nodes.back().neighbors.emplace_back(-1, position+numNodes, outTensor.dimensions[position+numNodes], true); }
            if(position != numNodes-1 ) { outTensor.nodes.back().neighbors.emplace_back(position+1, 0, nxtTensor->dimensions.back(), false); }

            const size_t leftIdxOffset = nxtTensor->size/nxtTensor->dimensions.front();
            const size_t extIdxOffset = nxtTensor->dimensions.back();
            const size_t myLeftIdxOffset = myNode.size/myNode.dimensions.front();
            const size_t myExtIdxOffset = myNode.dimensions.back();
            const size_t otherLeftIdxOffset = otherNode.size/otherNode.dimensions.front();
            const size_t otherExtIdxOffset = otherNode.dimensions.back();
            const size_t otherGeneralOffset = (position == 0 ? 0 : myNode.dimensions.front()*leftIdxOffset) + (position == numNodes-1 ? 0 : myNode.dimensions.back());
            
            
            
            // Copy own Tensor into place
            if(position == numNodes-1) {
                for(size_t leftIdx = 0; leftIdx < myNode.dimensions.front(); ++leftIdx) {
                    for(size_t extIdx = 0; extIdx < myNode.size/(myNode.dimensions.front()*myNode.dimensions.back()); ++extIdx) {
                        // RightIdx can be copied as one piece
                        array_scaled_copy(nxtTensor->data.get() + leftIdx*leftIdxOffset + extIdx*extIdxOffset, myNode.factor*realMe.tensorObjectReadOnly->factor, myNode.data.get() + leftIdx*myLeftIdxOffset + extIdx*myExtIdxOffset, myNode.dimensions.back());
                    }
                }
            } else {
                REQUIRE(!myNode.has_factor(), "Only Core node, which has to be the last node, is allowed to have a factor");
                for(size_t leftIdx = 0; leftIdx < myNode.dimensions.front(); ++leftIdx) {
                    for(size_t extIdx = 0; extIdx < myNode.size/(myNode.dimensions.front()*myNode.dimensions.back()); ++extIdx) {
                        // RightIdx can be copied as one piece
                        array_copy(nxtTensor->data.get() + leftIdx*leftIdxOffset + extIdx*extIdxOffset, myNode.data.get() + leftIdx*myLeftIdxOffset + extIdx*myExtIdxOffset, myNode.dimensions.back());
                    }
                }
            }
            
            
            // Copy other Tensor into place
            if(position == numNodes-1) {
                for(size_t leftIdx = 0; leftIdx < otherNode.dimensions.front(); ++leftIdx) {
                    for(size_t extIdx = 0; extIdx < otherNode.size/(otherNode.dimensions.front()*otherNode.dimensions.back()); ++extIdx) {
                        // RightIdx can be copied as one piece
                        array_scaled_copy(nxtTensor->data.get() + leftIdx*leftIdxOffset + extIdx*extIdxOffset + otherGeneralOffset, otherNode.factor*realOther.tensorObjectReadOnly->factor, otherNode.data.get() + leftIdx*otherLeftIdxOffset + extIdx*otherExtIdxOffset, otherNode.dimensions.back());
                    }
                }
            } else {
                REQUIRE(!otherNode.has_factor(), "Only Core node, which has to be the last node, is allowed to have a factor");
                for(size_t leftIdx = 0; leftIdx < otherNode.dimensions.front(); ++leftIdx) {
                    for(size_t extIdx = 0; extIdx < otherNode.size/(otherNode.dimensions.front()*otherNode.dimensions.back()); ++extIdx) {
                        // RightIdx can be copied as one piece
                        array_copy(nxtTensor->data.get() + leftIdx*leftIdxOffset + extIdx*extIdxOffset + otherGeneralOffset, otherNode.data.get() + leftIdx*otherLeftIdxOffset + extIdx*otherExtIdxOffset, otherNode.dimensions.back());
                    }
                }
            }
        }
        
        outTensor.cannonicalize_right();
        _out.assign(std::move(tmpOut));
        return true;
    }

    template<>
    void TTNetwork<false>::specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) {
        REQUIRE(_me.tensorObject == this, "Internal Error.");
        const std::vector<Index> myIndices = _me.get_assigned_indices();
        const std::vector<Index> otherIndices = _other.get_assigned_indices();
        
        if (myIndices == otherIndices) {
            // special case: _other is a TTStack with same indices
            const internal::TTStack<false> *otherTTs = dynamic_cast<const internal::TTStack<false> *>(_other.tensorObjectReadOnly);
            if (otherTTs) {
                *static_cast<TensorNetwork*>(this) = *_other.tensorObjectReadOnly; // TODO haesslich
                contract_stack(_me);
                REQUIRE(degree() == nodes.size(), "ie");
                
                dynamic_cast<TTTensor *>(_me.tensorObject)->cannonicalize_right();
                return;
            }
            
            // special case: _other is a TTTensor with same indices
            // NOTE as TTTensorStack inherits from TTTensor, the order is important
            const TTTensor *otherTT = dynamic_cast<const TTTensor *>(_other.tensorObjectReadOnly);
            if (otherTT) {
                *this = *otherTT; // assignment of tensorNetwork is sufficient
                return;
            }
        }
        
        // general case: not yet implemented... // TODO
        if (_other.tensorObjectReadOnly->nodes.size() > 1) {
            LOG(warning, "assigning a general tensor network to TTTensor not yet implemented. casting to fullTensor first");
        }
        // special case: _other is a fullTensor
        FullTensor rhs(_other.degree());
        // get _other with my index order
        rhs(myIndices) = _other; 
        // cast to TTTensor
        *_me.tensorObject = TTTensor(rhs);
    }
    
    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork(const TTNetwork & _cpy) : TensorNetwork(_cpy) { }
    
    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork(      TTNetwork&& _mov) : TensorNetwork(std::move(_mov)) { }

    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork(const TensorNetwork &_cpy, double _eps) : TensorNetwork(_cpy) {
        LOG(fatal, "cast of arbitrary tensor network to TT not yet supported");
    }

    template<bool isOperator>
    TTNetwork<isOperator>::TTNetwork(TensorNetwork &&_mov, double _eps) : TensorNetwork(std::move(_mov)) {
        LOG(fatal, "cast of arbitrary tensor network to TT not yet supported");
    }
    
    template<bool isOperator>
    TTNetwork<isOperator> TTNetwork<isOperator>::construct_identity(const std::vector<size_t>& _dimensions) {
        const size_t N = isOperator?2:1;
        REQUIRE(isOperator, "tttensor identity ill-defined");
        #ifndef DISABLE_RUNTIME_CHECKS_
        for (size_t d : _dimensions) {
            REQUIRE(d > 0, "trying to construct random TTTensor with dimension 0");
        }
        #endif
        
        TTNetwork result;
        if(_dimensions.size() == 0) {
//          std::shared_ptr<Tensor> newT(new FullTensor(0));
//          newT->data.get()[0] = dis(rnd);
//          result.nodes.emplace_back(std::move(newT));
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
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    template<bool isOperator>
    void TTNetwork<isOperator>::construct_train_from_full(TensorNetwork& _out, const FullTensor& _A, const double _eps) {
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
    
    
    template<bool isOperator>
    void TTNetwork<isOperator>::round_train(TensorNetwork& _me, const std::vector<size_t>& _maxRanks, const double _eps) {
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
    
    
    template<bool isOperator>
    void TTNetwork<isOperator>::contract_stack(const IndexedTensorWritable<TensorNetwork> &_me) {
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
    
    #ifndef DISABLE_RUNTIME_CHECKS_
        template<bool isOperator>
        bool TTNetwork<isOperator>::is_valid_tt() const {
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
                    if (node.tensorObject) {
                        REQUIRE(!node.tensorObject->has_factor(), "n="<<n);
                    }
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
        template<bool isOperator>
        bool TTNetwork<isOperator>::is_valid_tt() const {
            return true;
        }
    #endif
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    template<bool isOperator>
    TTNetwork<isOperator> TTNetwork<isOperator>::dyadic_product(const TTNetwork<isOperator> &_lhs, const TTNetwork<isOperator> &_rhs) {
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
    
    template<bool isOperator>
    TTNetwork<isOperator> TTNetwork<isOperator>::dyadic_product(const std::vector<std::reference_wrapper<TTNetwork<isOperator>>> &_tensors) {
        if (_tensors.size() == 0) {
            return TTNetwork();
        } 
        TTNetwork result(_tensors.front());
        for (size_t i=0; i<_tensors.size(); ++i) {
            dyadic_product(result, _tensors[i]);
        }
        return result;
    }
    
    
    template<bool isOperator>
    void TTNetwork<isOperator>::round(value_t _eps) {
        cannonicalize_left();
        const size_t N = isOperator?2:1;
        round_train(*this, std::vector<size_t>(degree()/N-1, size_t(-1)), _eps);
    }

    template<bool isOperator>
    void TTNetwork<isOperator>::round(size_t _maxRank) {
        cannonicalize_left();
        const size_t N = isOperator?2:1;
        round_train(*this, std::vector<size_t>(degree()/N-1 ,_maxRank), 1e-15);
    }

    template<bool isOperator>
    void TTNetwork<isOperator>::round(const std::vector<size_t> &_maxRanks) {
        cannonicalize_left();
        round_train(*this, _maxRanks, 1e-15);
    }
    
    template<bool isOperator>
    void TTNetwork<isOperator>::round(int _maxRank) {
        REQUIRE( _maxRank > 0, "MaxRank must be positive");
        round(size_t(_maxRank));
    }

    template<bool isOperator>
    std::vector<size_t> TTNetwork<isOperator>::ranks() const {
        std::vector<size_t> res;
        for (size_t n=0; n<nodes.size()-1; ++n) {
            res.push_back(nodes[n].neighbors.back().dimension);
        }
        return res;
    }
    
    template<bool isOperator>
    size_t TTNetwork<isOperator>::rank(size_t _i) const {
        REQUIRE(_i < nodes.size()-1, "requested illegal rank");
        return nodes[_i].neighbors.back().dimension;
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
    void TTNetwork<isOperator>::cannonicalize_left() {
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
    
    template<bool isOperator>
    void TTNetwork<isOperator>::cannonicalize_right() {
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
    
    template<bool isOperator>
    TensorNetwork* TTNetwork<isOperator>::get_copy() const {
        return new TTNetwork(*this);
    }
    
    template<bool isOperator>
    value_t TTNetwork<isOperator>::frob_norm() const {
        REQUIRE(is_valid_tt(), "frob_norm of illegal TT");
        return nodes.back().tensorObject->frob_norm();
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
    
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - TTOperator - - - - - - - - - - - - - - - - - - - - - - - - - - */

    template<>
    bool TTNetwork<true>::specialized_contraction(IndexedTensorWritable<TensorNetwork> &_out, const IndexedTensorReadOnly<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) const {
        REQUIRE(!_out.tensorObject, "Internal Error.");
        const std::vector<Index> myIndices = _me.get_assigned_indices(); // NOTE this is necessary to remove span-0 indices (and for REQUIREs)
        const std::vector<Index> otherIndices = _other.get_assigned_indices();
        
        bool otherTT = dynamic_cast<const TTTensor *>(_other.tensorObjectReadOnly) || dynamic_cast<const internal::TTStack<false> *>(_other.tensorObjectReadOnly);
        bool otherTO = !otherTT && (dynamic_cast<const TTOperator *>(_other.tensorObjectReadOnly) || dynamic_cast<const internal::TTStack<true> *>(_other.tensorObjectReadOnly));
        
        if (!otherTT && !otherTO) {
            return false;
        }
    // 	LOG(info, "If indices are right I can do a specialized contraction: " << myIndices << " vs. " << otherIndices);
        
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
            if (equal(myIndices.begin(), midIndexItr, otherIndices.begin(), otherIndices.end()) || equal(midIndexItr, myIndices.end(), otherIndices.begin(), otherIndices.end())) {
                TensorNetwork *res = new internal::TTStack<false>;
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
            if (   equal(myIndices.begin(), midIndexItr, otherIndices.begin(), otherMidIndexItr) 
                || equal(midIndexItr, myIndices.end(), otherIndices.begin(), otherMidIndexItr)
                || equal(myIndices.begin(), midIndexItr, otherMidIndexItr, otherIndices.end()) 
                || equal(midIndexItr, myIndices.end(), otherMidIndexItr, otherIndices.end())    ) 
            {
                TensorNetwork *res = new internal::TTStack<true>;
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

    template<>
    void TTNetwork<true>::specialized_evaluation(const IndexedTensorWritable<TensorNetwork> &_me, const IndexedTensorReadOnly<TensorNetwork> &_other) {
        const std::vector<Index> myIndices = _me.get_assigned_indices();
        const std::vector<Index> otherIndices = _other.get_assigned_indices();
        
        //TODO if TTTensor were a template class this would be so much simpler... TODO -> unify
        // transposed index order
        auto midIndexItr = myIndices.begin();
        size_t spanSum = 0;
        bool transposed = false;
        while (spanSum < _me.degree() / 2) {
            REQUIRE(midIndexItr != myIndices.end(), "ie");
            spanSum += midIndexItr->span;
            ++midIndexItr;
        }
        if (spanSum == _me.degree() / 2) {
            // tansposition possible on my end
            auto otherMidIndexItr = otherIndices.begin();
            spanSum = 0;
            while (spanSum < _other.degree() / 2) {
                REQUIRE(otherMidIndexItr != otherIndices.end(), "ie");
                spanSum += otherMidIndexItr->span;
                ++otherMidIndexItr;
            }
            if (spanSum == _other.degree() / 2) {
                // other tensor also transposable
                transposed = (equal(myIndices.begin(), midIndexItr, otherMidIndexItr, otherIndices.end())) 
                            && (equal(midIndexItr, myIndices.end(), otherIndices.begin(), otherMidIndexItr));
                
            }
        }
        
        // identical index order or transposed
        if (transposed || otherIndices == myIndices) {
            // special case: _other is a TTStack with same indices
            const internal::TTStack<true> *otherTTs = dynamic_cast<const internal::TTStack<true> *>(_other.tensorObjectReadOnly);
            if (otherTTs) {
                *_me.tensorObject = *_other.tensorObjectReadOnly;
                contract_stack(_me);
                REQUIRE(_me.tensorObject->nodes.size() == _me.degree()/2, "ie");
                dynamic_cast<TTOperator *>(_me.tensorObject)->cannonicalize_right();
                if (transposed) dynamic_cast<TTOperator *>(_me.tensorObject)->transpose();
                return;
            }
            
            // special case: _other is a TTOperator with same indices
            // NOTE as TTOperatorStack inherits from TTOperator, the order is important
            const TTOperator *otherTT = dynamic_cast<const TTOperator *>(_other.tensorObjectReadOnly);
            if (otherTT) {
                *_me.tensorObject = *otherTT; // assignment of tensorNetwork is sufficient
                if (transposed) dynamic_cast<TTOperator *>(_me.tensorObject)->transpose();
                return;
            }
        }
        
        // general case: not yet implemented... // TODO
        if (_other.tensorObjectReadOnly->nodes.size() > 1) {
            LOG(warning, "assigning a general tensor network to TTOperator not yet implemented. casting to fullTensor first");
        }
        // special case: _other is a fullTensor
        FullTensor rhs(_other.degree());
        // get _other with my index order
        rhs(myIndices) = _other; 
        // cast to TTOperator
        *_me.tensorObject = TTNetwork(rhs);
    }

    
    
    // Explicit instantiation of the two template parameters that will be implemented in the xerus library
    template class TTNetwork<false>;
    template class TTNetwork<true>;

}


