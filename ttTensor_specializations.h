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

#include "xerus.h"

namespace xerus {

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

// forward declaration of the two template instanciations that will be implemented in the xerus library
extern template class TTNetwork<false>;
extern template class TTNetwork<true>;


}


