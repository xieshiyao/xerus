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
    
    std::unique_ptr<Tensor> prepare_split(size_t& _lhsSize, size_t& _rhsSize, size_t& _rank, const IndexedTensorReadOnly<Tensor>& _base, const IndexedTensorWritable<Tensor>& _lhs, const IndexedTensorWritable<Tensor>& _rhs) {
        const std::vector<Index> baseIndices = _base.get_assigned_indices();
        
        // Calculate future spans of lhs and rhs.
        size_t lhsSpan = 1; // Start with 1 because there is a new dimension introduced in the split. 
        for(const Index& idx : baseIndices) { if(contains(_lhs.indices, idx)) { lhsSpan += idx.span; } }
        const size_t rhsSpan = _base.degree() - lhsSpan + 2; // +1 because of the extra dim and +1 because the extra dim is allready added to lhsSpan
        
        const std::vector<Index> lhsIndices = _lhs.get_assigned_indices(lhsSpan);
        const std::vector<Index> rhsIndices = _rhs.get_assigned_indices(rhsSpan);
        
        std::vector<Index> reorderedBaseIndices;
        reorderedBaseIndices.reserve(baseIndices.size());
        
        std::vector<size_t> reorderedBaseDimensions, lhsDims, rhsDims;
        reorderedBaseDimensions.reserve(_base.degree());
        lhsDims.reserve(lhsIndices.size());
        rhsDims.reserve(rhsIndices.size());
        
        _lhsSize=1;
        _rhsSize=1;
        
        
        // Work through the indices of lhs
        for(size_t i = 0; i < lhsIndices.size()-1; ++i) {
            REQUIRE(!contains(rhsIndices, lhsIndices[i]), "Left and right part of factorization must not share common indices.");
            
            // Find index in A and get dimension offset
            size_t j, dimOffset = 0;
            for(j = 0; lhsIndices[i] != baseIndices[j]; ++j) {
                REQUIRE(j < baseIndices.size(), " All but the last Index of the left part must be contained in the input Tensor of factorization.");
                dimOffset += baseIndices[j].span;
            }
            
            reorderedBaseIndices.push_back(baseIndices[j]);
            for(size_t k = 0; k < baseIndices[j].span; ++k) {
                reorderedBaseDimensions.push_back(_base.tensorObjectReadOnly->dimensions.at(dimOffset+k));
                lhsDims.push_back(_base.tensorObjectReadOnly->dimensions[dimOffset+k]);
                _lhsSize *= _base.tensorObjectReadOnly->dimensions[dimOffset+k];
            }
        }

        // Work through the indices of rhs
        for(size_t i = 1; i < rhsIndices.size(); ++i) {
            // Find index in A and get dimension offset
            size_t j, dimOffset = 0;
            for(j = 0; rhsIndices[i] != baseIndices[j]; ++j) {
                REQUIRE(j < baseIndices.size(), " All but the first Index of the right part must be contained in the input Tensor of factorization");
                dimOffset += baseIndices[j].span;
            }
            
            reorderedBaseIndices.push_back(baseIndices[j]);
            for(size_t k = 0; k < baseIndices[j].span; ++k) {
                reorderedBaseDimensions.push_back(_base.tensorObjectReadOnly->dimensions.at(dimOffset+k));
                rhsDims.push_back(_base.tensorObjectReadOnly->dimensions[dimOffset+k]);
                _rhsSize *= _base.tensorObjectReadOnly->dimensions[dimOffset+k];
            }
        }
        
        IndexedTensor<Tensor> reorderedBaseTensor(_base.tensorObjectReadOnly->construct_new(std::move(reorderedBaseDimensions), DONT_SET_ZERO()), std::move(reorderedBaseIndices), false);
        evaluate(reorderedBaseTensor, _base);
        reorderedBaseTensor.tensorObject->ensure_own_data();
        
        _rank = std::min(_lhsSize, _rhsSize);
        lhsDims.push_back(_rank);
        rhsDims.insert(rhsDims.begin(), _rank);
        
        _lhs.tensorObject->reset(std::move(lhsDims), DONT_SET_ZERO());
        _lhs.tensorObject->ensure_own_data_no_copy();
        _lhs.check_indices(false);
        
        _rhs.tensorObject->reset(std::move(rhsDims), DONT_SET_ZERO());
        _rhs.tensorObject->ensure_own_data_no_copy();
        _rhs.check_indices(false);
        
        return std::unique_ptr<Tensor>(reorderedBaseTensor.tensorObject);
    }
    
    void SVD::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
        REQUIRE(_output.size() == 3, "SVD requires two output tensors, not " << _output.size());
        const IndexedTensorReadOnly<Tensor>& A = input;
        const IndexedTensorWritable<Tensor>& U = *_output[0];
        const IndexedTensorWritable<Tensor>& S = *_output[1];
        const IndexedTensorWritable<Tensor>& Vt = *_output[2];
        
        REQUIRE(S.indices.size() == 2, "S must be a diagonal matrix and therefore must have two indices.");
        REQUIRE(U.indices.back() == S.indices[0], "The last index of U must conincide with the first of S");
        REQUIRE(Vt.indices[0] == S.indices[1], " The second index of S must coincide with the first of Vt");
        REQUIRE(!U.tensorObject->is_sparse() && !Vt.tensorObject->is_sparse(), "U and Vt have to be FullTensors, as they are defenitely not sparse.");
        
        size_t lhsSize, rhsSize, rank;
        
        std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, A, U, Vt);
        
        std::unique_ptr<value_t[]> tmpS(new value_t[rank]);
        
        if(reorderedBaseTensor->is_sparse()){
            LOG(fatal, "Sparse SVD not yet implemented.");
        } else {
            blasWrapper::svd(static_cast<FullTensor*>(U.tensorObject)->data.get(), tmpS.get(), static_cast<FullTensor*>(Vt.tensorObject)->data.get(), static_cast<FullTensor*>(reorderedBaseTensor.get())->data.get(), lhsSize, rhsSize);
        }
        
        // Determine the real rank
        for(size_t j = 1; j<rank; ++j) {
            REQUIRE(tmpS[j] >= 0, "Singular Value must be >= 0");
            if (tmpS[j] <= epsilon) {
                rank = j;
                break;
            }
        }
        
        //Apply factor to the diagonal matrix 
        array_scale(tmpS.get(), reorderedBaseTensor->factor, rank);
        
        // Create tensor from diagonal values
        S.tensorObject->reset(std::vector<size_t>(2, rank));
        if(S.tensorObject->is_sparse()) {
            for(size_t i = 0; i < rank; ++i) {
                static_cast<FullTensor&>(*S.tensorObject)[{i,i}] = tmpS[i];
            }
        } else {
            value_t* const dataPtr =  static_cast<FullTensor*>(S.tensorObject)->data.get();
            for(size_t i = 0; i < rank; ++i) {
                dataPtr[i*rank+i] = tmpS[i];
            }
        }
        
        static_cast<FullTensor*>(U.tensorObject)->resize_dimension(U.degree()-1, rank);
        static_cast<FullTensor*>(Vt.tensorObject)->resize_dimension(0, rank);
    }


    void QR::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
        REQUIRE(_output.size() == 2, "QR factorisation requires two output tensors, not " << _output.size());
        const IndexedTensorReadOnly<Tensor>& A = *input;
        const IndexedTensorWritable<Tensor>& Q = *_output[0];
        const IndexedTensorWritable<Tensor>& R = *_output[1];
         
        REQUIRE(Q.indices.back() == R.indices[0], "The last index of Q must coincide with the first index of R in QR factorization.");
        REQUIRE(!Q.tensorObject->is_sparse() && !R.tensorObject->is_sparse(), "Q and R have to be FullTensors, as they are defenitely not sparse.");
        
        size_t lhsSize, rhsSize, rank;
        
        std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, A, Q, R);
        
        REQUIRE(Q.degree()+R.degree() == A.degree()+2, "Q and R must have have the degree of the input tensor plus two in RQ factorization.");
        
        // R has to carry the constant factor
        R.tensorObject->factor = reorderedBaseTensor->factor;
        
        if(reorderedBaseTensor->is_sparse()){
            LOG(fatal, "Sparse QR not yet implemented.");
        } else {
            blasWrapper::qr_destructive(static_cast<FullTensor*>(Q.tensorObject)->data.get(), static_cast<FullTensor*>(R.tensorObject)->data.get(), static_cast<const FullTensor*>(reorderedBaseTensor.get())->data.get(), lhsSize, rhsSize);
        }
    }

    void RQ::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
        REQUIRE(_output.size() == 2, "RQ factorisation requires two output tensors, not " << _output.size());
        const IndexedTensorReadOnly<Tensor>& A = *input;
        const IndexedTensorWritable<Tensor>& R = *_output[0];
        const IndexedTensorWritable<Tensor>& Q = *_output[1];
        
        REQUIRE(Q.indices[0] == R.indices.back(), "The last index of R must coincide with the first index of Q in RQ factorization.");
        REQUIRE(!Q.tensorObject->is_sparse() && !R.tensorObject->is_sparse(), "Q and R have to be FullTensors, as they are defenitely not sparse.");
        
        size_t lhsSize, rhsSize, rank;
        
        std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, A, R, Q);
        
        REQUIRE(Q.degree()+R.degree() == A.degree()+2, "Q and R must have have the degree of A plus two in RQ factorization.");
        
        // R has to carry the constant factor
        R.tensorObject->factor = reorderedBaseTensor->factor;
        
        blasWrapper::rq_destructive(static_cast<FullTensor*>(R.tensorObject)->data.get(), static_cast<FullTensor*>(Q.tensorObject)->data.get(), static_cast<const FullTensor*>(reorderedBaseTensor.get())->data.get(), lhsSize, rhsSize);
    }
}
