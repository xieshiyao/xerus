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
 * @brief Implementation of the FullTensor factorizations.
 */

#include <xerus/indexedTensor_tensor_factorisations.h>
#include <xerus/index.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>
#include <xerus/indexedTensor_tensor_operators.h>
#include <xerus/blasLapackWrapper.h>

namespace xerus {
	
	std::unique_ptr<Tensor> prepare_split(size_t& _lhsSize, size_t& _rhsSize, size_t& _rank, std::vector<Index>& _lhsPreliminaryIndices, std::vector<Index>& _rhsPreliminaryIndices, const IndexedTensorReadOnly<Tensor>& _base, const IndexedTensorWritable<Tensor>& _lhs, const IndexedTensorWritable<Tensor>& _rhs) {
		const std::vector<Index> baseIndices = _base.get_assigned_indices();
		
		// Calculate the future order of lhs and rhs.
		size_t lhsOrder = 1, rhsOrder = 1; // Start with 1 because there is a new dimension introduced in the split. 
		
		for (const Index& idx : baseIndices) {
			if (idx.open()) {
				if(misc::contains(_lhs.indices, idx)) {
					lhsOrder += idx.span; 
				} else {
					REQUIRE(misc::contains(_rhs.indices, idx), "Every open index of factorisation base must be contained in one of the targets");
					rhsOrder += idx.span;
				}
			}
		}
		
		const std::vector<Index> lhsIndices = _lhs.get_assigned_indices(lhsOrder);
		const std::vector<Index> rhsIndices = _rhs.get_assigned_indices(rhsOrder);
		
		std::vector<Index> reorderedBaseIndices;
		reorderedBaseIndices.reserve(baseIndices.size());
		
// 		std::vector<Index> lhsPreliminaryIndices;
		_lhsPreliminaryIndices.reserve(lhsIndices.size());
		
// 		std::vector<Index> rhsPreliminaryIndices;
		_rhsPreliminaryIndices.reserve(rhsIndices.size());
		
		std::vector<size_t> reorderedBaseDimensions, lhsDims, rhsDims;
		reorderedBaseDimensions.reserve(_base.degree());
		lhsDims.reserve(lhsIndices.size());
		rhsDims.reserve(rhsIndices.size());
		
		_lhsSize=1;
		_rhsSize=1;
		
		Index auxiliaryIndex;

		// Work through the indices of lhs
		IF_CHECK(bool foundCommon = false;)
		for(size_t i = 0; i < lhsIndices.size(); ++i) {
			// Find index in A and get dimension offset
			size_t j, dimOffset = 0;
			for(j = 0; j < baseIndices.size() && lhsIndices[i] != baseIndices[j]; ++j) {
				dimOffset += baseIndices[j].span;
			}
			
			if(j < baseIndices.size()) {
				_lhsPreliminaryIndices.push_back(baseIndices[j]);
				reorderedBaseIndices.push_back(baseIndices[j]);
				for(size_t k = 0; k < baseIndices[j].span; ++k) {
					reorderedBaseDimensions.push_back(_base.tensorObjectReadOnly->dimensions.at(dimOffset+k));
					lhsDims.push_back(_base.tensorObjectReadOnly->dimensions[dimOffset+k]);
					_lhsSize *= _base.tensorObjectReadOnly->dimensions[dimOffset+k];
				}
			} else {
				REQUIRE(!foundCommon, "Left part of factorization must have exactly one index that is not contained in base. Here it is more than one.");
				IF_CHECK(foundCommon = true;)
				auxiliaryIndex = lhsIndices[i];
			}
		}
		_lhsPreliminaryIndices.push_back(auxiliaryIndex);

		// Work through the indices of rhs
		IF_CHECK(foundCommon = false;)
		for(size_t i = 0; i < rhsIndices.size(); ++i) {
			// Find index in A and get dimension offset
			size_t j, dimOffset = 0;
			for(j = 0; j < baseIndices.size() && rhsIndices[i] != baseIndices[j]; ++j) {
				dimOffset += baseIndices[j].span;
			}
			
			if(j < baseIndices.size()) {
				_rhsPreliminaryIndices.push_back(baseIndices[j]);
				reorderedBaseIndices.push_back(baseIndices[j]);
				for(size_t k = 0; k < baseIndices[j].span; ++k) {
					reorderedBaseDimensions.push_back(_base.tensorObjectReadOnly->dimensions.at(dimOffset+k));
					rhsDims.push_back(_base.tensorObjectReadOnly->dimensions[dimOffset+k]);
					_rhsSize *= _base.tensorObjectReadOnly->dimensions[dimOffset+k];
				}
			} else {
				REQUIRE(!foundCommon, "Right part of factorization must have exactly one index that is not contained in base. Here it is more than one.");
				IF_CHECK(foundCommon = true;)
				auxiliaryIndex = rhsIndices[i];
			}
		}
		_rhsPreliminaryIndices.insert(_rhsPreliminaryIndices.begin(), auxiliaryIndex);
		
		IndexedTensor<Tensor> reorderedBaseTensor(_base.tensorObjectReadOnly->construct_new(std::move(reorderedBaseDimensions), DONT_SET_ZERO()), std::move(reorderedBaseIndices), false);
		evaluate(reorderedBaseTensor, _base);
		reorderedBaseTensor.tensorObject->ensure_own_data();
		
		_rank = std::min(_lhsSize, _rhsSize);
		lhsDims.push_back(_rank);
		rhsDims.insert(rhsDims.begin(), _rank);
		
		_lhs.tensorObject->reset(std::move(lhsDims), DONT_SET_ZERO());
		_lhs.tensorObject->ensure_own_data_no_copy();
		IF_CHECK( _lhs.check_indices(false); )
		
		_rhs.tensorObject->reset(std::move(rhsDims), DONT_SET_ZERO());
		_rhs.tensorObject->ensure_own_data_no_copy();
		IF_CHECK( _rhs.check_indices(false); )
		
		return std::unique_ptr<Tensor>(reorderedBaseTensor.tensorObject);
	}
	
	void SVD::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
		REQUIRE(_output.size() == 3, "SVD requires two output tensors, not " << _output.size());
		const IndexedTensorReadOnly<Tensor>& A = input;
		const IndexedTensorWritable<Tensor>& U = *_output[0];
		const IndexedTensorWritable<Tensor>& S = *_output[1];
		const IndexedTensorWritable<Tensor>& Vt = *_output[2];
		
		IF_CHECK(S.check_indices(2, false));
		REQUIRE(!U.tensorObject->is_sparse() && !Vt.tensorObject->is_sparse(), "U and Vt have to be FullTensors, as they are defenitely not sparse.");
		REQUIRE(epsilon < 1, "Epsilon must be smaller than one.");
		REQUIRE(maxRank > 0, "maxRank must be larger than zero.");
		
		size_t lhsSize, rhsSize, rank;
		std::vector<Index> lhsPreliminaryIndices, rhsPreliminaryIndices;
		
		std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, lhsPreliminaryIndices, rhsPreliminaryIndices, A, U, Vt);
		
		std::unique_ptr<value_t[]> tmpS(new value_t[rank]);
		
		// Calculate the actual SVD
		if(reorderedBaseTensor->is_sparse()) {
			LOG(fatal, "Sparse SVD not yet implemented.");
		} else {
			blasWrapper::svd(static_cast<FullTensor*>(U.tensorObject)->data.get(), tmpS.get(), static_cast<FullTensor*>(Vt.tensorObject)->data.get(), static_cast<FullTensor*>(reorderedBaseTensor.get())->data.get(), lhsSize, rhsSize);
		}
		
		// Apply factor to the diagonal matrix
		misc::array_scale(tmpS.get(), reorderedBaseTensor->factor, rank);
		
		// Account for hard threshold
		rank = std::min(rank, maxRank);
		
		// Apply soft threshold and determine the real rank
		double realSoftThreshold;
		if(preventZero) {
			realSoftThreshold = std::min((1-EPSILON)*tmpS[0], softThreshold);
		} else {
			realSoftThreshold = softThreshold;
		}
		
		tmpS[0] = std::max(0.0, tmpS[0] - realSoftThreshold);
		for(size_t j = 1; j < rank; ++j) {
			tmpS[j] -= realSoftThreshold;
			if (tmpS[j] <= epsilon*tmpS[0]) {
				rank = j;
				break;
			}
		}
		
		// Create tensor from diagonal values
		S.tensorObject->reset(std::vector<size_t>(2, rank));
		if(S.tensorObject->is_sparse()) {
			for(size_t i = 0; i < rank; ++i) {
				static_cast<SparseTensor&>(*S.tensorObject)[i*rank+i] = tmpS[i];
			}
		} else {
			value_t* const dataPtr =  static_cast<FullTensor*>(S.tensorObject)->data.get();
			for(size_t i = 0; i < rank; ++i) {
				dataPtr[i*rank+i] = tmpS[i];
			}
		}
		
		static_cast<FullTensor*>(U.tensorObject)->resize_dimension(U.degree()-1, rank);
		static_cast<FullTensor*>(Vt.tensorObject)->resize_dimension(0, rank);
		
		// Post evaluate the results
		std::vector<Index> midPreliminaryIndices({lhsPreliminaryIndices.back(), rhsPreliminaryIndices.front()});
		U = (*U.tensorObjectReadOnly)(lhsPreliminaryIndices);
		S = (*S.tensorObjectReadOnly)(midPreliminaryIndices);
		Vt = (*Vt.tensorObjectReadOnly)(rhsPreliminaryIndices);
		
		REQUIRE(U.tensorObject->all_entries_valid(), "Internal Error");
		REQUIRE(S.tensorObject->all_entries_valid(), "Internal Error");
		REQUIRE(Vt.tensorObject->all_entries_valid(), "Internal Error");
	}


	void QR::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
		REQUIRE(_output.size() == 2, "QR factorisation requires two output tensors, not " << _output.size());
		const IndexedTensorReadOnly<Tensor>& A = *input;
		const IndexedTensorWritable<Tensor>& Q = *_output[0];
		const IndexedTensorWritable<Tensor>& R = *_output[1];
		 
		REQUIRE(!Q.tensorObject->is_sparse() && !R.tensorObject->is_sparse(), "Q and R have to be FullTensors, as they are defenitely not sparse.");
		
		size_t lhsSize, rhsSize, rank;
		std::vector<Index> lhsPreliminaryIndices, rhsPreliminaryIndices;
		
		std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, lhsPreliminaryIndices, rhsPreliminaryIndices, A, Q, R);
		
		// R has to carry the constant factor
		R.tensorObject->factor = reorderedBaseTensor->factor;
		
		if(reorderedBaseTensor->is_sparse()) {
			LOG(fatal, "Sparse QR not yet implemented.");
		} else {
			blasWrapper::qr_destructive(static_cast<FullTensor*>(Q.tensorObject)->data.get(), static_cast<FullTensor*>(R.tensorObject)->data.get(), static_cast<const FullTensor*>(reorderedBaseTensor.get())->data.get(), lhsSize, rhsSize);
		}
		
		// Post evaluate the results
		Q = (*Q.tensorObjectReadOnly)(lhsPreliminaryIndices);
		R = (*R.tensorObjectReadOnly)(rhsPreliminaryIndices);
	}

	void RQ::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
		REQUIRE(_output.size() == 2, "RQ factorisation requires two output tensors, not " << _output.size());
		const IndexedTensorReadOnly<Tensor>& A = *input;
		const IndexedTensorWritable<Tensor>& R = *_output[0];
		const IndexedTensorWritable<Tensor>& Q = *_output[1];
		
		REQUIRE(!Q.tensorObject->is_sparse() && !R.tensorObject->is_sparse(), "Q and R have to be FullTensors, as they are defenitely not sparse.");
		
		size_t lhsSize, rhsSize, rank;
		std::vector<Index> lhsPreliminaryIndices, rhsPreliminaryIndices;
		
		std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, lhsPreliminaryIndices, rhsPreliminaryIndices, A, R, Q);
		
		// R has to carry the constant factor
		R.tensorObject->factor = reorderedBaseTensor->factor;
		
		if(reorderedBaseTensor->is_sparse()) {
			LOG(fatal, "Sparse QR not yet implemented.");
		} else {
			blasWrapper::rq_destructive(static_cast<FullTensor*>(R.tensorObject)->data.get(), static_cast<FullTensor*>(Q.tensorObject)->data.get(), static_cast<const FullTensor*>(reorderedBaseTensor.get())->data.get(), lhsSize, rhsSize);
		}

		// Post evaluate the results
		R = (*R.tensorObjectReadOnly)(lhsPreliminaryIndices);
		Q = (*Q.tensorObjectReadOnly)(rhsPreliminaryIndices);
	}
	
	
	void QC::operator()(const std::vector<const IndexedTensorWritable<Tensor>*>& _output) const {
		REQUIRE(_output.size() == 2, "QC factorisation requires two output tensors, not " << _output.size());
		const IndexedTensorReadOnly<Tensor>& A = *input;
		const IndexedTensorWritable<Tensor>& Q = *_output[0];
		const IndexedTensorWritable<Tensor>& C = *_output[1];
		
		REQUIRE(!Q.tensorObject->is_sparse() && !C.tensorObject->is_sparse(), "Q and C have to be FullTensors, as they are definitely not sparse.");
		
		size_t lhsSize, rhsSize, rank;
		std::vector<Index> lhsPreliminaryIndices, rhsPreliminaryIndices;
		
		std::unique_ptr<Tensor> reorderedBaseTensor = prepare_split(lhsSize, rhsSize, rank, lhsPreliminaryIndices, rhsPreliminaryIndices, A, Q, C);
		
		// C has to carry the constant factor
		C.tensorObject->factor = reorderedBaseTensor->factor;
		
		if(reorderedBaseTensor->is_sparse()){
			LOG(fatal, "Sparse QC not yet implemented.");
		} else {
			std::unique_ptr<double[]> Qt, Ct;
			blasWrapper::qc(Qt, Ct, static_cast<const FullTensor*>(reorderedBaseTensor.get())->data.get(),  lhsSize, rhsSize, rank);
			
			// TODO there should be a reset function to use instead of directly accesing those values. -- There is, called reset
			static_cast<FullTensor*>(Q.tensorObject)->data.reset(Qt.release(), &internal::array_deleter_vt);
			Q.tensorObject->dimensions.back() = rank; 
			Q.tensorObject->size = misc::product(Q.tensorObject->dimensions);
			static_cast<FullTensor*>(C.tensorObject)->data.reset(Ct.release(), &internal::array_deleter_vt);
			C.tensorObject->dimensions.front() = rank;
			C.tensorObject->size = misc::product(C.tensorObject->dimensions);
		}
		
		// Post evaluate the results
		Q = (*Q.tensorObjectReadOnly)(lhsPreliminaryIndices);
		C = (*C.tensorObjectReadOnly)(rhsPreliminaryIndices);
	}
}
