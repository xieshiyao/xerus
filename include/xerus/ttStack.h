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
 * @brief Header file for the TTStack class.
 */

#pragma once

#include "misc/missingFunctions.h"
#include "misc/check.h"

#include "index.h"
#include "indexedTensor_tensor_factorisations.h"
#include "tensor.h"
#include "tensorNetwork.h"
#include "indexedTensorMoveable.h"
#include "indexedTensorList.h"

namespace xerus {
	namespace internal {
		template<bool isOperator>
		///@brief Internal class used to represent stacks of (possibly multiply) applications of TTOperators to either a TTTensor or TTOperator.
		class TTStack final : public TensorNetwork {
		public:
			const bool cannonicalization_required;
			
			const size_t futureCorePosition;
			
			explicit TTStack(const bool _canno, const size_t _corePos = 0) : cannonicalization_required(_canno), futureCorePosition(_corePos) {};
			
			virtual void operator*=(const value_t _factor) override {
				REQUIRE(nodes.size() > 0, "There must not be a TTNetwork without any node");
				
				if(cannonicalization_required) {
					*nodes[futureCorePosition+1].tensorObject *= _factor;
				} else if(degree() > 0) {
					*nodes[1].tensorObject *= _factor;
				} else {
					*nodes[0].tensorObject *= _factor;
				}
			}
			
			virtual void operator/=(const value_t _divisor) override {
				operator*=(1/_divisor);
			}
			
			/*- - - - - - - - - - - - - - - - - - - - - - - - - - Operator specializations - - - - - - - - - - - - - - - - - - - - - - - - - - */
			virtual void specialized_evaluation(IndexedTensorWritable<TensorNetwork>&& _me _unused_ , IndexedTensorReadOnly<TensorNetwork>&& _other _unused_) override {
				LOG(fatal, "TTStack not supported as a storing type");
			}
			virtual bool specialized_contraction(std::unique_ptr<IndexedTensorMoveable<TensorNetwork>>& _out, IndexedTensorReadOnly<TensorNetwork>&& _me, IndexedTensorReadOnly<TensorNetwork>&& _other) const override {
				return TTNetwork<isOperator>::specialized_contraction_f(_out, std::move(_me), std::move(_other));
			}
			virtual bool specialized_sum(std::unique_ptr<IndexedTensorMoveable<TensorNetwork>>& _out, IndexedTensorReadOnly<TensorNetwork>&& _me, IndexedTensorReadOnly<TensorNetwork>&& _other) const override {
				return TTNetwork<isOperator>::specialized_sum_f(_out, std::move(_me), std::move(_other));
			}
			
			virtual TensorNetwork* get_copy() const override {
				return new TTStack(*this);
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
