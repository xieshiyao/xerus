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

#include "indexedTensor_tensor_factorisations.h"

namespace xerus {
    class TensorNode {
    public:
        struct Link {
            /// The index of the otherNode this Link links to. (Deprecated: < 0 if this is an external index (==-1) )
            size_t other; 
            /// IndexPosition on the other node or index of external index
            size_t indexPosition;
            /// Always equals to other->tensorObject->dimensions[indexPosition]
            size_t dimension;
            /// Flag to mark Links that correspond to external indices TODO: use this
            bool external;
            
            Link() {}
            Link(const size_t _other, const size_t _indexPos, const size_t _dim, const bool _external) : other(_other), indexPosition(_indexPos), dimension(_dim), external(_external) {}
            
            bool links(const size_t _other) const { return !external && other == _other; }
        };

        std::shared_ptr<Tensor> tensorObject;
        std::vector<Link> neighbors;
        bool erased;
        
        explicit TensorNode() : erased(true) { }
        
        ALLOW_MOVE(std::shared_ptr<Tensor>, shared_tensor_ptr)
        explicit TensorNode(shared_tensor_ptr&& _tensorObject) : tensorObject(std::forward<shared_tensor_ptr>(_tensorObject)), erased(false) {}
        
        
        template<ADD_MOVE(std::shared_ptr<Tensor>, SPT), ADD_MOVE(std::vector<Link>, VL)>
        explicit TensorNode(SPT&& _tensorObject, VL&& _neighbors) : tensorObject(std::forward<SPT>(_tensorObject)), neighbors(std::forward<VL>(_neighbors)), erased(false) {}
        
        TensorNode strippped_copy() const {
            return TensorNode(std::shared_ptr<Tensor>(), neighbors);
        }
        
        void ensure_own_tensor() {
            if(!tensorObject.unique()) {
                tensorObject.reset(tensorObject->get_copy());
            }
        }
        
        void add_factor(const value_t _factor) {
            ensure_own_tensor();
            tensorObject->factor *= _factor;
        }
        
        // All getters are written without the use of tensorObject so that they also work for empty nodes
        
        size_t size() const {
            size_t s = 1;
            for (const Link &l : neighbors) {
                s *= l.dimension;
            }
            return s;
        }
        
        size_t degree() const {
            return neighbors.size();
        }
        
        void erase() {
            erased = true;
            neighbors.clear();
            tensorObject.reset();
        }
    };
}

std::ostream &operator<<(std::ostream &_out, const xerus::TensorNode::Link &_rhs);

