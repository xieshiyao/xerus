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

#include <xerus/tensorNode.h>
#include <xerus/tensor.h>


namespace xerus {
    
    TensorNode::TensorNode() : erased(true) { }
    
    TensorNode::TensorNode(const std::shared_ptr<Tensor>&  _tensorObject) : tensorObject(_tensorObject), erased(false) {}
    
    TensorNode::TensorNode(      std::shared_ptr<Tensor>&& _tensorObject) : tensorObject(std::move(_tensorObject)), erased(false) {}
    
    TensorNode TensorNode::strippped_copy() const {
        return TensorNode(std::shared_ptr<Tensor>(), neighbors);
    }
    
    void TensorNode::ensure_own_tensor() {
        if(!tensorObject.unique()) {
            tensorObject.reset(tensorObject->get_copy());
        }
    }
    
    void TensorNode::add_factor(const value_t _factor) {
        ensure_own_tensor();
        tensorObject->factor *= _factor;
    }
    
    
    size_t TensorNode::size() const {
        size_t s = 1;
        for (const Link &l : neighbors) {
            s *= l.dimension;
        }
        return s;
    }
    
    size_t TensorNode::degree() const {
        return neighbors.size();
    }
    
    void TensorNode::erase() {
        erased = true;
        neighbors.clear();
        tensorObject.reset();
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - External functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    std::ostream &operator<<(std::ostream &_out, const xerus::TensorNode::Link &_rhs) {
        _out << "L{" << _rhs.other << " (" << _rhs.indexPosition << "), " << _rhs.dimension << "}";
        return _out;
    }
}
