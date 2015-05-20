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
    
    TensorNode::TensorNode(const TensorNode&  _other) : tensorObject(_other.tensorObject ? _other.tensorObject->get_copy() : nullptr), neighbors(_other.neighbors), erased(_other.erased) { }
    
    TensorNode::TensorNode(      TensorNode&& _other) : tensorObject(std::move(_other.tensorObject)), neighbors(std::move(_other.neighbors)), erased(_other.erased) { }
    
    TensorNode::TensorNode(      std::unique_ptr<Tensor>&& _tensorObject) : tensorObject(std::move(_tensorObject)), neighbors(), erased(false) {}
    
    TensorNode::TensorNode(std::unique_ptr<Tensor>&& _tensorObject, const std::vector<Link>& _neighbors) : tensorObject(std::move(_tensorObject)), neighbors(_neighbors), erased(false) {}
    
    TensorNode::TensorNode(std::unique_ptr<Tensor>&& _tensorObject,       std::vector<Link>&& _neighbors) : tensorObject(std::move(_tensorObject)), neighbors(std::move(_neighbors)), erased(false) {}
    
    TensorNode& TensorNode::operator=(const TensorNode&  _other) {
        if(_other.tensorObject) {
            tensorObject.reset(_other.tensorObject->get_copy());
        }
        neighbors = _other.neighbors;
        erased = _other.erased;
        return *this;
    }
    
    TensorNode& TensorNode::operator=(      TensorNode&& _other) {
        tensorObject = std::move(_other.tensorObject);
        neighbors = std::move(_other.neighbors);
        erased = _other.erased;
        return *this;
    }
    
    TensorNode TensorNode::strippped_copy() const {
        return TensorNode(std::unique_ptr<Tensor>(), neighbors);
    }
    
    void TensorNode::add_factor(const value_t _factor) {
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
        _out << "L{";
		if (_rhs.external) {
			_out << "ext";
		} else {
			_out << _rhs.other;
		}
		_out << " (" << _rhs.indexPosition << "), dim " << _rhs.dimension << "}";
        return _out;
    }
}
