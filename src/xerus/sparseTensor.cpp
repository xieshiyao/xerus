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

#include <xerus/sparseTensor.h>

namespace xerus {
    SparseTensor::SparseTensor() : Tensor(), entries(new std::map<size_t, value_t>()) {}

    SparseTensor::SparseTensor( const SparseTensor&  _other) : Tensor(_other), entries(_other.entries) { }
    
    SparseTensor::SparseTensor(       SparseTensor&& _other) : Tensor(std::move(_other)), entries(std::move(_other.entries)) { }
    
    SparseTensor::SparseTensor(const std::vector<size_t> & _dimensions) : Tensor(_dimensions), entries(new std::map<size_t, value_t>()) { }
        
    SparseTensor::SparseTensor(      std::vector<size_t>&& _dimensions) : Tensor(std::move(_dimensions)), entries(new std::map<size_t, value_t>()) { }
    
    SparseTensor::SparseTensor(std::initializer_list<size_t>&& _dimensions) : Tensor(std::move(_dimensions)), entries(new std::map<size_t, value_t>()) { }
    
    SparseTensor::SparseTensor(const FullTensor& _full, const double _eps) : Tensor(_full), entries(new std::map<size_t, value_t>()) {
        const value_t* const fullDataPtr = _full.data.get();
        
        for(size_t i = 0; i < _full.size; ++i) {
            if(std::abs(fullDataPtr[i]) >= _eps) {
                entries->insert({i, fullDataPtr[i]});
//                 entries->emplace_hint(entries->end(), i, fullDataPtr[i]); TODO use this instead
            }
        }
    }
    
    Tensor* SparseTensor::get_copy() const {
        return new SparseTensor(*this);
    }
        
    Tensor* SparseTensor::get_moved_copy() {
        return new SparseTensor(std::move(*this));
    }
    
    Tensor* SparseTensor::construct_new() const {
        return new SparseTensor();
    }
    
    Tensor* SparseTensor::construct_new(const std::vector<size_t>&  _dimensions) const {
        return new SparseTensor(_dimensions);
    }
    
    Tensor* SparseTensor::construct_new(      std::vector<size_t>&& _dimensions) const {
        return new SparseTensor(std::move(_dimensions));
    }
    
    Tensor* SparseTensor::construct_new(const std::vector<size_t>&  _dimensions, DONT_SET_ZERO) const {
        return new SparseTensor(_dimensions);
    }
    
    Tensor* SparseTensor::construct_new(      std::vector<size_t>&& _dimensions, DONT_SET_ZERO) const {
        return new SparseTensor(std::move(_dimensions));
    }
        
        
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
    void SparseTensor::ensure_own_data() {
        if(!entries.unique()) {
            entries.reset(new std::map<size_t, value_t>(*entries));
        }
    }


    void SparseTensor::ensure_own_data_no_copy() {
        if(!entries.unique()) {
            entries.reset(new std::map<size_t, value_t>());
        }
    }

    void SparseTensor::apply_factor() {
        if(has_factor()) {
            ensure_own_data();
            for(std::pair<const size_t, value_t>& entry : *entries) {
                entry.second *= factor;
            }
            factor = 1.0;
        }
    }

    void SparseTensor::ensure_own_data_and_apply_factor() {
        ensure_own_data();
        if(has_factor()) {
            for(std::pair<const size_t, value_t>& entry : *entries) {
                entry.second *= factor;
            }
            factor = 1.0;
        }
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    SparseTensor& SparseTensor::operator=(const SparseTensor&  _other) {
        assign(_other);
        entries = _other.entries;
        return *this;
    }
    
    SparseTensor& SparseTensor::operator=(      SparseTensor&& _other) {
        assign(std::move(_other));
        entries = _other.entries;
        return *this;
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    value_t& SparseTensor::operator[](const size_t _position) {
        apply_factor();
        return (*entries)[_position];
    }
    
    value_t SparseTensor::operator[](const size_t _position) const {
        const std::map<size_t, value_t>::const_iterator entry = entries->find(_position);
        if(entry == entries->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    value_t& SparseTensor::operator[](const std::vector<size_t>& _indices) {
        apply_factor();
        CHECK(_indices.size() == dimensions.size(), fatal, "Wrong number of indices given");
        size_t finalIndex = 0;
        for(size_t i = 0; i < _indices.size(); ++i) {
            REQUIRE(_indices[i] < dimensions[i], "Index "<< i <<" out of bounds " << _indices[i] << " >=! " << dimensions[i]);
            finalIndex *= dimensions[i];
            finalIndex += _indices[i];
        }
        return (*entries)[finalIndex];
    }
    
    value_t SparseTensor::operator[](const std::vector<size_t>& _indices) const {
        CHECK(_indices.size() == dimensions.size(), fatal, "Wrong number of indices given");
        size_t finalIndex = 0;
        for(size_t i = 0; i < _indices.size(); ++i) {
            REQUIRE(_indices[i] < dimensions[i], "Index "<< i <<" out of bounds " << _indices[i] << " >=! " << dimensions[i]);
            finalIndex *= dimensions[i];
            finalIndex += _indices[i];
        }
        
        const std::map<size_t, value_t>::const_iterator entry = entries->find(finalIndex);
        if(entry == entries->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    value_t SparseTensor::at(const size_t _position) const {
        const std::map<size_t, value_t>::const_iterator entry = entries->find(_position);
        if(entry == entries->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    value_t SparseTensor::at(const std::vector<size_t>& _indices) const {
        CHECK(_indices.size() == dimensions.size(), fatal, "Wrong number of indices given");
        size_t finalIndex = 0;
        for(size_t i = 0; i < _indices.size(); ++i) {
            REQUIRE(_indices[i] < dimensions[i], "Index "<< i <<" out of bounds " << _indices[i] << " >=! " << dimensions[i]);
            finalIndex *= dimensions[i];
            finalIndex += _indices[i];
        }
        
        const std::map<size_t, value_t>::const_iterator entry = entries->find(finalIndex);
        if(entry == entries->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    SparseTensor& SparseTensor::operator+=(const SparseTensor& _other) {
		PA_START;
        // We want a*x+b*y and turn it into b*((a/b)*x+y)
        factor /= _other.factor;
        ensure_own_data_and_apply_factor();
        factor = _other.factor;
        
        for(const std::pair<size_t, value_t>& entry : *_other.entries) {
            std::pair<std::map<size_t, value_t>::iterator, bool> result = entries->insert(entry);
            if(!result.second) {
                result.first->second += entry.second;
            }
        }
        PA_END("ADD/SUB", "SparseTensor ADD/SUB SparseTensor", misc::to_string(size));
        return *this;
    }
    
    SparseTensor  SparseTensor::operator+(const SparseTensor& _other) const {
        SparseTensor ret(*this);
        ret += _other;
        return ret;
    }
    
    SparseTensor& SparseTensor::operator-=(const SparseTensor& _other) {
		PA_START;
        // We want a*x-b*y and turn it into -b*((-a/b)*x+y)
        factor /= -_other.factor;
        ensure_own_data_and_apply_factor();
        factor = -_other.factor;
        
        for(const std::pair<size_t, value_t>& entry : *_other.entries) {
            std::pair<std::map<size_t, value_t>::iterator, bool> result = entries->insert(entry);
            if(!result.second) {
                result.first->second += entry.second;
            }
        }
        PA_END("ADD/SUB", "SparseTensor ADD/SUB SparseTensor", misc::to_string(size));
        return *this;
    }
    
    SparseTensor SparseTensor::operator-(const SparseTensor& _other) const {
        SparseTensor ret(*this);
        ret -= _other;
        return ret;
    }

    
    SparseTensor SparseTensor::operator*(const value_t _prod) const {
        SparseTensor ret(*this);
        ret.factor *= _prod;
        return ret;
    }

    
    SparseTensor SparseTensor::operator/(const value_t _div) const {
        SparseTensor ret(*this);
        ret.factor /= _div;
        return ret;
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Higher functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    void SparseTensor::reset(const std::vector<size_t>&  _newDim, DONT_SET_ZERO) {
        change_dimensions(_newDim);
        factor = 1.0;
        if(entries.unique()) {
            entries->clear();
        } else {
            entries.reset(new std::map<size_t, value_t>());
        }
    }
    
    void SparseTensor::reset(      std::vector<size_t>&& _newDim, DONT_SET_ZERO) {
        change_dimensions(std::move(_newDim));
        factor = 1.0;
        if(entries.unique()) {
            entries->clear();
        } else {
            entries.reset(new std::map<size_t, value_t>());
        }
    }
    
    void SparseTensor::reset(const std::vector<size_t>&  _newDim) {
        change_dimensions(_newDim);
        factor = 1.0;
        if(entries.unique()) {
            entries->clear();
        } else {
            entries.reset(new std::map<size_t, value_t>());
        }
    }
    
    void SparseTensor::reset(      std::vector<size_t>&& _newDim) {
        change_dimensions(std::move(_newDim));
        factor = 1.0;
        if(entries.unique()) {
            entries->clear();
        } else {
            entries.reset(new std::map<size_t, value_t>());
        }
    }
    
    bool SparseTensor::is_sparse() const {
        return true;
    }
    

    size_t SparseTensor::count_non_zero_entries(const value_t _eps) const {
        size_t count = 0;
        for(const std::pair<size_t, value_t>& entry : *entries) {
            if(std::abs(entry.second) > _eps) { count++; } 
        }
        return count;
    }
    
    value_t SparseTensor::frob_norm() const {
        value_t norm = 0;
        for(const std::pair<size_t, value_t>& entry : *entries) {
            norm += misc::sqr(entry.second);
        }
        return std::abs(factor)*sqrt(norm);
    }
    
    std::string SparseTensor::to_string() const {
        if (degree() == 0) return xerus::misc::to_string(at(0));
        std::string result;
        for (size_t i=0; i<size; ++i) {
            result += xerus::misc::to_string(at(i)) + " ";
            if ((i+1) % (size / dimensions[0]) == 0) {
                result += '\n';
            } else if (degree() > 1 && (i+1) % (size / dimensions[0] / dimensions[1]) == 0) {
                result += '\t';
            } else if (degree() > 2 && (i+1) % (size / dimensions[0] / dimensions[1] / dimensions[2]) == 0) {
                result += "/ ";
            }
        }
        return result;
    }
    
    bool SparseTensor::compare_to_data(std::vector<value_t> _values, const double _eps) const {
        REQUIRE(entries, "Internal Error");
        if(size != _values.size()) { return false; }
        for(size_t i=0; i < size; ++i) {
            if(std::abs(at(i)-_values[i]) > _eps) { return false; }
        }
        return true;
    }

    bool SparseTensor::compare_to_data(const value_t* _values, const double _eps) const {
        for(size_t i=0; i < size; ++i) {
            if(std::abs(at(i)-_values[i]) > _eps) { return false; }
        }
        return true;
    }
}
