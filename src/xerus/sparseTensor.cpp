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
 * @brief Implementation of the SparseTensor class.
 */

#include <xerus/sparseTensor.h>

namespace xerus {
	SparseTensor::SparseTensor() : Tensor(Representation::Sparse) { }

    SparseTensor::SparseTensor( const SparseTensor&  _other) : Tensor(_other) { }
    
    SparseTensor::SparseTensor(       SparseTensor&& _other) : Tensor(std::move(_other)) { }
    
    SparseTensor::SparseTensor(const std::vector<size_t> & _dimensions) : Tensor(_dimensions, Representation::Sparse) { }
        
	SparseTensor::SparseTensor(      std::vector<size_t>&& _dimensions) : Tensor(std::move(_dimensions), Representation::Sparse) { }
    
    SparseTensor::SparseTensor(std::initializer_list<size_t>&& _dimensions) : Tensor(std::move(_dimensions), Representation::Sparse) { }
    
    SparseTensor::SparseTensor(const FullTensor& _full, const double _eps) : Tensor(_full) {
		denseData.reset();
		representation = Representation::Sparse;
		sparseData.reset(new std::map<size_t, value_t>());
        for(size_t i = 0; i < _full.size; ++i) {
            if(std::abs(_full[i]) >= _eps) {
                sparseData->insert({i, _full[i]});
//                 sparseData->emplace_hint(sparseData->end(), i, fullDataPtr[i]); TODO use this instead
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
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    SparseTensor& SparseTensor::operator=(const SparseTensor&  _other) {
        assign(_other);
		sparseData = _other.sparseData;
        return *this;
    }
    
    SparseTensor& SparseTensor::operator=(      SparseTensor&& _other) {
        assign(std::move(_other));
		sparseData = _other.sparseData;
        return *this;
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    value_t& SparseTensor::operator[](const size_t _position) {
		REQUIRE(_position < size, "position " << _position << " does not exist in (sparse)tensor of dimensions " << dimensions);
        apply_factor();
        return (*sparseData)[_position];
    }
    
    value_t SparseTensor::operator[](const size_t _position) const {
		REQUIRE(_position < size, "position " << _position << " does not exist in (sparse)tensor of dimensions " << dimensions);
		const std::map<size_t, value_t>::const_iterator entry = sparseData->find(_position);
		if(entry == sparseData->end()) {
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
        return (*sparseData)[finalIndex];
    }
    
    value_t SparseTensor::operator[](const std::vector<size_t>& _indices) const {
        CHECK(_indices.size() == dimensions.size(), fatal, "Wrong number of indices given");
        size_t finalIndex = 0;
        for(size_t i = 0; i < _indices.size(); ++i) {
            REQUIRE(_indices[i] < dimensions[i], "Index "<< i <<" out of bounds " << _indices[i] << " >=! " << dimensions[i]);
            finalIndex *= dimensions[i];
            finalIndex += _indices[i];
        }
        
        const std::map<size_t, value_t>::const_iterator entry = sparseData->find(finalIndex);
        if(entry == sparseData->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    value_t SparseTensor::at(const size_t _position) const {
        const std::map<size_t, value_t>::const_iterator entry = sparseData->find(_position);
        if(entry == sparseData->end()) {
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
        
        const std::map<size_t, value_t>::const_iterator entry = sparseData->find(finalIndex);
        if(entry == sparseData->end()) {
            return 0;
        } else {
            return factor*entry->second;
        }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    SparseTensor& SparseTensor::operator+=(const SparseTensor& _other) {
        ensure_own_data_and_apply_factor();
        
		PA_START;
        for(const std::pair<size_t, value_t>& entry : *_other.sparseData) {
            std::pair<std::map<size_t, value_t>::iterator, bool> result = sparseData->emplace(entry.first, _other.factor*entry.second);
            if(!result.second) {
                result.first->second += _other.factor*entry.second;
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
        ensure_own_data_and_apply_factor();
        
		PA_START;
        for(const std::pair<size_t, value_t>& entry : *_other.sparseData) {
            std::pair<std::map<size_t, value_t>::iterator, bool> result = sparseData->emplace(entry.first, -_other.factor*entry.second);
            if(!result.second) {
                result.first->second -= _other.factor*entry.second;
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
	
	void SparseTensor::fix_slate(const size_t _dimension, const size_t _slatePosition) {
		LOG(fatal, "Not yet implemented"); // TODO 
	}
    
    
    bool SparseTensor::is_sparse() const {
        return true;
    }
    

    size_t SparseTensor::count_non_zero_entries(const value_t _eps) const {
        size_t count = 0;
        for(const std::pair<size_t, value_t>& entry : *sparseData) {
            if(std::abs(entry.second) > _eps) { count++; } 
        }
        return count;
    }
    
    bool SparseTensor::all_entries_valid() const {
        for(const std::pair<size_t, value_t>& entry : *sparseData) {
            if(!std::isfinite(entry.second)) {return false; } 
        }
        return true;
    }
    
    value_t SparseTensor::frob_norm() const {
        value_t norm = 0;
        for(const std::pair<size_t, value_t>& entry : *sparseData) {
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
    
    bool SparseTensor::compare_to_data(const std::vector<value_t>& _values, const double _eps) const {
		REQUIRE(sparseData, "Internal Error");
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
