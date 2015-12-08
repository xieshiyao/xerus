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
 * @brief Implementation of the Tensor class.
 */

#include <xerus/tensor.h>
#include <xerus/misc/check.h>
#include <xerus/misc/missingFunctions.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	Tensor::Tensor(const Representation _representation) : Tensor(std::vector<size_t>({}), _representation) { } 
	
	Tensor::Tensor( Tensor&& _other ) : Tensor(_other) {} // TODO improve
	
	Tensor::Tensor(const std::vector<size_t>& _dimensions, const Representation _representation, const Initialisation _init) 
		: dimensions(_dimensions), size(misc::product(dimensions)), representation(_representation)
	{
        REQUIRE(size != 0, "May not create tensors with an dimension == 0.");
		
		if(representation == Representation::Dense) {
			denseData.reset(new value_t[size], internal::array_deleter_vt);
		} else {
			sparseData.reset(new std::map<size_t, value_t>());
		}
    }
    
    Tensor::Tensor(std::vector<size_t>&& _dimensions, const Representation _representation, const Initialisation _init) 
	: dimensions(std::move(_dimensions)), size(misc::product(dimensions)), representation(_representation)
	{
        REQUIRE(size != 0, "May not create tensors with an dimension == 0.");
		
		if(representation == Representation::Dense) {
			denseData.reset(new value_t[size], internal::array_deleter_vt);
		} else {
			sparseData.reset(new std::map<size_t, value_t>());
		}
    }
    
    Tensor::Tensor(std::initializer_list<size_t>&& _dimensions, const Representation _representation, const Initialisation _init) 
	: Tensor(std::vector<size_t>(_dimensions), _representation, _init) {}
    
    
    
    FullTensor Tensor::ones(const std::vector<size_t>& _dimensions) {
		FullTensor ret(_dimensions, DONT_SET_ZERO());
		value_t* const data = ret.data_pointer();
		for(size_t i = 0; i < ret.size; ++i) {
			data[i] = 1.0;
		}
		return ret;
	}
		
	SparseTensor Tensor::identity(const std::vector<size_t>& _dimensions) {
		REQUIRE(_dimensions.size()%2 == 0, "Identity tensor must have even degree, here: " << _dimensions.size());
		const size_t d = _dimensions.size();
		
		SparseTensor ret(_dimensions);
		std::vector<size_t> position(_dimensions.size(), 0);
		
		if(d == 0) {
			ret[position] = 1.0;
		} else {
			bool notMultiLevelBreak = true;
			while(notMultiLevelBreak) {
				ret[position] = 1.0;
				
				position[0]++; position[d/2]++;
				size_t node = 0;
				while(position[node] == std::min(_dimensions[node], _dimensions[d/2+node])) {
					position[node] = position[d/2+node] = 0;
					node++;
					if(node == d/2) {notMultiLevelBreak = false; break;}
					position[node]++; position[d/2+node]++;
				}
			}
		}
		
		return ret;
	}
		
	SparseTensor Tensor::kronecker(const std::vector<size_t>& _dimensions) {
		SparseTensor ret(_dimensions);
		if(_dimensions.empty()) {
			ret[{}] = 1.0;
		} else {
			for(size_t i = 0; i < misc::min(ret.dimensions); ++i) {
				ret[std::vector<size_t>(ret.degree(), i)] = 1.0;
			}
		}
		return ret;
	}
		
	SparseTensor Tensor::dirac(const std::vector<size_t>& _dimensions, const std::vector<size_t>& _position) {
		SparseTensor ret(_dimensions);
		ret[_position] = 1.0;
		return ret;
	}
		
	SparseTensor Tensor::dirac(const std::vector<size_t>& _dimensions, const size_t _position) {
		SparseTensor ret(_dimensions);
		ret[_position] = 1.0;
		return ret;
	}
    
    Tensor::~Tensor() {}
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Information - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	/// @brief Returns whether the current representation is dense.
	bool Tensor::dense() const {
		REQUIRE((representation == Representation::Dense && denseData && !sparseData) || (representation == Representation::Sparse && sparseData && !denseData), "Internal Error: " << bool(representation) << bool(denseData) << bool(sparseData));
		return representation == Representation::Dense;
	}
	
	/// @brief Returns whether the current representation is sparse.
	bool Tensor::sparse() const {
		REQUIRE((representation == Representation::Dense && denseData && !sparseData) || (representation == Representation::Sparse && sparseData && !denseData), "Internal Error: " << bool(representation) << bool(denseData) << bool(sparseData));
		return representation == Representation::Sparse;
	}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
	size_t Tensor::multiIndex_to_position(const std::vector<size_t>& _multiIndex, const std::vector<size_t>& _dimensions) {
		REQUIRE(_multiIndex.size() == _dimensions.size(), "MultiIndex has wrong degree. Given " << _multiIndex.size() << ",  expected " << _dimensions.size());
		
		size_t finalIndex = 0;
		for(size_t i = 0; i < _multiIndex.size(); ++i) {
			REQUIRE(_multiIndex[i] < _dimensions[i], "Index "<< i <<" out of bounds " << _multiIndex[i] << " >=! " << _dimensions[i]);
			finalIndex *= _dimensions[i];
			finalIndex += _multiIndex[i];
		}
		
		return finalIndex;
	}
	
	void Tensor::ensure_own_data() {
		if(dense()) {
			if(!denseData.unique()) {
				value_t* const oldDataPtr = denseData.get();
				denseData.reset(new value_t[size], internal::array_deleter_vt);
				misc::array_copy(denseData.get(), oldDataPtr, size);
			}
		} else {
			if(!sparseData.unique()) {
				sparseData.reset(new std::map<size_t, value_t>(*sparseData));
			}
		}
	}
	
	
	void Tensor::ensure_own_data_no_copy() {
		if(dense()) {
			if(!denseData.unique()) {
				denseData.reset(new value_t[size], internal::array_deleter_vt);
			}
		} else {
			if(!sparseData.unique()) {
				sparseData.reset(new std::map<size_t, value_t>());
			}
		}
	}
	
	void Tensor::apply_factor() {
		if(has_factor()) {
			if(dense()) {
				if(denseData.unique()) {
					misc::array_scale(denseData.get(), factor, size);
				} else {
					value_t* const oldDataPtr = denseData.get();
					denseData.reset(new value_t[size], internal::array_deleter_vt);
					misc::array_scaled_copy(denseData.get(), factor, oldDataPtr, size);
				}
			} else {
				if(!sparseData.unique()) {
					sparseData.reset(new std::map<size_t, value_t>(*sparseData));
				}
				
				for(std::pair<const size_t, value_t>& entry : *sparseData) {
					entry.second *= factor;
				}
			}
			
			factor = 1.0;
		}
	}
	
	void Tensor::ensure_own_data_and_apply_factor() {
		if(has_factor()) {
			apply_factor(); 
		} else {
			ensure_own_data();
		}
	}
	
	
	void Tensor::reset(const std::vector<size_t>&  _newDim, const Representation _representation, const Initialisation _init) {
		const size_t oldDataSize = size;
		change_dimensions(_newDim);
		factor = 1.0;
		
		if(_representation == Representation::Dense) {
			if(dense()) {
				if(oldDataSize != size || !denseData.unique()) {
					denseData.reset(new value_t[size], internal::array_deleter_vt);
				}
			} else {
				sparseData.reset();
				denseData.reset(new value_t[size], internal::array_deleter_vt);
				representation = _representation;
			}
			
			if(_init == Initialisation::Zero) {
				memset(denseData.get(), 0, size*sizeof(value_t));
			}
		} else {
			if(sparse()) {
				if(sparseData.unique()) {
					sparseData->clear();
				} else {
					sparseData.reset(new std::map<size_t, value_t>());
				}
			} else {
				denseData.reset();
				sparseData.reset(new std::map<size_t, value_t>());
				representation = _representation;
			}
		}
	}
	
	void Tensor::reset(const std::vector<size_t>&  _newDim, const Initialisation _init) {
		const size_t oldDataSize = size;
		change_dimensions(_newDim);
		factor = 1.0;
		
		if(dense()) {
			if(oldDataSize != size || !denseData.unique()) {
				denseData.reset(new value_t[size], internal::array_deleter_vt);
			}
			
			if(_init == Initialisation::Zero) {
				memset(denseData.get(), 0, size*sizeof(value_t));
			}
		} else {
			if(sparseData.unique()) {
				sparseData->clear();
			} else {
				sparseData.reset(new std::map<size_t, value_t>());
			}
		}
	}
	
	void Tensor::reset(const std::vector<size_t>& _newDim, const std::shared_ptr<value_t>& _newData) {
		change_dimensions(_newDim);
		factor = 1.0;
		
		if(sparse()) {
			sparseData.reset();
			representation = Representation::Dense;
		}
		
		denseData = _newData;
	}
	
	void Tensor::reset(const std::vector<size_t>& _newDim, std::unique_ptr<value_t[]>&& _newData) {
		change_dimensions(_newDim);
		factor = 1.0;
		
		if(sparse()) {
			sparseData.reset();
			representation = Representation::Dense;
		}
		
		denseData = std::move(_newData);
	}
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
    Tensor& Tensor::operator*=(const value_t _factor) {
        factor *= _factor;
        return *this;
    }
    

    Tensor& Tensor::operator/=(const value_t _divisor) {
        factor /= _divisor;
        return *this;
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Access - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	value_t& Tensor::operator[](const size_t _position) {
		REQUIRE(_position < size, "Position " << _position << " does not exist in Tensor of dimensions " << dimensions);
		
		ensure_own_data_and_apply_factor();
		
		if(dense()) {
			return denseData.get()[_position];
		} else {
			return (*sparseData)[_position];
		}
	}

	value_t Tensor::operator[](const size_t _position) const {
		REQUIRE(_position < size, "Position " << _position << " does not exist in Tensor of dimensions " << dimensions);
		
		if(dense()) {
			return factor*denseData.get()[_position];
		} else {
			const std::map<size_t, value_t>::const_iterator entry = sparseData->find(_position);
			if(entry == sparseData->end()) {
				return 0.0;
			} else {
				return factor*entry->second;
			}
		}
	}

	value_t& Tensor::operator[](const std::vector<size_t>& _positions) {
		return operator[](multiIndex_to_position(_positions, dimensions));
	}
	
	value_t Tensor::operator[](const std::vector<size_t>& _positions) const {
		return operator[](multiIndex_to_position(_positions, dimensions));
	}
	
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    
    /// Indexes the tensor.
    IndexedTensor<Tensor> Tensor::operator()(const std::vector<Index>&  _indices) {
		dense(); 
        return IndexedTensor<Tensor>(this, _indices, false);
    }
    
    /// Indexes the tensor.
    IndexedTensor<Tensor> Tensor::operator()(      std::vector<Index>&& _indices) {
		dense(); 
        return IndexedTensor<Tensor>(this, std::move(_indices), false);
    }
    
    /// Indexes the tensor.
    IndexedTensorReadOnly<Tensor> Tensor::operator()(const std::vector<Index>&  _indices) const {
		dense(); 
        return IndexedTensorReadOnly<Tensor>(this, _indices);
    }
    
    /// Indexes the tensor.
    IndexedTensorReadOnly<Tensor> Tensor::operator()(      std::vector<Index>&& _indices) const {
		dense(); 
        return IndexedTensorReadOnly<Tensor>(this, std::move(_indices));
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    size_t Tensor::degree() const {
        return dimensions.size();
    }
    
    bool Tensor::has_factor() const {
        #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wfloat-equal"
            return (factor != 1.0);
        #pragma GCC diagnostic pop
    }
    
    size_t Tensor::reorder_costs() const {
        return is_sparse() ? 10*count_non_zero_entries() : size;
    }
    
    
    void Tensor::reinterpret_dimensions(const std::vector<size_t>& _newDimensions) {
        REQUIRE(misc::product(_newDimensions) == size, "New dimensions must not change the size of the tensor in reinterpretation: " << misc::product(_newDimensions) << " != " << size);
        dimensions = _newDimensions;
    }
    
    void Tensor::reinterpret_dimensions(      std::vector<size_t>&& _newDimensions) {
        REQUIRE(misc::product(_newDimensions) == size, "New dimensions must not change the size of the tensor in reinterpretation: " << misc::product(_newDimensions) << " != " << size);
        dimensions = std::move(_newDimensions);
    }
    
    void Tensor::reinterpret_dimensions(std::initializer_list<size_t>&& _newDimensions) {
        reinterpret_dimensions(std::vector<size_t>(std::move(_newDimensions)));
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
    void Tensor::assign(const Tensor& _other) {
        dimensions = _other.dimensions;
        size = _other.size;
        factor = _other.factor;
    }
    
    void Tensor::assign(Tensor&& _other) {
        dimensions = std::move(_other.dimensions);
        size = _other.size;
        factor = _other.factor;
    }
    
    void Tensor::change_dimensions(const std::vector<size_t>& _newDimensions) {
        dimensions = _newDimensions;
        size = misc::product(dimensions);
    }
    
    void Tensor::change_dimensions(      std::vector<size_t>&& _newDimensions) {
        dimensions = std::move(_newDimensions);
        size = misc::product(dimensions);
    }
    
    bool approx_equal(const xerus::Tensor& _a, const xerus::Tensor& _b, const xerus::value_t _eps) {
        REQUIRE(_a.dimensions == _b.dimensions, 
                "The dimensions of the compared tensors don't match: " << _a.dimensions <<" vs. " << _b.dimensions << " and " << _a.size << " vs. " << _b.size);
        
		const double avgNorm = (_a.frob_norm() + _b.frob_norm())/2.0;
        if (!_a.is_sparse() && !_b.is_sparse()) { // Special treatment for two fullTensors because blas has better accuracy
            return frob_norm(static_cast<const FullTensor&>(_a) - static_cast<const FullTensor&>(_b)) <= _eps*avgNorm;
		} else if(_a.is_sparse()) {
			FullTensor bma(_b);
			bma -= _a;
			return frob_norm(bma) <= _eps*avgNorm;
		} else if(_b.is_sparse()) {
			FullTensor amb(_a);
			amb -= _b;
			return frob_norm(amb) <= _eps*avgNorm;
        } else { // Special treatment if both are sparse, because better asyptotic is possible.
            return frob_norm(static_cast<const SparseTensor&>(_a) - static_cast<const SparseTensor&>(_b)) <= _eps*avgNorm;
        }
    }
    
bool approx_entrywise_equal(const Tensor& _a, const Tensor& _b, const value_t _eps) {
        REQUIRE(_a.dimensions == _b.dimensions, 
                "The dimensions of the compared tensors don't match: " << _a.dimensions <<" vs. " << _b.dimensions << " and " << _a.size << " vs. " << _b.size);
        
        PA_START;
        if (_a.is_sparse() && _b.is_sparse()) { // Special treatment if both are sparse, because better asyptotic is possible.
            const SparseTensor& A = static_cast<const SparseTensor&>(_a);
            const SparseTensor& B = static_cast<const SparseTensor&>(_b);
            const SparseTensor C = A-B;
            const double factorizedEps = _eps/std::abs(C.factor);
            
            for(const std::pair<size_t, value_t>& entry : *C.sparseData) {
                if(std::abs(entry.second)/(std::abs(A[entry.first])+std::abs(B[entry.first])) > factorizedEps) { return false; }
            }
        } else {
            for(size_t i=0; i < _a.size; ++i) {
                if(std::abs(_a[i]-_b[i])/(std::abs(_a[i])+std::abs(_b[i])) > _eps) { return false; }
            }
        }
        PA_END("Approx entrywise equal", "", misc::to_string(_a.size));
        return true;
    }
}
