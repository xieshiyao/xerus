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
#include <xerus/selectedFunctions.h>
 
#include <xerus/blasLapackWrapper.h>
#include <xerus/cs_wrapper.h>
#include <xerus/sparseTimesFullContraction.h>

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	Tensor::Tensor(const Representation _representation) : Tensor(std::vector<size_t>({}), _representation) { } 
	
	Tensor::Tensor(const std::vector<size_t>& _dimensions, const Representation _representation, const Initialisation _init) 
		: dimensions(_dimensions), size(misc::product(dimensions)), representation(_representation)
	{
        REQUIRE(size != 0, "May not create tensors with an dimension == 0.");
		
		if(representation == Representation::Dense) {
			denseData.reset(new value_t[size], internal::array_deleter_vt);
			if(_init == Initialisation::Zero) {
				misc::array_set_zero(denseData.get(), size);
			}
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
			if(_init == Initialisation::Zero) {
				misc::array_set_zero(denseData.get(), size);
			}
		} else {
			sparseData.reset(new std::map<size_t, value_t>());
		}
    }
    
    Tensor::Tensor(const std::vector<size_t>& _dimensions, std::unique_ptr<value_t[]>&& _data)
	: dimensions(_dimensions), size(misc::product(dimensions)), representation(Representation::Dense), denseData(std::move(_data)) { }
    
    
    Tensor::Tensor(const std::vector<size_t>& _dimensions, const std::function<value_t()>& _f) : Tensor(_dimensions, Representation::Dense, Initialisation::None) {
		value_t* const realData = denseData.get();
		for (size_t i=0; i < size; ++i) {
			realData[i] = _f();
		}
	}
	
	Tensor::Tensor(const std::vector<size_t>& _dimensions, const std::function<value_t(const size_t)>& _f) : Tensor(_dimensions, Representation::Dense, Initialisation::None) {
		value_t* const realData = denseData.get();
		for (size_t i=0; i < size; ++i) {
			realData[i] = _f(i);
		}
	}
	
	Tensor::Tensor(const std::vector<size_t>& _dimensions, const std::function<value_t(const std::vector<size_t>&)>& _f) : Tensor(_dimensions, Representation::Dense, Initialisation::None) {
		value_t* const realData = denseData.get();
		std::vector<size_t> multIdx(degree(), 0);
		size_t idx = 0;
		while (true) {
			realData[idx] = _f(multIdx);
			// Increasing indices
			idx++;
			size_t changingIndex = degree()-1;
			multIdx[changingIndex]++;
			while(multIdx[changingIndex] == dimensions[changingIndex]) {
				multIdx[changingIndex] = 0;
				changingIndex--;
				// Return on overflow 
				if(changingIndex >= degree()) { return; }
				multIdx[changingIndex]++;
			}
		}
	}
	
	Tensor::Tensor(const std::vector<size_t>& _dimensions, std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f, const size_t _N) : Tensor(_dimensions, Representation::Sparse, Initialisation::Zero) {
		REQUIRE(_N <= size, "Cannot create more non zero entries that the dimension of the Tensor.");
		for (size_t i=0; i < _N; ++i) {
			std::pair<size_t, value_t> entry = _f(i, size);
			REQUIRE(entry.first < size, "Postion is out of bounds " << entry.first);
			REQUIRE(!misc::contains(*sparseData, entry.first), "Allready contained " << entry.first);
			sparseData->insert(std::move(entry));
		}
	}
    
    
    
    
    Tensor Tensor::ones(const std::vector<size_t>& _dimensions) {
		Tensor ret(_dimensions, Representation::Dense, Initialisation::None);
		value_t* const data = ret.get_dense_data();
		for(size_t i = 0; i < ret.size; ++i) {
			data[i] = 1.0;
		}
		return ret;
	}
		
	Tensor Tensor::identity(const std::vector<size_t>& _dimensions) {
		REQUIRE(_dimensions.size()%2 == 0, "Identity tensor must have even degree, here: " << _dimensions.size());
		const size_t d = _dimensions.size();
		
		Tensor ret(_dimensions, Representation::Sparse);
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
		
	Tensor Tensor::kronecker(const std::vector<size_t>& _dimensions) {
		Tensor ret(_dimensions, Representation::Sparse);
		if(_dimensions.empty()) {
			ret[{}] = 1.0;
		} else {
			for(size_t i = 0; i < misc::min(ret.dimensions); ++i) {
				ret[std::vector<size_t>(ret.degree(), i)] = 1.0;
			}
		}
		return ret;
	}
		
	Tensor Tensor::dirac(const std::vector<size_t>& _dimensions, const std::vector<size_t>& _position) {
		Tensor ret(_dimensions, Representation::Sparse);
		ret[_position] = 1.0;
		return ret;
	}
		
	Tensor Tensor::dirac(const std::vector<size_t>& _dimensions, const size_t _position) {
		Tensor ret(_dimensions, Representation::Sparse);
		ret[_position] = 1.0;
		return ret;
	}
    
    
    Tensor Tensor::dense_copy() const {
		Tensor ret(*this);
		ret.use_dense_representation();
		return ret;
	}
	
	
	Tensor Tensor::sparse_copy() const {
		Tensor ret(*this);
		ret.use_sparse_representation();
		return ret;
	}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Standard operators - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	Tensor& Tensor::operator=( Tensor&& _other) {
		std::swap(dimensions, _other.dimensions);
		std::swap(size, _other.size);
		std::swap(representation, _other.representation);
		factor = _other.factor;
		std::swap(denseData, _other.denseData);
		std::swap(sparseData, _other.sparseData);
		return *this;
	}
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Information - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	size_t Tensor::degree() const {
		return dimensions.size();
	}
	
	bool Tensor::has_factor() const {
		#pragma GCC diagnostic push
			#pragma GCC diagnostic ignored "-Wfloat-equal"
			return (factor != 1.0);
		#pragma GCC diagnostic pop
	}
	
	bool Tensor::is_dense() const {
		REQUIRE((representation == Representation::Dense && denseData && !sparseData) || (representation == Representation::Sparse && sparseData && !denseData), "Internal Error: " << bool(representation) << bool(denseData) << bool(sparseData));
		return representation == Representation::Dense;
	}
	
	bool Tensor::is_sparse() const {
		REQUIRE((representation == Representation::Dense && denseData && !sparseData) || (representation == Representation::Sparse && sparseData && !denseData), "Internal Error: " << bool(representation) << bool(denseData) << bool(sparseData));
		return representation == Representation::Sparse;
	}
	
	size_t Tensor::sparsity() const {
		REQUIRE(is_sparse(), "IE");
		return sparseData->size();
	}
	
	
	size_t Tensor::count_non_zero_entries(const value_t _eps) const {
		if(is_dense()) {
			size_t count = 0;
			for(size_t i = 0; i < size; ++i) {
				if(std::abs(denseData.get()[i]) > _eps ) { count++; }
			}
			return count;
		} else {
			size_t count = 0;
			for(const std::pair<size_t, value_t>& entry : *sparseData) {
				if(std::abs(entry.second) > _eps) { count++; } 
			}
			return count;
		}
	}
	
	bool Tensor::all_entries_valid() const {
		if(is_dense()) {
			for(size_t i = 0; i < size; ++i) {
				if(!std::isfinite(denseData.get()[i])) { return false; } 
			}
		} else {
			for(const std::pair<size_t, value_t>& entry : *sparseData) {
				if(!std::isfinite(entry.second)) {return false; } 
			}
		}
		return true;
	}
	
	size_t Tensor::reorder_costs() const {
		return is_sparse() ? 10*sparsity() : size;
	}
	
	value_t Tensor::frob_norm() const {
		if(is_dense()) {
			return std::abs(factor)*blasWrapper::two_norm(denseData.get(), size);
		} else {
			value_t norm = 0;
			for(const std::pair<size_t, value_t>& entry : *sparseData) {
				norm += misc::sqr(entry.second);
			}
			return std::abs(factor)*sqrt(norm);
		}
	}
	


    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	Tensor& Tensor::operator+=(const Tensor& _other) {
		plus_minus_equal<1>(*this, _other);
		return *this;
	}
	
	Tensor& Tensor::operator-=(const Tensor& _other) {
		plus_minus_equal<-1>(*this, _other);
		return *this;
	}
	
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
		
		if(is_dense()) {
			return denseData.get()[_position];
		} else {
			return (*sparseData)[_position];
		}
	}

	value_t Tensor::operator[](const size_t _position) const {
		REQUIRE(_position < size, "Position " << _position << " does not exist in Tensor of dimensions " << dimensions);
		
		if(is_dense()) {
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
	
	
	value_t& Tensor::at(const size_t _position) {
		REQUIRE(_position < size, "Position " << _position << " does not exist in Tensor of dimensions " << dimensions);
		REQUIRE(!has_factor(), "at() must not be called if there is a factor.");
		REQUIRE((is_dense() && denseData.unique()) || (is_sparse() && sparseData.unique()) , "Data must be unique to call at().");
		
		if(is_dense()) {
			return denseData.get()[_position];
		} else {
			return (*sparseData)[_position];
		}
	}
	
	value_t Tensor::cat(const size_t _position) const {
		REQUIRE(_position < size, "Position " << _position << " does not exist in Tensor of dimensions " << dimensions);
		REQUIRE(!has_factor(), "at() must not be called if there is a factor.");
		
		if(is_dense()) {
			return denseData.get()[_position];
		} else {
			const std::map<size_t, value_t>::const_iterator entry = sparseData->find(_position);
			if(entry == sparseData->end()) {
				return 0.0;
			} else {
				return entry->second;
			}
		}
	}
	
	
	
	value_t* Tensor::get_dense_data() {
		use_dense_representation();
		ensure_own_data_and_apply_factor();
		return denseData.get();
	}
	
	value_t* Tensor::get_unsanitized_dense_data() {
		REQUIRE(is_dense(), "Unsanitized dense data requested, but representation is not dense!");
		return denseData.get();
	}
	
	const value_t* Tensor::get_unsanitized_dense_data() const  {
		REQUIRE(is_dense(), "Unsanitized dense data requested, but representation is not dense!");
		return denseData.get();
	}
	
	value_t* Tensor::override_dense_data()  {
		factor = 1.0;
		if(!denseData.unique()) {
			sparseData.reset();
			denseData.reset(new value_t[size], internal::array_deleter_vt);
			representation = Representation::Dense;
		}
		return denseData.get();
	}
	
	const std::shared_ptr<value_t>& Tensor::get_internal_dense_data() {
		REQUIRE(is_dense(), "Internal dense data requested, but representation is not dense!");
		return denseData;
	}
	
	
	std::map<size_t, value_t>& Tensor::get_sparse_data() {
		CHECK(is_sparse(), warning, "Request for sparse data although the Tensor is not sparse.");
		use_sparse_representation();
		ensure_own_data_and_apply_factor();
		return *sparseData.get();
	}
	
	std::map<size_t, value_t>& Tensor::get_unsanitized_sparse_data() {
		REQUIRE(is_sparse(), "Unsanitized sparse data requested, but representation is not sparse!");
		return *sparseData.get();
	}
	
	const std::map<size_t, value_t>& Tensor::get_unsanitized_sparse_data() const  {
		REQUIRE(is_sparse(), "Unsanitized sparse data requested, but representation is not sparse!");
		return *sparseData.get();
	}
	
	std::map<size_t, value_t>& Tensor::override_sparse_data() {
		factor = 1.0;
		if(sparseData.unique()) {
			REQUIRE(is_sparse(), "Internal Error");
			sparseData->clear();
		} else {
			denseData.reset();
			sparseData.reset(new std::map<size_t, value_t>());
			representation = Representation::Sparse;
		}
		
		return *sparseData.get();
	}
	
	const std::shared_ptr<std::map<size_t, value_t>>& Tensor::get_internal_sparse_data() {
		REQUIRE(is_sparse(), "Internal sparse data requested, but representation is not sparse!");
		return sparseData;
	}

	
	
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// Indexes the tensor.
    IndexedTensor<Tensor> Tensor::operator()(const std::vector<Index>&  _indices) {
		is_dense(); 
        return IndexedTensor<Tensor>(this, _indices, false);
    }
    
    /// Indexes the tensor.
    IndexedTensor<Tensor> Tensor::operator()(      std::vector<Index>&& _indices) {
		is_dense(); 
        return IndexedTensor<Tensor>(this, std::move(_indices), false);
    }
    
    /// Indexes the tensor.
    IndexedTensorReadOnly<Tensor> Tensor::operator()(const std::vector<Index>&  _indices) const {
		is_dense(); 
        return IndexedTensorReadOnly<Tensor>(this, _indices);
    }
    
    /// Indexes the tensor.
    IndexedTensorReadOnly<Tensor> Tensor::operator()(      std::vector<Index>&& _indices) const {
		is_dense(); 
        return IndexedTensorReadOnly<Tensor>(this, std::move(_indices));
    }
    
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Modififiers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	void Tensor::reset(const std::vector<size_t>&  _newDim, const Representation _representation, const Initialisation _init) {
		const size_t oldDataSize = size;
		dimensions = _newDim;
		size = misc::product(dimensions);
		factor = 1.0;
		
		if(_representation == Representation::Dense) {
			if(is_dense()) {
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
			if(is_sparse()) {
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
		dimensions = _newDim;
		size = misc::product(dimensions);
		factor = 1.0;
		
		if(is_dense()) {
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
		dimensions = _newDim;
		size = misc::product(dimensions);
		factor = 1.0;
		
		if(is_sparse()) {
			sparseData.reset();
			representation = Representation::Dense;
		}
		
		denseData = _newData;
	}
	
	void Tensor::reset(const std::vector<size_t>& _newDim, std::unique_ptr<value_t[]>&& _newData) {
		dimensions = _newDim;
		size = misc::product(dimensions);
		factor = 1.0;
		
		if(is_sparse()) {
			sparseData.reset();
			representation = Representation::Dense;
		}
		
		denseData = std::move(_newData);
	}
	
	void Tensor::reinterpret_dimensions(const std::vector<size_t>& _newDimensions) {
		REQUIRE(misc::product(_newDimensions) == size, "New dimensions must not change the size of the tensor in reinterpretation: " << misc::product(_newDimensions) << " != " << size);
		dimensions = _newDimensions;
	}
	
	void Tensor::reinterpret_dimensions(      std::vector<size_t>&& _newDimensions) {
		REQUIRE(misc::product(_newDimensions) == size, "New dimensions must not change the size of the tensor in reinterpretation: " << misc::product(_newDimensions) << " != " << size);
		dimensions = std::move(_newDimensions);
	}
	
	void Tensor::resize_dimension(const size_t _n, const size_t _newDim, const size_t _cutPos) {
		REQUIRE(_n < degree(), "Can't resize dimension " << _n << " as the tensor is only order " << degree());
		REQUIRE(_newDim > 0, "Dimension must be larger than 0! Is " << _newDim);
		REQUIRE(_newDim >= dimensions[_n] || _cutPos >= dimensions[_n] -_newDim -1, "Cannot remove " << dimensions[_n] -_newDim << " slates starting at cut position " << _cutPos);
		
		if (dimensions[_n] == _newDim) { return; }  // Trivial case: Nothing to do
		
		const size_t dimStepSize = misc::product(dimensions, _n+1, degree());
		const size_t oldStepSize = dimensions[_n]*dimStepSize;
		const size_t newStepSize = _newDim*dimStepSize;
		const size_t blockCount = size / oldStepSize; //  == misc::product(dimensions, 0, _n);
		const size_t newsize = blockCount*newStepSize;
		
		REQUIRE(newsize == size/dimensions[_n] * _newDim, 
					dimensions[_n] << " " << _newDim << " " << oldStepSize << " " << newStepSize 
					<< " " << blockCount << size << " " << newsize);
		
		if(is_dense()) {
			value_t* const tmp = new value_t[newsize];
			const size_t cutOffset = _cutPos*dimStepSize;
			
			if (_newDim > dimensions[_n]) {
				if (_cutPos < dimensions[_n]) {
					const size_t numInsert = (newStepSize-oldStepSize);
					for (size_t i = 0; i < blockCount; ++i) {
						misc::array_copy(tmp+i*newStepSize, denseData.get()+i*oldStepSize, cutOffset);
						misc::array_set_zero(tmp+i*newStepSize+cutOffset, numInsert);
						misc::array_copy(tmp+i*newStepSize+cutOffset+numInsert, denseData.get()+i*oldStepSize+cutOffset, (oldStepSize-cutOffset));
					}
				} else {
					for (size_t i = 0; i < blockCount; ++i) {
						misc::array_copy(tmp+i*newStepSize, denseData.get()+i*oldStepSize, oldStepSize);
						misc::array_set_zero(tmp+i*newStepSize+oldStepSize, (newStepSize-oldStepSize));
					}
				}
			} else { // _newDim < dimensions[_n]
				if (_cutPos < dimensions[_n]) {
					const size_t removedSize = (dimensions[_n]-_newDim)*dimStepSize;
					const size_t firstSize = cutOffset+dimStepSize-removedSize;
					const size_t secondSize = oldStepSize-cutOffset-dimStepSize;
					for (size_t i = 0; i < blockCount; ++i) {
						misc::array_copy(tmp+i*newStepSize, denseData.get()+i*oldStepSize, firstSize);
						misc::array_copy(tmp+i*newStepSize+firstSize, denseData.get()+i*oldStepSize+cutOffset+dimStepSize, secondSize);
					}
				} else {
					for (size_t i = 0; i < blockCount; ++i) {
						misc::array_copy(tmp+i*newStepSize, denseData.get()+i*oldStepSize, newStepSize);
					}
				}
			}
			
			denseData.reset(tmp, internal::array_deleter_vt);
			
		} else {
			REQUIRE(_cutPos >= dimensions[_n], "Support for dimension resize at arbitary position not yet ther efor Sparse Tensor");
			
			std::map<size_t, value_t>* newData = new std::map<size_t, value_t>();
			
			for(const std::pair<size_t, value_t>& entry : *sparseData.get()) {
				// Decode the position as i*oldStepSize + j*DimStepSize + k
				const size_t k = entry.first%dimStepSize;
				const size_t j = (entry.first/dimStepSize)%dimensions[_n];
				const size_t i = (entry.first/dimStepSize)/dimensions[_n];
				
				// Only add if within the new range
				if(j <= _newDim) {
					newData->emplace(i*newStepSize+j*dimStepSize+k, entry.second);
				}
			}
			
			sparseData.reset(newData);
		}
		
		dimensions[_n] = _newDim;
		size = newsize;
	}
	
	void Tensor::fix_slate(const size_t _dimension, const size_t _slatePosition) {
		REQUIRE(_slatePosition < dimensions[_dimension], "The given slatePosition must be smaller than the corresponding dimension. Here " << _slatePosition << " >= " << dimensions[_dimension]);
		
		if(is_dense()) {
			size_t stepCount = 1, blockSize = 1;
			for(size_t i = 0; i < _dimension; ++i) { stepCount *= dimensions[i]; }
			for(size_t i = _dimension+1; i < dimensions.size(); ++i) { blockSize *= dimensions[i]; }
			
			const size_t stepSize = dimensions[_dimension]*blockSize;
			size_t inputPosition = _slatePosition*blockSize;
			
			value_t * const newData = new value_t[stepCount*blockSize];
			
			// Copy data
			for(size_t i = 0; i < stepCount; ++i) {
				misc::array_copy(newData+i*blockSize, denseData.get()+inputPosition, blockSize);
				inputPosition += stepSize;
			}
			
			// Set data
			denseData.reset(newData, &internal::array_deleter_vt);
			
			// Adjust dimensions
			dimensions.erase(dimensions.begin()+long(_dimension));
			size = stepCount*blockSize;
		} else {
			// TODO implement!
			LOG(fatal, "fix_slate is not yet implemented for sparse representations."); 
		}
	}
	
	void Tensor::remove_slate(const size_t _indexNb, const size_t _pos) {
		REQUIRE(_indexNb < degree(), "");
		REQUIRE(_pos < dimensions[_indexNb], _pos << " " << dimensions[_indexNb]);
		REQUIRE(dimensions[_indexNb] > 1, "");
		
		resize_dimension(_indexNb, dimensions[_indexNb]-1, _pos);
	}
	
	
	//TODO Allow more 2d
	void Tensor::modify_diag_elements(const std::function<void(value_t&)>& _f) {
		REQUIRE(degree() == 2, "Diagonal elements are only well defined if degree equals two. Here: "  << degree());
		ensure_own_data_and_apply_factor();
		const size_t numDiags = std::min(dimensions[0], dimensions[1]);
		const size_t N = dimensions[1];
		for(size_t i = 0; i < numDiags; ++i){
			_f(at(i+i*N));
		}
	}
	
	//TODO Allow more 2d
	void Tensor::modify_diag_elements(const std::function<void(value_t&, const size_t)>& _f) {
		REQUIRE(degree() == 2, "Diagonal elements are only well defined if degree equals two. Here: "  << degree());
		ensure_own_data_and_apply_factor();
		const size_t numDiags = std::min(dimensions[0], dimensions[1]);
		const size_t N = dimensions[1];
		for(size_t i = 0; i < numDiags; ++i){
			_f(at(i+i*N), i);
		}
	}
	
	void Tensor::modify_elements(const std::function<void(value_t&)>& _f) {
		ensure_own_data_and_apply_factor();
		if(is_dense()) {
			for(size_t i = 0; i < size; ++i) { _f(at(i)); }
		} else {
			for(size_t i = 0; i < size; ++i) {
				value_t val = cat(i);
				const value_t oldVal = val;
				_f(val);
				if(val != 0.0) {
					at(i) = val;
				} else if( oldVal != 0.0) {
					IF_CHECK( size_t succ =) sparseData->erase(i);
					REQUIRE(succ == 1, "Internal Error");
				}
			}
		}
	}
	
	#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 7) || defined(__clang__)
		void Tensor::modify_elements(const std::function<void(value_t&, const size_t)>& _f) {
			ensure_own_data_and_apply_factor();
			if(is_dense()) {
				for(size_t i = 0; i < size; ++i) { _f(at(i), i); }
			} else {
				for(size_t i = 0; i < size; ++i) {
					value_t val = cat(i);
					const value_t oldVal = val;
					_f(val, i);
					if(val != 0.0) {
						at(i) = val;
					} else if( oldVal != 0.0) {
						IF_CHECK( size_t succ =) sparseData->erase(i);
						REQUIRE(succ == 1, "Internal Error");
					}
				}
			}
		}
		
		void Tensor::modify_elements(const std::function<void(value_t&, const std::vector<size_t>&)>& _f) {
			ensure_own_data_and_apply_factor();
			
			std::vector<size_t> multIdx(degree(), 0);
			size_t idx = 0;
			while (true) {
				if(is_dense()) {
					_f(at(idx), multIdx);
				} else {
					value_t val = cat(idx);
					const value_t oldVal = val;
					_f(val, multIdx);
					if(val != 0.0) {
						at(idx) = val;
					} else if( oldVal != 0.0) {
						IF_CHECK( size_t succ =) sparseData->erase(idx);
						REQUIRE(succ == 1, "Internal Error");
					}
				}
				
				// Increase indices
				idx++;
				size_t changingIndex = degree()-1;
				multIdx[changingIndex]++;
				while(multIdx[changingIndex] == dimensions[changingIndex]) {
					multIdx[changingIndex] = 0;
					changingIndex--;
					// Return on overflow 
					if(changingIndex >= degree()) { return; }
					multIdx[changingIndex]++;
				}
			}
		}
	#endif
    
    void Tensor::use_dense_representation() {
		if(is_sparse()) {
			denseData.reset(new value_t[size], internal::array_deleter_vt);
			misc::array_set_zero(denseData.get(), size);
			for(const std::pair<size_t, value_t>& entry : *sparseData) {
				denseData.get()[entry.first] = entry.second;
			}
			sparseData.reset();
			representation = Representation::Dense;
		}
	}
	
	void Tensor::use_sparse_representation(const value_t _eps) {
		if(is_dense()) {
			sparseData.reset(new std::map<size_t, value_t>());
			for(size_t i = 0; i < size; ++i) {
				if(std::abs(factor*denseData.get()[i]) >= _eps) {
					sparseData->insert({i, factor*denseData.get()[i]}); // TODO use emplace_hint
				}
			}
			
			denseData.reset();
			representation = Representation::Sparse;
		}
	}
	
	
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Miscellaneous - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
	std::string Tensor::to_string() const {
		if (degree() == 0) return xerus::misc::to_string(operator[](0));
		
		std::string result;
		for (size_t i = 0; i < size; ++i) {
			result += xerus::misc::to_string(operator[](i)) + " ";
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
	
	
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Internal Helper functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	template<int sign>
	void Tensor::plus_minus_equal(Tensor& _me, const Tensor& _other) {
		REQUIRE(_me.dimensions == _other.dimensions, "In Tensor sum the dimensions must coincide.");
		
		if(_me.is_dense()) {
			_me.ensure_own_data_and_apply_factor();
			if(_other.is_dense()) {
				misc::array_add(_me.denseData.get(), sign*_other.factor, _other.denseData.get(), _me.size);
			} else {
				add_sparse_to_full(_me.denseData, sign*_other.factor, _other.sparseData);
			}
		} else {
			if(_other.is_dense()) {
				_me.denseData.reset(new value_t[_me.size], internal::array_deleter_vt);
				misc::array_scaled_copy(_me.denseData.get(), sign*_other.factor, _other.denseData.get(), _me.size);
				add_sparse_to_full(_me.denseData, _me.factor, _me.sparseData);
				_me.factor = 1.0;
				_me.representation = Representation::Dense;
				_me.sparseData.reset();
			} else {
				_me.ensure_own_data_and_apply_factor();
				add_sparse_to_sparse(_me.sparseData, sign*_other.factor, _other.sparseData);
			}
		}
	}
	
	void Tensor::add_sparse_to_full(const std::shared_ptr<value_t>& _denseData, const value_t _factor, const std::shared_ptr<const std::map<size_t, value_t>>& _sparseData) {
		for(const std::pair<size_t, value_t>& entry : *_sparseData) {
			_denseData.get()[entry.first] += _factor*entry.second;
		}
	}
	
	
	void Tensor::add_sparse_to_sparse(const std::shared_ptr<std::map<size_t, value_t>>& _sum, const value_t _factor, const std::shared_ptr<const std::map<size_t, value_t>>& _summand) {
		for(const std::pair<size_t, value_t>& entry : *_summand) {
			std::pair<std::map<size_t, value_t>::iterator, bool> result = _sum->emplace(entry.first, _factor*entry.second);
			if(!result.second) {
				result.first->second += _factor*entry.second;
			}
		}
	}
	
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
		if(is_dense()) {
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
		if(is_dense()) {
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
			if(is_dense()) {
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
	
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - External functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	
	/*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	
	Tensor operator+(const Tensor& _lhs, const Tensor& _rhs) {
		Tensor result(_lhs);
		result += _rhs;
		return result;
	}
	
	Tensor operator-(const Tensor& _lhs, const Tensor& _rhs) {
		Tensor result(_lhs);
		result -= _rhs;
		return result;
	}
	
	Tensor operator*(const value_t _factor, const Tensor& _tensor) {
		Tensor result(_tensor);
		result *= _factor;
		return result;
	}
	
	Tensor operator/(const Tensor& _tensor, const value_t _divisor) {
		Tensor result(_tensor);
		result /= _divisor;
		return result;
	}
	
	void contract(Tensor& _result, const Tensor& _lhs, const bool _lhsTrans, const Tensor& _rhs, const bool _rhsTrans, const size_t _numIndices) {
		bool sparseLhs = _lhs.is_sparse();
		bool sparseRhs = _rhs.is_sparse();
		bool sparseResult = _result.is_sparse();
		
		REQUIRE(_numIndices <= _lhs.degree() && _numIndices <= _rhs.degree(), "_numIndices must not be larger than one of the degrees. here " << _numIndices << " > " << _lhs.degree() << "/" << _rhs.degree());
		const size_t leftDim = _lhsTrans ? misc::product(_lhs.dimensions, _numIndices, _lhs.degree()) : misc::product(_lhs.dimensions, 0, _lhs.degree() - _numIndices);
		const size_t midDim = _lhsTrans ? misc::product(_lhs.dimensions, 0, _numIndices) : misc::product(_lhs.dimensions, _lhs.degree()-_numIndices, _lhs.degree());
		const size_t rightDim = _rhsTrans ? misc::product(_rhs.dimensions, 0, _rhs.degree()-_numIndices) : misc::product(_rhs.dimensions, _numIndices, _rhs.degree());
		
		IF_CHECK(
			std::vector<size_t> resultDims;
			if(_lhsTrans) {
				for(size_t i = _numIndices; i < _lhs.degree(); ++i) {
					resultDims.emplace_back(_lhs.dimensions[i]);
				}
			} else {
				for(size_t i = 0; i < _lhs.degree()-_numIndices; ++i) {
					resultDims.emplace_back(_lhs.dimensions[i]);
				}
			}
			if(_rhsTrans) {
				for(size_t i = 0; i < _rhs.degree()-_numIndices; ++i) {
					resultDims.emplace_back(_rhs.dimensions[i]);
				}
			} else {
				for(size_t i = _numIndices; i < _rhs.degree(); ++i) {
					resultDims.emplace_back(_rhs.dimensions[i]);
				}
			}
			REQUIRE(resultDims == _result.dimensions, "The given results has wrong dimensions " << _result.dimensions << " should be " << resultDims);
			
			for(size_t i = 0; i < _numIndices; ++i) {
				REQUIRE((_lhsTrans ? _lhs.dimensions[i] : _lhs.dimensions[_lhs.degree()-_numIndices+i]) == (_rhsTrans ? _rhs.dimensions[_rhs.degree()-_numIndices+i] : _rhs.dimensions[i]), "Dimensions of the be contracted indices do not coincide.");
			}
		)
		
		
		if(!sparseLhs && !sparseRhs && !sparseResult) { // Full * Full => Full
			blasWrapper::matrix_matrix_product(_result.override_dense_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
											   _lhs.get_unsanitized_dense_data(), _lhsTrans, midDim, 
											   _rhs.get_unsanitized_dense_data(), _rhsTrans);
		} else if(sparseLhs && !sparseRhs && !sparseResult) { // Sparse * Full => Full
			matrix_matrix_product(_result.override_dense_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_sparse_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_dense_data(), _rhsTrans);
		} else if(!sparseLhs && sparseRhs && !sparseResult) { // Full * Sparse => Full
			matrix_matrix_product(_result.override_dense_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_dense_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_sparse_data(), _rhsTrans);
		} else if(sparseLhs && !sparseRhs && sparseResult) { // Sparse * Full => Sparse
			matrix_matrix_product(_result.override_sparse_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_sparse_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_dense_data(), _rhsTrans);
		} else if(!sparseLhs && sparseRhs && sparseResult) { // Full * Sparse => Sparse
			matrix_matrix_product(_result.override_sparse_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_dense_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_sparse_data(), _rhsTrans);
		} else {
			LOG(fatal, "Sparse times Sparse contraction not implemented (as Tensro function)");
		}
	}
	
	
	void svd(Tensor& _U, Tensor& _S, Tensor& _Vt, const Tensor& _input, const size_t _splitPos, const size_t _maxRank, const value_t _eps) {
		REQUIRE(_eps < 1, "Epsilon must be smaller than one.");
		
		const size_t lhsSize = misc::product(_input.dimensions, 0, _splitPos);
		const size_t rhsSize = misc::product(_input.dimensions, _splitPos, _input.degree());
		size_t rank = std::min(lhsSize, rhsSize);
		
		std::unique_ptr<value_t[]> tmpS(new value_t[rank]);
		
		
		REQUIRE(_U.size == lhsSize*rank, "Err");
		REQUIRE(_Vt.size == rhsSize*rank, "Err");
		
		// Calculate the actual SVD
		if(_input.is_sparse()) {
			LOG(fatal, "Sparse SVD not yet implemented.");
		} else {
			blasWrapper::svd(_U.override_dense_data(), tmpS.get(), _Vt.override_dense_data(), _input.get_unsanitized_dense_data(), lhsSize, rhsSize);
		}
		
		// Apply factor to the diagonal matrix
		misc::array_scale(tmpS.get(), std::abs(_input.factor), rank);
		
		// Account for hard threshold
		if(_maxRank != 0) {
			rank = std::min(rank, _maxRank);
		}
		
		// Find rank due to the Epsilon (NOTE the scaling factor can be ignored, as it does not change the ratios).
		for(size_t j = 1; j < rank; ++j) {
			if (tmpS[j] <= _eps*tmpS[0]) {
				rank = j;
				break;
			}
		}
		
		// Create tensor from diagonal values
		_S.reset(std::vector<size_t>(2, rank), Tensor::Representation::Sparse);
		for(size_t i = 0; i < rank; ++i) {
			_S[i*rank+i] = std::abs(_input.factor)*tmpS[i];
		}
		
		// Account for a negative factor
		if(_input.factor < 0.0) {
			_U *= -1;
		}
		
		_U.resize_dimension(_U.degree()-1, rank);
		_Vt.resize_dimension(0, rank);
	}
	
	Tensor entrywise_product(const Tensor &_A, const Tensor &_B) {
		REQUIRE(_A.dimensions == _B.dimensions, "Entrywise product ill-defined for non-equal dimensions.");
		if(_A.is_dense() && _B.is_dense()) {
			Tensor result(_A);
			result.ensure_own_data();
			value_t* const dataPtrA = result.get_unsanitized_dense_data();
			const value_t* const dataPtrB = _B.get_unsanitized_dense_data();
			result *= _B.factor;
			for(size_t i = 0; i < result.size; ++i) {
				dataPtrA[i] *= dataPtrB[i];
			}
			return result;
		} else if(_A.is_sparse()) {
			Tensor result(_A);
			
			for(std::pair<const size_t, value_t>& entry : result.get_unsanitized_sparse_data()) {
				entry.second *= _B[entry.first];
			}
			return result;
		} else { // _B.is_sparse()
			Tensor result(_B);
			
			for(std::pair<const size_t, value_t>& entry : result.get_unsanitized_sparse_data()) {
				entry.second *= _A[entry.first];
			}
			return result;
		}
	}
    
    bool approx_equal(const xerus::Tensor& _a, const xerus::Tensor& _b, const xerus::value_t _eps) {
        REQUIRE(_a.dimensions == _b.dimensions, 
                "The dimensions of the compared tensors don't match: " << _a.dimensions <<" vs. " << _b.dimensions << " and " << _a.size << " vs. " << _b.size);
        
		return frob_norm(_a - _b) <= _eps*(_a.frob_norm() + _b.frob_norm())/2.0;
    }
    
	bool approx_entrywise_equal(const Tensor& _a, const Tensor& _b, const value_t _eps) {
        REQUIRE(_a.dimensions == _b.dimensions, 
                "The dimensions of the compared tensors don't match: " << _a.dimensions <<" vs. " << _b.dimensions << " and " << _a.size << " vs. " << _b.size);
        
        if (_a.is_sparse() && _b.is_sparse()) { // Special treatment if both are sparse, because better asyptotic is possible.
			for(const std::pair<size_t, value_t>& entry : _a.get_unsanitized_sparse_data()) {
				if(!approx_equal(_a.factor*entry.second, _b[entry.first], _eps)) { return false; }
			}
			
			for(const std::pair<size_t, value_t>& entry : _b.get_unsanitized_sparse_data()) {
				if(!approx_equal(_a[entry.first], _b.factor*entry.second, _eps)) { return false; }
			}
        } else {
            for(size_t i=0; i < _a.size; ++i) {
                if(!approx_equal(_a[i], _b[i], _eps)) { return false; }
            }
        }
        return true;
    }
    
    bool approx_entrywise_equal(const xerus::Tensor& _tensor, const std::vector<value_t>& _values, const xerus::value_t _eps) {
		if(_tensor.size != _values.size()) { return false; }
		for(size_t i = 0; i < _tensor.size; ++i) {
			if(!approx_equal(_tensor[i], _values[i], _eps)) { return false; }
		}
		return true;
	}
}
