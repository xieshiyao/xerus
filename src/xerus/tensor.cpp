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
	
	Tensor::Tensor(const Representation _representation) : Tensor(DimensionTuple({}), _representation) { } 
	
	Tensor::Tensor(const DimensionTuple& _dimensions, const Representation _representation, const Initialisation _init) 
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
    
    Tensor::Tensor(DimensionTuple&& _dimensions, const Representation _representation, const Initialisation _init) 
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
    
    Tensor::Tensor(const DimensionTuple& _dimensions, std::unique_ptr<value_t[]>&& _data)
	: dimensions(_dimensions), size(misc::product(dimensions)), representation(Representation::Dense), denseData(std::move(_data)) { }
    
    
    Tensor::Tensor(const DimensionTuple& _dimensions, const std::function<value_t()>& _f) : Tensor(_dimensions, Representation::Dense, Initialisation::None) {
		value_t* const realData = denseData.get();
		for (size_t i=0; i < size; ++i) {
			realData[i] = _f();
		}
	}
	
	Tensor::Tensor(const DimensionTuple& _dimensions, const std::function<value_t(const size_t)>& _f) : Tensor(_dimensions, Representation::Dense, Initialisation::None) {
		value_t* const realData = denseData.get();
		for (size_t i=0; i < size; ++i) {
			realData[i] = _f(i);
		}
	}
	
	Tensor::Tensor(const DimensionTuple& _dimensions, const std::function<value_t(const MultiIndex&)>& _f) : Tensor(_dimensions, Representation::Dense, Initialisation::None) {
		value_t* const realData = denseData.get();
		MultiIndex multIdx(degree(), 0);
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
	
	Tensor::Tensor(const DimensionTuple& _dimensions, const size_t _N, const std::function<std::pair<size_t, value_t>(size_t, size_t)>& _f) : Tensor(_dimensions, Representation::Sparse, Initialisation::Zero) {
		REQUIRE(_N <= size, "Cannot create more non zero entries that the dimension of the Tensor.");
		for (size_t i=0; i < _N; ++i) {
			std::pair<size_t, value_t> entry = _f(i, size);
			REQUIRE(entry.first < size, "Postion is out of bounds " << entry.first);
			REQUIRE(!misc::contains(*sparseData, entry.first), "Allready contained " << entry.first);
			sparseData->insert(std::move(entry));
		}
	}
    
    
    
    
    Tensor Tensor::ones(const DimensionTuple& _dimensions) {
		Tensor ret(_dimensions, Representation::Dense, Initialisation::None);
		value_t* const data = ret.get_dense_data();
		for(size_t i = 0; i < ret.size; ++i) {
			data[i] = 1.0;
		}
		return ret;
	}
		
	Tensor Tensor::identity(const DimensionTuple& _dimensions) {
		REQUIRE(_dimensions.size()%2 == 0, "Identity tensor must have even degree, here: " << _dimensions.size());
		const size_t d = _dimensions.size();
		
		Tensor ret(_dimensions, Representation::Sparse);
		MultiIndex position(_dimensions.size(), 0);
		
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
		
	Tensor Tensor::kronecker(const DimensionTuple& _dimensions) {
		Tensor ret(_dimensions, Representation::Sparse);
		if(_dimensions.empty()) {
			ret[{}] = 1.0;
		} else {
			for(size_t i = 0; i < misc::min(ret.dimensions); ++i) {
				ret[MultiIndex(ret.degree(), i)] = 1.0;
			}
		}
		return ret;
	}
		
	Tensor Tensor::dirac(const DimensionTuple& _dimensions, const MultiIndex& _position) {
		Tensor ret(_dimensions, Representation::Sparse);
		ret[_position] = 1.0;
		return ret;
	}
		
	Tensor Tensor::dirac(const DimensionTuple& _dimensions, const size_t _position) {
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

	value_t& Tensor::operator[](const MultiIndex& _positions) {
		return operator[](multiIndex_to_position(_positions, dimensions));
	}
	
	value_t Tensor::operator[](const MultiIndex& _positions) const {
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
	
	void Tensor::reset(const DimensionTuple&  _newDim, const Representation _representation, const Initialisation _init) {
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
	
	void Tensor::reset(const DimensionTuple&  _newDim, const Initialisation _init) {
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
	
	void Tensor::reset(const DimensionTuple& _newDim, const std::shared_ptr<value_t>& _newData) {
		dimensions = _newDim;
		size = misc::product(dimensions);
		factor = 1.0;
		
		if(is_sparse()) {
			sparseData.reset();
			representation = Representation::Dense;
		}
		
		denseData = _newData;
	}
	
	void Tensor::reset(const DimensionTuple& _newDim, std::unique_ptr<value_t[]>&& _newData) {
		dimensions = _newDim;
		size = misc::product(dimensions);
		factor = 1.0;
		
		if(is_sparse()) {
			sparseData.reset();
			representation = Representation::Dense;
		}
		
		denseData = std::move(_newData);
	}
	
	void Tensor::reinterpret_dimensions(const DimensionTuple& _newDimensions) {
		REQUIRE(misc::product(_newDimensions) == size, "New dimensions must not change the size of the tensor in reinterpretation: " << misc::product(_newDimensions) << " != " << size);
		dimensions = _newDimensions;
	}
	
	void Tensor::reinterpret_dimensions(      DimensionTuple&& _newDimensions) {
		REQUIRE(misc::product(_newDimensions) == size, "New dimensions must not change the size of the tensor in reinterpretation: " << misc::product(_newDimensions) << " != " << size);
		dimensions = std::move(_newDimensions);
	}
	
	void Tensor::resize_dimension(const size_t _dimPos, const size_t _newDim, size_t _cutPos) {
		REQUIRE(_dimPos < degree(), "Can't resize dimension " << _dimPos << " as the tensor is only order " << degree());

		if (dimensions[_dimPos] == _newDim) { return; }  // Trivial case: Nothing to do

		const size_t oldDim = dimensions[_dimPos];
		_cutPos = std::min(_cutPos, oldDim);
		
		REQUIRE(_newDim > 0, "Dimension must be larger than 0! Is " << _newDim);
		REQUIRE(_newDim > oldDim || _cutPos >= oldDim -_newDim, "Cannot remove " << oldDim -_newDim << " slates starting (exclusivly) at position " << _cutPos);
		
		const size_t dimStepSize = misc::product(dimensions, _dimPos+1, degree());
		const size_t oldStepSize = oldDim*dimStepSize;
		const size_t newStepSize = _newDim*dimStepSize;
		const size_t blockCount = size / oldStepSize; //  == misc::product(dimensions, 0, _n);
		const size_t newsize = blockCount*newStepSize;
		
		if(is_dense()) {
			std::unique_ptr<value_t[]> tmpData(new value_t[newsize]);
			
			if (_newDim > oldDim) { // Add new slates
				const size_t insertBlockSize = (_newDim-oldDim)*dimStepSize;
				if(_cutPos == 0) {
					for ( size_t i = 0; i < blockCount; ++i ) {
						value_t* const currData = tmpData.get()+i*newStepSize;
						misc::array_set_zero(currData, insertBlockSize);
						misc::array_copy(currData+insertBlockSize, denseData.get()+i*oldStepSize, oldStepSize);
					}
				} else if(_cutPos == oldDim) {
					for ( size_t i = 0; i < blockCount; ++i ) {
						value_t* const currData = tmpData.get()+i*newStepSize;
						misc::array_copy(currData, denseData.get()+i*oldStepSize, oldStepSize);
						misc::array_set_zero(currData+oldStepSize, insertBlockSize);
					}
				} else {
					const size_t preBlockSize = _cutPos*dimStepSize;
					const size_t postBlockSize = (oldDim-_cutPos)*dimStepSize;
					for (size_t i = 0; i < blockCount; ++i) {
						value_t* const currData = tmpData.get()+i*newStepSize;
						misc::array_copy(currData, denseData.get()+i*oldStepSize, preBlockSize);
						misc::array_set_zero(currData+preBlockSize, insertBlockSize);
						misc::array_copy(currData+preBlockSize+insertBlockSize, denseData.get()+i*oldStepSize+preBlockSize, postBlockSize);
					}
				}
			} else { // Remove slates
				if (_cutPos < oldDim) {
					const size_t removedSlates = (oldDim-_newDim);
					const size_t removedBlockSize = removedSlates*dimStepSize;
					const size_t preBlockSize = (_cutPos-removedSlates)*dimStepSize;
					const size_t postBlockSize = (oldDim-_cutPos)*dimStepSize;
					
					REQUIRE(removedBlockSize+preBlockSize+postBlockSize == oldStepSize && preBlockSize+postBlockSize == newStepSize, "IE");
					
					for (size_t i = 0; i < blockCount; ++i) {
						value_t* const currData = tmpData.get()+i*newStepSize;
						misc::array_copy(currData, denseData.get()+i*oldStepSize, preBlockSize);
						misc::array_copy(currData+preBlockSize, denseData.get()+i*oldStepSize+preBlockSize+removedBlockSize, postBlockSize);
					}
				} else {
					for (size_t i = 0; i < blockCount; ++i) {
						misc::array_copy(tmpData.get()+i*newStepSize, denseData.get()+i*oldStepSize, newStepSize);
					}
				}
			}
			denseData.reset(tmpData.release(), internal::array_deleter_vt);
		
		} else {
			std::unique_ptr<std::map<size_t, value_t>> tmpData(new std::map<size_t, value_t>());
			
			if (_newDim > oldDim) { // Add new slates
				const size_t slatesAdded = _newDim-oldDim;
				for(const std::pair<size_t, value_t>& entry : *sparseData.get()) {
					// Decode the position as i*oldStepSize + j*DimStepSize + k
					const size_t k = entry.first%dimStepSize;
					const size_t j = (entry.first%oldStepSize)/dimStepSize;
					const size_t i = entry.first/oldStepSize;
					
					if(j < _cutPos) { // Entry remains at current position j
						tmpData->emplace(i*newStepSize+j*dimStepSize+k, entry.second);
					} else { // Entry moves to position j+slatesAdded
						tmpData->emplace(i*newStepSize+(j+slatesAdded)*dimStepSize+k, entry.second);
					}
				}
			} else { // Remove slates
				const size_t slatesRemoved = oldDim-_newDim;
				for(const std::pair<size_t, value_t>& entry : *sparseData.get()) {
					// Decode the position as i*oldStepSize + j*DimStepSize + k
					const size_t k = entry.first%dimStepSize;
					const size_t j = (entry.first%oldStepSize)/dimStepSize;
					const size_t i = entry.first/oldStepSize;
					
					if(j < _cutPos-slatesRemoved) { // Entry remains at current position j
						tmpData->emplace(i*newStepSize+j*dimStepSize+k, entry.second);
					} else if(j >= _cutPos) { // Entry advances in position j
						tmpData->emplace(i*newStepSize+(j-slatesRemoved)*dimStepSize+k, entry.second);
					}
				}
			}
			sparseData.reset(tmpData.release());
		}
		
		dimensions[_dimPos] = _newDim;
		size = newsize;
	}
	
	void Tensor::fix_slate(const size_t _dimPos, const size_t _slatePosition) {
		REQUIRE(_slatePosition < dimensions[_dimPos], "The given slatePosition must be smaller than the corresponding dimension. Here " << _slatePosition << " >= " << dimensions[_dimPos]);
		
		const size_t stepCount = misc::product(dimensions, 0, _dimPos);
		const size_t blockSize = misc::product(dimensions, _dimPos+1, degree());
		const size_t stepSize = dimensions[_dimPos]*blockSize;
		const size_t slateOffset = _slatePosition*blockSize;
		
		if(is_dense()) {
			std::unique_ptr<value_t[]> tmpData(new value_t[stepCount*blockSize]);
			
			// Copy data
			for(size_t i = 0; i < stepCount; ++i) {
				misc::array_copy(tmpData.get()+i*blockSize, denseData.get()+i*stepSize+slateOffset, blockSize);
			}
			
			denseData.reset(tmpData.release(), &internal::array_deleter_vt);
		} else {
			std::unique_ptr<std::map<size_t, value_t>> tmpData(new std::map<size_t, value_t>());
			
			for(const std::pair<size_t, value_t>& entry : *sparseData.get()) {
				// Decode the position as i*stepSize + j*blockSize + k
				const size_t k = entry.first%blockSize;
				const size_t j = (entry.first%stepSize)/blockSize;
				const size_t i = entry.first/stepSize;
				
				if(j == _slatePosition) {
					tmpData->emplace(i*blockSize+k, entry.second);
				}
			}
			
			sparseData.reset(tmpData.release());
		}
		
		
		// Adjust dimensions
		dimensions.erase(dimensions.begin()+long(_dimPos));
		size = stepCount*blockSize;
	}
	
	void Tensor::remove_slate(const size_t _indexNb, const size_t _pos) {
		REQUIRE(_indexNb < degree(), "");
		REQUIRE(_pos < dimensions[_indexNb], _pos << " " << dimensions[_indexNb]);
		REQUIRE(dimensions[_indexNb] > 1, "");
		
		resize_dimension(_indexNb, dimensions[_indexNb]-1, _pos+1);
	}
	
	
	void Tensor::modify_diag_elements(const std::function<void(value_t&)>& _f) {
		ensure_own_data_and_apply_factor();
		
		if(degree() == 0) {
			_f(at(0)); 
		} else {
			size_t stepSize = 1;
			for(size_t i = 1; i < degree(); ++i) {
				stepSize *= dimensions[i];
				stepSize += 1;
			}
			
			const size_t numDiags = misc::min(dimensions);
			
			for(size_t i = 0; i < numDiags; ++i){
				_f(at(i*stepSize));
			}
		}
	}
	
	void Tensor::modify_diag_elements(const std::function<void(value_t&, const size_t)>& _f) {
		ensure_own_data_and_apply_factor();
		
		if(degree() == 0) {
			_f(at(0), 0); 
		} else {
			size_t stepSize = 1;
			for(size_t i = 1; i < degree(); ++i) {
				stepSize *= dimensions[i];
				stepSize += 1;
			}
			
			const size_t numDiags = misc::min(dimensions);
			
			for(size_t i = 0; i < numDiags; ++i){
				_f(at(i*stepSize), i);
			}
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
				if(misc::hard_not_equal(val, 0.0)) {
					at(i) = val;
				} else if( misc::hard_not_equal(oldVal, 0.0)) {
					IF_CHECK( size_t succ =) sparseData->erase(i);
					REQUIRE(succ == 1, "Internal Error");
				}
			}
		}
	}

	void Tensor::modify_elements(const std::function<void(value_t&, const size_t)>& _f) {
		ensure_own_data_and_apply_factor();
		if(is_dense()) {
			for(size_t i = 0; i < size; ++i) { _f(at(i), i); }
		} else {
			for(size_t i = 0; i < size; ++i) {
				value_t val = cat(i);
				const value_t oldVal = val;
				_f(val, i);
				if(misc::hard_not_equal(val, 0.0)) {
					at(i) = val;
				} else if( misc::hard_not_equal(oldVal, 0.0)) {
					IF_CHECK( size_t succ =) sparseData->erase(i);
					REQUIRE(succ == 1, "Internal Error");
				}
			}
		}
	}
	
	void Tensor::modify_elements(const std::function<void(value_t&, const MultiIndex&)>& _f) {
		ensure_own_data_and_apply_factor();
		
		MultiIndex multIdx(degree(), 0);
		size_t idx = 0;
		while (true) {
			if(is_dense()) {
				_f(at(idx), multIdx);
			} else {
				value_t val = cat(idx);
				const value_t oldVal = val;
				_f(val, multIdx);
				if(misc::hard_not_equal(val, 0.0)) {
					at(idx) = val;
				} else if( misc::hard_not_equal(oldVal, 0.0)) {
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
					sparseData->emplace_hint(sparseData->end(), i, factor*denseData.get()[i]);
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
	
	size_t Tensor::multiIndex_to_position(const MultiIndex& _multiIndex, const DimensionTuple& _dimensions) {
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
		REQUIRE(_numIndices <= _lhs.degree() && _numIndices <= _rhs.degree(), "Cannot contract more indices than each involdes has.");
		
		const size_t lhsRemainOrder = _lhs.degree() - _numIndices;
		const size_t lhsRemainStart = _lhsTrans ? _numIndices : 0;
		const size_t lhsContractStart = _lhsTrans ? 0 : lhsRemainOrder;
		const size_t lhsRemainEnd = lhsRemainStart + lhsRemainOrder;
		
		const size_t rhsRemainOrder = _rhs.degree() - _numIndices;
		const size_t rhsRemainStart = _rhsTrans ? 0 : _numIndices;
		const size_t rhsContractStart = _rhsTrans ? rhsRemainOrder : 0;
		const size_t rhsRemainEnd = rhsRemainStart + rhsRemainOrder;
		
		REQUIRE(std::equal(_lhs.dimensions.begin() + lhsContractStart, _lhs.dimensions.begin() + lhsContractStart + _numIndices, _rhs.dimensions.begin() + rhsContractStart), "Dimensions of the be contracted indices do not coincide.");
		
		
		// Check whether _result has to be reset
		if( _result.degree() != lhsRemainOrder + rhsRemainOrder 
			|| !std::equal(_lhs.dimensions.begin() + lhsRemainStart, _lhs.dimensions.begin() + lhsRemainEnd, _result.dimensions.begin())
			|| !std::equal(_rhs.dimensions.begin() + rhsRemainStart, _rhs.dimensions.begin() + rhsRemainEnd, _result.dimensions.begin())
		) {
			Tensor::DimensionTuple resultDim;
			resultDim.reserve(lhsRemainOrder + rhsRemainOrder);
			resultDim.insert(resultDim.end(), _lhs.dimensions.begin() + lhsRemainStart, _lhs.dimensions.begin() + lhsRemainEnd);
			resultDim.insert(resultDim.end(), _rhs.dimensions.begin() + rhsRemainStart, _rhs.dimensions.begin() + rhsRemainEnd);
			_result.reset(std::move(resultDim), Tensor::Initialisation::None);
		}
		
		const size_t leftDim = misc::product(_lhs.dimensions, lhsRemainStart, lhsRemainEnd);
		const size_t midDim = misc::product(_lhs.dimensions, lhsContractStart, lhsContractStart+_numIndices);
		const size_t rightDim = misc::product(_rhs.dimensions, rhsRemainStart, rhsRemainEnd);
		
		if(!_lhs.is_sparse() && !_rhs.is_sparse()) { // Full * Full => Full
			blasWrapper::matrix_matrix_product(_result.override_dense_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
											   _lhs.get_unsanitized_dense_data(), _lhsTrans, midDim, 
											   _rhs.get_unsanitized_dense_data(), _rhsTrans);
			
		} else if(_lhs.is_sparse() && _rhs.is_sparse() ) { // Sparse * Sparse => Sparse
			internal::matrix_matrix_product(_result.override_sparse_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_sparse_data(), _lhsTrans, midDim,
								  _rhs.get_unsanitized_sparse_data(), _rhsTrans);
			
		} else if(_lhs.is_sparse() && !_rhs.is_sparse() && !_result.is_sparse()) { // Sparse * Full => Full
			matrix_matrix_product(_result.override_dense_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_sparse_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_dense_data(), _rhsTrans);
			
		} else if(!_lhs.is_sparse() && _rhs.is_sparse() && !_result.is_sparse()) { // Full * Sparse => Full
			matrix_matrix_product(_result.override_dense_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_dense_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_sparse_data(), _rhsTrans);
			
		} else if(_lhs.is_sparse() && !_rhs.is_sparse() && _result.is_sparse()) { // Sparse * Full => Sparse
			matrix_matrix_product(_result.override_sparse_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_sparse_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_dense_data(), _rhsTrans);
			
		} else if(!_lhs.is_sparse() && _rhs.is_sparse() && _result.is_sparse()) { // Full * Sparse => Sparse
			matrix_matrix_product(_result.override_sparse_data(), leftDim, rightDim, _lhs.factor*_rhs.factor, 
								  _lhs.get_unsanitized_dense_data(), _lhsTrans, midDim, 
								  _rhs.get_unsanitized_sparse_data(), _rhsTrans);
		}
	}
	
	_inline_ std::tuple<size_t, size_t, size_t> prepare_factorization_output(Tensor& _lhs, Tensor& _rhs, const Tensor& _input, const size_t _splitPos) {
		REQUIRE(_splitPos < _input.degree(), "_splitPos must smaller than the degree.");
		
		const size_t lhsSize = misc::product(_input.dimensions, 0, _splitPos);
		const size_t rhsSize = misc::product(_input.dimensions, _splitPos, _input.degree());
		const size_t rank = std::min(lhsSize, rhsSize);
		
		_lhs.factor = 1.0;
		_rhs.factor = 1.0;
		
		if( !_lhs.is_dense()
			|| _lhs.degree() != _splitPos+1
			|| rank != _lhs.dimensions.back() 
			|| !std::equal(_input.dimensions.begin(), _input.dimensions.begin() + _splitPos, _lhs.dimensions.begin())) 
		{
			Tensor::DimensionTuple newDimU;
			newDimU.insert(newDimU.end(), _input.dimensions.begin(), _input.dimensions.begin() + _splitPos);
			newDimU.push_back(rank);
			_lhs.reset(newDimU, Tensor::Representation::Dense, Tensor::Initialisation::None);
		}
		
		if(!_rhs.is_dense()
			|| _rhs.degree() != _input.degree()-_splitPos+1 
			|| rank != _rhs.dimensions.front() 
			|| !std::equal(_input.dimensions.begin()+_splitPos, _input.dimensions.end(), _rhs.dimensions.begin()+1)) 
		{
			Tensor::DimensionTuple newDimVt;
			newDimVt.push_back(rank);
			newDimVt.insert(newDimVt.end(), _input.dimensions.begin() + _splitPos, _input.dimensions.end());
			_rhs.reset(newDimVt, Tensor::Representation::Dense, Tensor::Initialisation::None);
		}
		
		return std::make_tuple(lhsSize, rhsSize, rank);
	}
	
	void calculate_svd(Tensor& _U, Tensor& _S, Tensor& _Vt, const Tensor& _input, const size_t _splitPos, const size_t _maxRank, const value_t _eps) {
		REQUIRE(_eps < 1, "Epsilon must be smaller than one.");
		
		size_t lhsSize, rhsSize, rank;
		std::tie(lhsSize, rhsSize, rank) = prepare_factorization_output(_U, _Vt, _input, _splitPos);
		
		std::unique_ptr<value_t[]> tmpS(new value_t[rank]);
		
		// Calculate the actual SVD
		if(_input.is_sparse()) {
			LOG(fatal, "Sparse SVD not yet implemented.");
		} else {
			blasWrapper::svd(_U.override_dense_data(), tmpS.get(), _Vt.override_dense_data(), _input.get_unsanitized_dense_data(), lhsSize, rhsSize);
		}
		
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
		_S.reset(Tensor::DimensionTuple(2, rank), Tensor::Representation::Sparse);
		for(size_t i = 0; i < rank; ++i) {
			_S[i*rank+i] = std::abs(_input.factor)*tmpS[i];
		}
		
		// Account for a negative factor
		if(_input.factor < 0.0) {
			_Vt *= -1;
		}
		
		_U.resize_dimension(_U.degree()-1, rank);
		_Vt.resize_dimension(0, rank);
	}
	
	
	void calculate_qr(Tensor& _Q, Tensor& _R, const Tensor& _input, const size_t _splitPos) {
		size_t lhsSize, rhsSize, rank;
		std::tie(lhsSize, rhsSize, rank) = prepare_factorization_output(_Q, _R, _input, _splitPos);
		
		if(_input.is_sparse()) {
			LOG(fatal, "Sparse SVD not yet implemented.");
		} else {
			blasWrapper::qr(_Q.override_dense_data(), _R.override_dense_data(), _input.get_unsanitized_dense_data(), lhsSize, rhsSize);
		}
		
		_R.factor = _input.factor;
	}
	
	
	void calculate_rq(Tensor& _R, Tensor& _Q, const Tensor& _input, const size_t _splitPos) {
		size_t lhsSize, rhsSize, rank;
		std::tie(lhsSize, rhsSize, rank) = prepare_factorization_output(_R, _Q, _input, _splitPos);
		
		if(_input.is_sparse()) {
			LOG(fatal, "Sparse SVD not yet implemented.");
		} else {
			blasWrapper::rq(_R.override_dense_data(), _Q.override_dense_data(), _input.get_unsanitized_dense_data(), lhsSize, rhsSize);
		}
		
		_R.factor = _input.factor;
	}
	
	
	void calculate_qc(Tensor& _Q, Tensor& _C, const Tensor& _input, const size_t _splitPos, const value_t _eps) {
		LOG(fatal, "Not yet Implemented."); // TODO
	}
	
	
	void calculate_cq(Tensor& _C, Tensor& _Q, const Tensor& _input, const size_t _splitPos, const value_t _eps) {
		LOG(fatal, "Not yet Implemented."); // TODO
	}
	
	
	void solve_least_squares(Tensor& _x, const Tensor& _A, const Tensor& _b) {
		LOG(fatal, "Not yet Implemented."); // TODO
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
