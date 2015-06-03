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

#include <xerus/tensor.h>
#include <xerus/misc/check.h>
#include <xerus/misc/missingFunctions.h>
#include <xerus/fullTensor.h>
#include <xerus/sparseTensor.h>

namespace xerus {
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Constructors - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    Tensor::Tensor() : size(1), factor(1.0) {}

    Tensor::Tensor(const  Tensor&  _other ) : dimensions(_other.dimensions), size(_other.size), factor(_other.factor) { }
    
    Tensor::Tensor(       Tensor&& _other ) : dimensions(_other.dimensions), size(_other.size), factor(_other.factor) { }
    
    Tensor::Tensor(const std::vector<size_t>& _dimensions, const value_t _factor) : dimensions(_dimensions), size(misc::product(dimensions)), factor(_factor) {
        REQUIRE(size != 0, "May not create tensors with an dimension == 0.");
    }
    
    Tensor::Tensor(std::vector<size_t>&& _dimensions, const value_t _factor) : dimensions(std::move(_dimensions)), size(misc::product(dimensions)), factor(_factor) {
        REQUIRE(size != 0, "May not create tensors with an dimension == 0.");
    }
    
    Tensor::Tensor(std::initializer_list<size_t>&& _dimensions, const value_t _factor) : dimensions(std::move(_dimensions)), size(misc::product(dimensions)), factor(_factor) {
        REQUIRE(size != 0, "May not create tensors with an dimension == 0.");
    }
    
    Tensor::~Tensor() {}
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Basic arithmetics - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
        
    Tensor& Tensor::operator*=(const value_t _factor) {
        factor *= _factor;
        return *this;
    }
    

    Tensor& Tensor::operator/=(const value_t _divisor) {
        factor /= _divisor;
        return *this;
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - Indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    
    /// Indexes the tensor.
    IndexedTensor<Tensor> Tensor::operator()(const std::vector<Index>&  _indices) {
        return IndexedTensor<Tensor>(this, _indices, false);
    }
    
    /// Indexes the tensor.
    IndexedTensor<Tensor> Tensor::operator()(      std::vector<Index>&& _indices) {
        return IndexedTensor<Tensor>(this, std::move(_indices), false);
    }
    
    /// Indexes the tensor.
    IndexedTensorReadOnly<Tensor> Tensor::operator()(const std::vector<Index>&  _indices) const {
        return IndexedTensorReadOnly<Tensor>(this, _indices);
    }
    
    /// Indexes the tensor.
    IndexedTensorReadOnly<Tensor> Tensor::operator()(      std::vector<Index>&& _indices) const {
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
    
    void Tensor::reinterpret_dimensions( std::initializer_list<size_t> _newDimensions) {
        reinterpret_dimensions(std::vector<size_t>(_newDimensions));
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
        
        
        if (!_a.is_sparse() && !_b.is_sparse()) { // Special treatment for two fullTensors because blas has better accuracy
            return frob_norm(static_cast<const FullTensor&>(_a) - static_cast<const FullTensor&>(_b))/(_a.frob_norm() + _b.frob_norm()) < _eps;
        } else if (_a.is_sparse() && _b.is_sparse()) { // Special treatment if both are sparse, because better asyptotic is possible.
            return frob_norm(static_cast<const SparseTensor&>(_a) - static_cast<const SparseTensor&>(_b))/(_a.frob_norm() + _b.frob_norm()) < _eps;
        } else {
            PA_START;
            double sqrNorm = 0;
            for(size_t i=0; i < _a.size; ++i) {
                sqrNorm += misc::sqr(_a[i]-_b[i]);
            }
            PA_END("Mixed Norm calculation", "", misc::to_string(_a.size));
            return sqrt(sqrNorm)/(_a.frob_norm() + _b.frob_norm()) < _eps;
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
            
            for(const std::pair<size_t, value_t>& entry : *C.entries) {
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
