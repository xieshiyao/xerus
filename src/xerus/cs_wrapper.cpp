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
 * @brief Implementation of the suitesparse wrapper functions.
 */

#include <xerus/cs_wrapper.h>
#include <xerus/index.h>
#include <xerus/tensor.h>
#include <xerus/misc/performanceAnalysis.h>

namespace xerus {
    bool sparse_result(const size_t _lhsDim, const size_t _midDim, const size_t _rhsDim, const size_t _lhsEntries, const size_t _rhsEntries) {
        const size_t lhsSize = _lhsDim*_midDim;
        const size_t rhsSize = _midDim*_rhsDim;
        const size_t upperBound = std::min(_lhsDim*_rhsEntries, _rhsDim*_lhsEntries);
        const size_t finalSize = _lhsDim*_rhsDim;
        return finalSize > 100*upperBound || (lhsSize > 100*_lhsEntries && rhsSize > 100*_rhsEntries);
    }
    
    
    CsUniquePtr create_cs(const size_t _m, const size_t _n, const size_t _N) {
        REQUIRE(_m < std::numeric_limits<int>::max() && _n < std::numeric_limits<int>::max() && _N < std::numeric_limits<int>::max(), "Sparse Tensor is to large for SuiteSparse (" << _m << " x " << _n << ", " << _N << ")");
        return CsUniquePtr(cs_spalloc((int) _m, (int) _n, (int) _N, 1, 0), &cs_spfree);
    }
    
    
    // Converts an Indexed SparseTensor and an given matrification to the CSparse sparse matrix format
    CsUniquePtr to_cs_format(IndexedTensorReadOnly<Tensor>&& _tensor, const std::vector<Index>& _lhsIndices, const std::vector<Index>& _rhsIndices) {
        REQUIRE(_tensor.tensorObjectReadOnly->is_sparse(), "Only sparse Tensors can be converted to CS format.");
        
        size_t m = 1;
        size_t n = 1;
        const std::vector<Index> indices = _tensor.get_assigned_indices();
        for(size_t i = 0; i < indices.size(); ++i) {
            if(misc::contains(_lhsIndices, indices[i])) {
                REQUIRE(indices[i].open(), "Internal Error.");
                m *= indices[i].dimension();
            } else if(misc::contains(_rhsIndices, indices[i])) {
                REQUIRE(indices[i].open(), "Internal Error.");
                n *= indices[i].dimension();
            }   
        }
        
        std::vector<Index> inverseIndexOrder(_rhsIndices);
        inverseIndexOrder.insert(inverseIndexOrder.end(), _lhsIndices.begin(), _lhsIndices.end());
        
        Tensor reorderedTensor( Tensor::Representation::Sparse );
        reorderedTensor(inverseIndexOrder) = std::move(_tensor);
        
        REQUIRE(reorderedTensor.size == m*n, "Internal Error");
        
        // We want A (m x n) in compressed coloum storage. We have reordered the tensor to A^T (n x m) and 
        // transform it to compressed row storage. This results in A in compressed coloum storage as demanded.
        
        CsUniquePtr cs_format = create_cs(m, n, reorderedTensor.sparseData->size());
        
        int entryPos = 0;
        int currRow = -1;
        cs_format->i[0] = 0;
        
        for(const std::pair<size_t, value_t>& entry : *reorderedTensor.sparseData.get()) {
            cs_format->x[entryPos] = entry.second;
            cs_format->i[entryPos] = (int) (entry.first%m);
            while(currRow < (int) (entry.first/m)) {
                cs_format->p[++currRow] = entryPos;
            }
            entryPos++;
        }
        REQUIRE(currRow < (int) n && entryPos == (int) reorderedTensor.sparseData->size(), "Internal Error " << currRow << ", " << (int) n << " | " << entryPos << ", " << (int) reorderedTensor.sparseData->size());
        while(currRow < (int) n) {
            cs_format->p[++currRow] = entryPos;
        }
        
        return cs_format;
    }
    
    
    Tensor from_cs_format(const CsUniquePtr& _cs_format, const std::vector<size_t>& _dimensions) {
        REQUIRE(_cs_format, "NullPtr cannot be converted to Tensor.");
        
        Tensor reconstructedTensor(_dimensions, Tensor::Representation::Sparse);
        
        for(int i = 0; i < _cs_format->n; ++i) {
            for(int j = _cs_format->p[i]; j < _cs_format->p[i+1]; ++j) {
                IF_CHECK( auto ret = ) reconstructedTensor.sparseData->insert(std::pair<size_t, value_t>(_cs_format->i[j]*_cs_format->n+i, _cs_format->x[j]));
                REQUIRE(ret.second, "Internal Error");
            }
        }
        
        return reconstructedTensor;
    }
    
    
    CsUniquePtr matrix_matrix_product(const CsUniquePtr& _lhs, const CsUniquePtr& _rhs) {
		PA_START;
		CsUniquePtr result(cs_multiply(_lhs.get(), _rhs.get()), &cs_spfree);
		PA_END("Sparse BLAS", "Matrix-Matrix-Multiplication", misc::to_string(_lhs->m)+"x"+misc::to_string(_lhs->n)+" * "+misc::to_string(_rhs->m)+"x"+misc::to_string(_rhs->n));
        return CsUniquePtr(cs_multiply(_lhs.get(), _rhs.get()), &cs_spfree);
    }
    
    
    void print_cs(const CsUniquePtr& _cs_format) {
        std::cout << "Sparse Matrix parameters: n = " << _cs_format->n << ", m = " << _cs_format->m << ", Max Entries = " << _cs_format->nzmax << std::endl;
        std::cout << "Values =  {";
        std::cout << _cs_format->x[0];
        for(int i = 1; i < _cs_format->nzmax; ++i) {
            std::cout << ", " << _cs_format->x[i];
        }
        std::cout << "}" << std::endl;
        
        std::cout << "Row idx = {";
        std::cout << _cs_format->i[0];
        for(int i = 1; i < _cs_format->nzmax; ++i) {
            std::cout << ", " << _cs_format->i[i];
        }
        std::cout << "}" << std::endl;
        
        std::cout << "Col Pos = {";
        std::cout << _cs_format->p[0];
        for(int i = 1; i < _cs_format->n + 1; ++i) {
            std::cout << ", " << _cs_format->p[i];
        }
        std::cout << "}" << std::endl;
    }
}
