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

#include "xerus.h"

namespace xerus {
    CsUniquePtr create_cs(const size_t _m, const size_t _n, const size_t _N) {
        REQUIRE(_m < std::numeric_limits<int>::max() && _n < std::numeric_limits<int>::max() && _N < std::numeric_limits<int>::max(), "Sparse Tensor is to large for SuiteSparse");
        return CsUniquePtr(cs_spalloc((int) _m, (int) _n, (int) _N, 1, 0), &cs_spfree);
    }
    
    
    // Converts an Indexed SparseTensor and an given matrification to the CSparse sparse matrix format
    CsUniquePtr to_cs_format(const IndexedTensorReadOnly<Tensor>& _tensor, const std::vector<Index>& _lhsIndices, const std::vector<Index>& _rhsIndices) {
        REQUIRE(_tensor.tensorObjectReadOnly->is_sparse(), "Only sparse Tensors can be converted to CS format.");
        
        size_t rowDim = 1;
        size_t colDim = 1;
        const AssignedIndices assIndices = _tensor.assign_indices();
        for(size_t i = 0; i < assIndices.numIndices; ++i) {
            if(contains(_lhsIndices, assIndices.indices[i])) {
                REQUIRE(assIndices.indexOpen[i], "Internal Error.");
                rowDim *= assIndices.indexDimensions[i];
            } else if(contains(_rhsIndices, assIndices.indices[i])) {
                REQUIRE(assIndices.indexOpen[i], "Internal Error.");
                colDim *= assIndices.indexDimensions[i];
            }   
        }
        
        std::vector<Index> inverseIndexOrder(_rhsIndices);
        inverseIndexOrder.insert(inverseIndexOrder.end(), _lhsIndices.begin(), _lhsIndices.end());
        
        SparseTensor reorderedTensor;
        reorderedTensor(inverseIndexOrder) = _tensor;
        
        REQUIRE(reorderedTensor.size == rowDim*colDim, "Internal Error");
        
        CsUniquePtr cs_format = create_cs(rowDim, colDim, reorderedTensor.entries->size());
        
        int entryPos = 0;
        int currRow = -1;
        cs_format->i[0] = 0;
        
        for(const std::pair<size_t, value_t>& entry : *reorderedTensor.entries.get() ) {
            cs_format->x[entryPos] = entry.second;
            cs_format->i[entryPos] = (int) (entry.first%rowDim);
            while(currRow < (int) (entry.first/rowDim)) {
                cs_format->p[++currRow] = entryPos;
            }
            entryPos++;
        }
        REQUIRE(currRow <= (int) rowDim && entryPos == (int) reorderedTensor.entries->size(), "Internal Error " << currRow << ", " << (int) rowDim << " | " << entryPos << ", " << (int) reorderedTensor.entries->size());
        while(currRow < (int) colDim) {
            cs_format->p[++currRow] = entryPos;
        }
            
        cs_format->p[currRow] = entryPos;
        
        return cs_format;
    }
    
    
    SparseTensor from_cs_format(const CsUniquePtr& _cs_format, const std::vector<size_t>& _dimensions) {
        SparseTensor reconstructedTensor(_dimensions);
        
        for(int i = 0; i < _cs_format->n; ++i) {
            for(int j = _cs_format->p[i]; j < _cs_format->p[i+1]; ++j) {
                auto ret = reconstructedTensor.entries->insert(std::pair<size_t, value_t>(_cs_format->i[j]*_cs_format->n+i, _cs_format->x[j]));
                REQUIRE(ret.second, "Internal Error");
            }
        }
        
        return reconstructedTensor;
    }
    
    
    CsUniquePtr matrix_matrix_product(const CsUniquePtr& _lhs, const CsUniquePtr& _rhs) {
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
