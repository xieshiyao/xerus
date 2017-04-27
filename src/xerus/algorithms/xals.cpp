// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
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
 * @brief Implementation of the ADF variants. 
 */

#include <xerus/algorithms/xals.h>
 
#include <xerus/indexedTensorMoveable.h>
#include <xerus/misc/basicArraySupport.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace xerus {
    
     class InternalSolver {
        const size_t d;
        
         std::vector<Tensor> leftStack;
         std::vector<Tensor> rightStack;
         
         
        TTTensor& x;
        const TTOperator& A;
        const TTTensor& b;
         
     public:
         InternalSolver(TTTensor& _x, const TTOperator& _A, const TTTensor& _b) : d(_x.degree()), x(_x), A(_A), b(_b) {
             
         }
         
         void solve() {
             
         }
     };

        void xals(TTTensor& _x, const TTOperator& _A, const TTTensor& _b)  {
            InternalSolver solver(_x, _A, _b);
            return solver.solve();
        }
}

