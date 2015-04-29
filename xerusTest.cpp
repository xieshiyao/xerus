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

// File containing the unit tests for the xerus library

#include "xerus.h"
#include <test.h>    
   
//------------------------ Unit Tests --------------------------------
using namespace xerus;
    #include "unitTests/sparseTensor.hxx"
    #include "unitTests/sparseTensor_add_sub_cpy.hxx"
    #include "unitTests/sparseTensor_arithmetic_cpy.hxx"
    #include "unitTests/sparseTensor_assignment_cpy.hxx"
    #include "unitTests/sparseTensor_contraction_cpy.hxx"
    #include "unitTests/sparseTensor_product_cpy.hxx"
    #include "unitTests/sparseTensor_trace_cpy.hxx"
    #include "unitTests/sparseTensor_assignment.hxx"
//         #include "unitTests/sparseTensor_contraction.hxx"
    #include "unitTests/tensorNetwork.hxx"
    #include "unitTests/ttCreation.hxx"
    #include "unitTests/ttRounding.hxx"
    #include "unitTests/ttArithmetic.hxx"
    #include "unitTests/xBlasAnalysis.hxx"
    #include "unitTests/als.hxx"
    #include "unitTests/misc.hxx"
//      #include "unitTests/strassen.hxx"


    
int main() {
    LOG(fatal, "It looks like test was not active, although it had to be.");
    return 0;
}
