// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2017 Benjamin Huber and Sebastian Wolf. 
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
 * @brief Header file for the IHT algorithm and its variants.
 */

#pragma once

#include "../ttNetwork.h"

namespace xerus {
    
    TTTensor randomTTSVD(const Tensor& _x, const std::vector<size_t>& _ranks) {
        std::random_device rd;
        std::mt19937 rnd(rd());
        std::normal_distribution<double> dist(0, 1);
        
        const size_t d = _x.degree();
        TTTensor u(d);
        Tensor b = _x;
        
        for(long j = d; j--; j >= 2) {
            std::vector<size_t> gDims(_b.dimensions.cbegin(), _b.dimensions.begin()+(d-1));
            Tensor g(gDims, Tensor::Representation::Sparse);
            
            const auto& data = b.get_unsanitized_sparse_data();
            
            for(const auto& entry : data) {
                auto pos = Tensor::position_to_multiIndex(entry.first, b.dimensions);
                pos.pop_back();
                g[pos] = dist(rnd);
            }
            Tensor a;
            contract(a, g, false, b, false, j-1);
            
            Tensor R, Q;
            calculate_rq(R, Q, a, 1);
            
            u.set_component(j, Q);
            
            if(j == d) {
                contract(b, b, false, Q, true, 1);
            } else {
                contract(b, b, false, Q, true, 2);
            }
        }
        
        u.set_component(1, b);
    }
    
}


