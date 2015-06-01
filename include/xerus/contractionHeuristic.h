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

#include "tensorNetwork.h"

namespace xerus {

    //TODO individuelle schranken
    
    /// Class representing heuristics used to find good contraction orders of TensorNetworks.
    class ContractionHeuristic {
    public:
        static std::vector<ContractionHeuristic> *list;
        
        struct AddToVector {
            template<typename... args>
            AddToVector(args... _args) {
                if (!list) list = new std::vector<ContractionHeuristic>();
                list->emplace_back(_args...);
            }
        };
        
		/// Name of the heuristic
        std::string name;
        
        double score;
        
        std::vector<std::pair<size_t,size_t>> contractions;
        
        std::function<void(double &, std::vector<std::pair<size_t,size_t>> &, TensorNetwork &)> scoreFct;
        
        ContractionHeuristic(std::string _name, std::function<void(double &, std::vector<std::pair<size_t,size_t>> &, TensorNetwork &)> _scoreFct);
        
        double rescore(TensorNetwork _tn); // NOTE take as value to get a deep copy instead of reference!
    };
}

