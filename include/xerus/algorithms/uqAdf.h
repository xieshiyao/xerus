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
 * @brief Header file for the ADF algorithm and its variants.
 */

#pragma once

#include "../index.h"
#include "../ttNetwork.h"
#include "../performanceData.h"
#include "../measurments.h"

namespace xerus {
    
	
	class UQMeasurementSet {
	public:
		std::vector<std::vector<double>> randomVectors;
		std::vector<Tensor> solutions;
		
		std::vector<std::vector<double>> initialRandomVectors;
		std::vector<Tensor> initialSolutions;
		
		UQMeasurementSet() = default;
		UQMeasurementSet(const UQMeasurementSet&  _other) = default;
		UQMeasurementSet(      UQMeasurementSet&& _other) = default;
		
		void add(const std::vector<double>& _rndvec, const Tensor& _solution);
		
		void add_initial(const std::vector<double>& _rndvec, const Tensor& _solution);
	};
	
	
	
	void uq_adf(TTTensor& _x, const std::vector<std::vector<double>>& _randomVariables, const std::vector<Tensor>& _solutions);
	
	TTTensor uq_adf(const UQMeasurementSet& _measurments, const TTTensor& _guess);
	
	Tensor uq_avg(const TTTensor& _x, const size_t _N, const size_t _numSpecial);
	
	std::pair<std::vector<std::vector<double>>, std::vector<Tensor>> uq_mc(const TTTensor& _x, const size_t _N, const size_t _numSpecial);
}

