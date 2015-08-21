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


#include<xerus.h>

#include "../../include/xerus/misc/test.h"
using namespace xerus;

#ifdef PERFORMANCE_ANALYSIS
    UNIT_TEST(X_PerformanceAnalysis_X, Analysis,
		using xma = xerus::misc::AllocatorStorage;
        std::cout << misc::performanceAnalysis::get_analysis();
        LOG(Indices, "A total of " << Index().valueId << " indices were used (in this thread).");
		LOG(allocator, "");
		size_t totalStorage=0;
		for (size_t i=0; i<xma::NUM_BUCKETS; ++i) {
			if (xma::allocCount[i] == 0 && xma::currAlloc[i] == 0) continue;
			totalStorage += i * xma::BUCKET_SIZE * (size_t)xma::maxAlloc[i];
			LOG(allocator, i * xma::BUCKET_SIZE << " \tx\t " << xma::allocCount[i] << "\tmax: " << xma::maxAlloc[i] << '\t' << xma::currAlloc[i] << '\t' << totalStorage);
		}
		LOG(storageNeeded, totalStorage << " storage used: " << misc::astore.pools.size()*xma::POOL_SIZE);
		LOG(index, sizeof(Index));
		LOG(node, sizeof(TensorNetwork::TensorNode));
		LOG(link, sizeof(TensorNetwork::Link));
		LOG(tn, sizeof(TensorNetwork));
		LOG(fulltensor, sizeof(FullTensor));
		LOG(SparseTensor, sizeof(SparseTensor));
		LOG(tt, sizeof(TTTensor) << " | " << sizeof(TTOperator));
		LOG(measurement, sizeof(SinglePointMeasurment));
    )
#endif
