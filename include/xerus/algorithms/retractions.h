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
 * @brief Header file for several types of retractions.
 */

#pragma once

#include <vector>
#include "../fullTensor.h"

namespace xerus {

	///@brief class to compactly represent tangent vectors of the manifold of constant TT-rank
	class TTTangentVector {
	public:
		std::vector<FullTensor> components;
		///@brief creates a tangent vector by projecting @a _direction onto the tangent plane located at @a _base
		TTTangentVector(const TTTensor &_base, const TTTensor &_direction);
		TTTangentVector &operator+=(const TTTangentVector &_rhs);
		TTTangentVector &operator*=(value_t _alpha);
		TTTangentVector operator*(value_t _alpha);
		value_t scalar_product(const TTTensor &_base, const TTTangentVector &_other);
	};
	
	/// retraction that performs a HOSVD to project back onto the Manifold
	struct HOSVDRetraction {
		bool roundByVector;
		size_t rank;
		std::vector<size_t> rankVector;
		void operator()(TTTensor &_U, const TTTensor &_change) const;
		HOSVDRetraction(size_t _rank) : roundByVector(false), rank(_rank) {}
		HOSVDRetraction(const std::vector<size_t> &_rank) : roundByVector(true), rank(~0ul), rankVector(_rank) {}
	};
	
	/// retraction that performs an ALS half-sweep to project back onto the Manifold. Automatically retains the ranks of @a _U
	void ALSRetraction(TTTensor &_U, const TTTensor &_change);
	
	/// retraction that performs componentwise addition of @f$ U_i @f$ and @f$ W_i @f$ where @f$ W_i @f$ is the i-th component of the riemannian tangential vector representation
	void SubmanifoldRetraction(TTTensor &_U, const TTTensor &_change);
	
}

