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
 * @brief Implementation of the retractions declared in retractions.h.
 */

#include <xerus.h>

namespace xerus {

	void HOSVDRetraction::operator()(TTTensor &_U, const TTTensor &_change) const {
		Index i;
		_U(i&0) = _U(i&0) + _change(i&0);
		if (roundByVector) {
			_U.round(rankVector);
		} else {
			_U.round(rank);
		}
	}
	
	void HOSVDRetraction::operator()(TTTensor &_U, const TTTangentVector &_change) const {
		Index i;
		_U = _change.added_to_base();
		if (roundByVector) {
			_U.round(rankVector);
		} else {
			_U.round(rank);
		}
	}
	
	void ALSRetractionII(TTTensor &_U, const TTTensor &_change) {
		static const ALSVariant roundingALS(1, 2, 0, ALSVariant::lapack_solver);
		TTTensor target = _U + _change;
		roundingALS(_U, target);
// 		_U.move_core(0);
	}
	
	void ALSRetractionI(TTTensor &_U, const TTTangentVector &_change) {
		static const ALSVariant roundingALS(1, 2, 0, ALSVariant::lapack_solver);
		TTTensor target = _change.added_to_base();
		roundingALS(_U, target);
// 		_U.move_core(0);
	}
	
	//TODO this should use low-level calls and then be part of FullTensor
// 	static FullTensor pseudoInverse(const FullTensor &_A) {
// 		Index i1,i2,i3,i4;
// 		FullTensor U,S,V;
// 		(U(i1,i2), S(i2,i3), V(i3,i4)) = SVD(_A(i1,i4));
// 		S.modify_diag_elements([](value_t &a){a = 1/a;});
// 		FullTensor res;
// 		res(i1,i4) = V(i2,i1) * S(i2,i3) * U(i4,i3);
// 		return res;
// 	}
	
	void TTTangentVector::set_base(const TTTensor &_newBase) {
		REQUIRE(_newBase.dimensions == baseL.dimensions, "");
		baseL = _newBase;
		baseL.move_core(0, true);
		baseR = _newBase;
		baseR.move_core(baseR.degree()-1, true);
	}
	
	TTTangentVector::TTTangentVector(const TTTensor& _base, const TTTensor& _direction) {
		REQUIRE(_base.is_valid_tt(), "");
		REQUIRE(_base.dimensions == _direction.dimensions, "");
		Index i1,i2,j1,j2,r,s;
		
		baseL = _base;
		baseL.move_core(0, true);
		baseR = _base;
		baseR.move_core(baseR.degree()-1, true);
		
		std::vector<FullTensor> rightStackUV;
		rightStackUV.push_back(Tensor::ones({1,1}));
		for (size_t i=baseL.degree()-1; i>0; --i) {
			FullTensor newRight;
			newRight(i1,i2) = baseL.get_component(i)(i1,r,j1) * _direction.get_component(i)(i2,r,j2) * rightStackUV.back()(j1,j2);
			rightStackUV.emplace_back(std::move(newRight));
		}
		FullTensor left(Tensor::ones({1,1}));
		
		// project onto the single components
		for (size_t i=0; i<baseR.degree(); ++i) {
			FullTensor V;
			V(i1,r,j1) =  left(i1,i2) * _direction.get_component(i)(i2,r,j2) * rightStackUV.back()(j1,j2);
			if (i>0) {
				// P^- part
				V(i1,r,j1) = V(i1,r,j1) - V(i1,s,j2) * baseL.get_component(i)(i2,s,j2) * baseL.get_component(i)(i2,r,j1);
			}
			components.emplace_back(std::move(V));
			if (i<baseR.degree()-1) {
				left(j1,j2) = left(i1,i2) * baseR.get_component(i)(i1,r,j1) * _direction.get_component(i)(i2,r,j2);
			}
			rightStackUV.pop_back();
		}
	}
	
	TTTangentVector& TTTangentVector::operator+=(const TTTangentVector& _rhs) {
		REQUIRE(components.size() == _rhs.components.size(), "");
		for (size_t i=0; i<components.size(); ++i) {
			components[i] += _rhs.components[i];
		}
		return *this;
	}
	
	TTTangentVector& TTTangentVector::operator-=(const TTTangentVector& _rhs) {
		REQUIRE(components.size() == _rhs.components.size(), "");
		for (size_t i=0; i<components.size(); ++i) {
			components[i] -= _rhs.components[i];
		}
		return *this;
	}
	
	TTTangentVector& TTTangentVector::operator*=(value_t _alpha) {
		for (size_t i=0; i<components.size(); ++i) {
			components[i] *= _alpha;
		}
		return *this;
	}
	
	TTTangentVector TTTangentVector::operator*(value_t _alpha) const {
		TTTangentVector result(*this);
		result *= _alpha;
		return result;
	}
	
	value_t TTTangentVector::scalar_product(const TTTangentVector& _other) const {
		Index j;
		value_t result = 0;
		for (size_t i=0; i<components.size(); ++i) {
			result += value_t(components[i](j^3)*_other.components[i](j^3));
		}
		return result;
	}
	
	TTTensor TTTangentVector::change_direction_incomplete() const {
		TTTensor result(baseL.degree());
		Index i1,i2,n,r1,r2;
		REQUIRE(components.size()>1, "ie: not meant to be called in this case");
		for (size_t i=0; i<components.size()-1; ++i) {
			if (i == 0) {
				FullTensor newComponent(4);
				newComponent(r1,n,i2,r2) = Tensor::dirac({2},{0})(i2) * baseR.get_component(i)(r1,n,r2)
											+Tensor::dirac({2},{1})(i2) * components[i](r1,n,r2);
				newComponent.reinterpret_dimensions({1, components[i].dimensions[1], components[i].dimensions[2]*2});
				result.set_component(i, newComponent);
			} else if (i < components.size()-1) {
				FullTensor newComponent(5);
				newComponent(i1,r1,n,i2,r2) = Tensor::dirac({2,2},{0,0})(i1,i2) * baseR.get_component(i)(r1,n,r2)
											+Tensor::dirac({2,2},{0,1})(i1,i2) * components[i](r1,n,r2)
											+Tensor::dirac({2,2},{1,1})(i1,i2) * baseL.get_component(i)(r1,n,r2);
				newComponent.reinterpret_dimensions({components[i].dimensions[0]*2, components[i].dimensions[1], components[i].dimensions[2]*2});
				result.set_component(i, newComponent);
			} else {
				LOG(fatal, "ie");
			}
			/*else {
				FullTensor newComponent(4);
				newComponent(i1,r1,n,r2) = Tensor::dirac({2},{0})(i1) * components[i](r1,n,r2)
											+Tensor::dirac({2},{1})(i1) * baseR.get_component(i)(r1,n,r2);
				newComponent.reinterpret_dimensions({components[i].dimensions[0]*2, components[i].dimensions[1], 1});
				result.set_component(i, newComponent);
			}*/ // will be done in the calling functions
		}
		return result;
	}
	
	TTTangentVector::operator TTTensor() const {
// 		if (components.size() == 1) {
// 			TTTensor result(1);
// 			result.set_component(0, components[0]);
// 			return result;
// 		}
// 		TTTensor result = change_direction_incomplete();
// 		FullTensor newComponent(4);
// 		Index i1,i2,n,r1,r2;
// 		newComponent(i1,r1,n,r2) = Tensor::dirac({2},{0})(i1) * components[components.size()-1](r1,n,r2)
// 									+Tensor::dirac({2},{1})(i1) * baseL.get_component(components.size()-1)(r1,n,r2);
// 		newComponent.reinterpret_dimensions({components[components.size()-1].dimensions[0]*2, components[components.size()-1].dimensions[1], 1});
// 		result.set_component(components.size()-1, newComponent);
// 		result.move_core(0);
// 		return result;
		TTTensor result(baseL);
		result.set_component(0, components[0]);
		for (size_t i=1; i<components.size(); ++i) {
			TTTensor tmp(baseL);
			tmp.move_core(i, true);
			tmp.set_component(i, components[i]);
			result += tmp;
		}
		result.move_core(0);
		return result;
	}
	
	TTTensor TTTangentVector::added_to_base() const {
// 		if (components.size() == 1) {
// 			TTTensor result(1);
// 			result.set_component(0, components[0]);
// 			result = result + baseL;
// 			return result;
// 		}
// 		TTTensor result = change_direction_incomplete();
// 		FullTensor newComponent(4);
// 		Index i1,i2,n,r1,r2;
// 		newComponent(i1,r1,n,r2) = Tensor::dirac({2},{0})(i1) * (components[components.size()-1](r1,n,r2) + baseR.get_component(components.size()-1)(r1,n,r2))
// 									+Tensor::dirac({2},{1})(i1) * baseL.get_component(components.size()-1)(r1,n,r2);
// 		newComponent.reinterpret_dimensions({components[components.size()-1].dimensions[0]*2, components[components.size()-1].dimensions[1], 1});
// 		result.set_component(components.size()-1, newComponent);
// 		result.move_core(0);
// 		return result;
		TTTensor result(baseL);
// 		result.set_component(0, components[0]);
		for (size_t i=0; i<components.size(); ++i) {
			TTTensor tmp(baseL);
			tmp.move_core(i, true);
			tmp.set_component(i, components[i]);
			result += tmp;
		}
		result.move_core(0);
		return result;
	}
	
	
	
	void SubmanifoldRetractionII(TTTensor &_U, const TTTensor &_change) {
		TTTangentVector W(_U, _change);
		SubmanifoldRetractionI(_U, W);
	}
	
	void SubmanifoldRetractionI(TTTensor &_U, const TTTangentVector &_change) {
		static const Index i1,j1,r;
		for (size_t i=0; i<_U.degree(); ++i) {
			std::unique_ptr<FullTensor> newComponent(new FullTensor);
			(*newComponent)(i1,r,j1) = _U.get_component(i)(i1,r,j1) + _change.components[i](i1,r,j1);
			_U.set_component(i, std::move(newComponent));
		}
		_U.move_core(0, true);
	}
	
	
	
	
	// TODO do this without creating the change_direction tensor?
	void ProjectiveVectorTransport(const TTTensor &_oldBase, const TTTensor &_newBase, TTTangentVector &_tangentVector) {
		REQUIRE(_oldBase.cannonicalized && _oldBase.corePosition == 0, "tangent vectors only implemented for core position 0 atm");
		REQUIRE(_newBase.cannonicalized && _newBase.corePosition == 0, "tangent vectors only implemented for core position 0 atm");
		
		_tangentVector.set_base(_oldBase);
		_tangentVector = TTTangentVector(_newBase, TTTensor(_tangentVector));
	}
}
