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
		_U = _change.added_to_base(_U);
		if (roundByVector) {
			_U.round(rankVector);
		} else {
			_U.round(rank);
		}
	}
	
	void ALSRetractionII(TTTensor &_U, const TTTensor &_change) {
		static const ALSVariant roundingALS(1, 2, 0, ALSVariant::lapack_solver);
		Index i;
		TTTensor target;
		target(i&0) = _U(i&0) + _change(i&0);
		roundingALS(_U, target);
// 		_U.move_core(0);
	}
	
	void ALSRetractionI(TTTensor &_U, const TTTangentVector &_change) {
		static const ALSVariant roundingALS(1, 2, 0, ALSVariant::lapack_solver);
		Index i;
		TTTensor target = _change.added_to_base(_U);
		roundingALS(_U, target);
// 		_U.move_core(0);
	}
	
	//TODO this should use low-level calls and then be part of FullTensor
	static FullTensor pseudoInverse(const FullTensor &_A) {
		Index i1,i2,i3,i4;
		FullTensor U,S,V;
		(U(i1,i2), S(i2,i3), V(i3,i4)) = SVD(_A(i1,i4));
		S.modify_diag_elements([](value_t &a){a = 1/a;});
		FullTensor res;
		res(i1,i4) = V(i2,i1) * S(i2,i3) * U(i4,i3);
		return res;
	}
	
	
	TTTangentVector::TTTangentVector(const TTTensor& _base, const TTTensor& _direction) {
		REQUIRE(_base.cannonicalized && _base.corePosition == 0, "projection onto tangent plane is only implemented for core position 0 at the moment");
		REQUIRE(_base.dimensions == _direction.dimensions, "");
		Index i1,i2,j1,j2,r,s;
		std::vector<FullTensor> leftStackUV;
		std::vector<FullTensor> leftStackUU;
		FullTensor tmp({1,1}, [](){return 1.0;});
		leftStackUV.push_back(tmp);
		leftStackUU.push_back(tmp);
		for (size_t i=0; i<_base.degree()-1; ++i) {
			FullTensor newLeft;
			newLeft(j1,j2) = leftStackUV.back()(i1,i2) * _base.get_component(i)(i1,r,j1) * _direction.get_component(i)(i2,r,j2);
			leftStackUV.emplace_back(std::move(newLeft));
			newLeft(j1,j2) = leftStackUU.back()(i1,i2) * _base.get_component(i)(i1,r,j1) * _base.get_component(i)(i2,r,j2);
			leftStackUU.emplace_back(std::move(newLeft));
		}
		FullTensor right(tmp);
		FullTensor UTV;
		std::vector<FullTensor> tmpComponents;
		for (size_t i=_base.degree(); i>0; --i) {
			const size_t currIdx = i-1;
			std::unique_ptr<FullTensor> newComponent(new FullTensor);
			const Tensor &UComp = _base.get_component(currIdx);
			FullTensor V;
			FullTensor uuInv = pseudoInverse(leftStackUU.back());
			V(i1,r,j1) =  uuInv(i1,s) * leftStackUV.back()(s,i2) * _direction.get_component(currIdx)(i2,r,j2) * right(j1,j2);
			if (i!=_base.degree()) {
				V(i1,r,j1) = V(i1,r,j1) + UComp(i1,r,s)*UTV(s,j1);
			}
			if (currIdx!=0) {
				UTV(i1,i2) = V(i1,r,j1) * UComp(i2,r,j1);
				V(i1,r,j1) = V(i1,r,j1) - UTV(i1,s) * UComp(s,r,j1);
			}
			tmpComponents.emplace_back(std::move(V));
			if (currIdx != 0) {
				right(j1,j2) = UComp(j1,r,i1) * _direction.get_component(currIdx)(j2,r,i2) * right(i1,i2);
			}
			leftStackUV.pop_back();
			leftStackUU.pop_back();
		}
		while (!tmpComponents.empty()) {
			components.emplace_back(std::move(tmpComponents.back()));
			tmpComponents.pop_back();
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
	
	value_t TTTangentVector::scalar_product(const TTTensor& _base, const TTTangentVector& _other) const {
		REQUIRE(_base.degree() == components.size(), "");
		REQUIRE(components.size() == _other.components.size(), "");
		REQUIRE(_base.cannonicalized && _base.corePosition == 0, "tangent vectors only implemented for core position 0 atm");
		Index i1,i2,r,j1,j2;
		value_t result = 0;
		FullTensor left({1,1}, [](){return 1.0;});
		for (size_t i=0; i<components.size(); ++i) {
			result += value_t(left(i1,i2)*components[i](i1,r,j1)*_other.components[i](i2,r,j1));
			if (i < components.size()-1) {
				left(j1,j2) = left(i1,i2) * _base.get_component(i)(i1,r,j1) * _base.get_component(i)(i2,r,j2);
			}
		}
		return result;
	}
	
	TTTensor TTTangentVector::change_direction_not_orthogonalized(const TTTensor &_base) const {
		REQUIRE(_base.cannonicalized && _base.corePosition == 0, "tangent vectors only implemented for core position 0 atm");
		REQUIRE(_base.degree() == components.size(), "");
		TTTensor result(_base.degree());
		Index i1,i2,n,r1,r2;
		if (components.size() == 1) {
			result.set_component(0, components[0]);
			return result;
		}
		for (size_t i=0; i<components.size(); ++i) {
			const Tensor &baseComponent = _base.get_component(i);
			REQUIRE(baseComponent.dimensions == components[i].dimensions, "illegal base tensor for this tangent vector (ranks or external dimension do not coincide)");
			if (i == 0) {
				FullTensor newComponent(4);
				newComponent(r1,n,i2,r2) = Tensor::dirac({2},{0})(i2) * baseComponent(r1,n,r2)
											+Tensor::dirac({2},{1})(i2) * components[i](r1,n,r2);
				newComponent.reinterpret_dimensions({1, components[i].dimensions[1], components[i].dimensions[2]*2});
				result.set_component(i, newComponent);
			} else if (i < components.size()-1) {
				FullTensor newComponent(5);
				newComponent(i1,r1,n,i2,r2) = Tensor::dirac({2,2},{0,0})(i1,i2) * baseComponent(r1,n,r2)
											+Tensor::dirac({2,2},{0,1})(i1,i2) * components[i](r1,n,r2)
											+Tensor::dirac({2,2},{1,1})(i1,i2) * baseComponent(r1,n,r2);
				newComponent.reinterpret_dimensions({components[i].dimensions[0]*2, components[i].dimensions[1], components[i].dimensions[2]*2});
				result.set_component(i, newComponent);
			} else {
				FullTensor newComponent(4);
				newComponent(i1,r1,n,r2) = Tensor::dirac({2},{0})(i1) * components[i](r1,n,r2)
											+Tensor::dirac({2},{1})(i1) * baseComponent(r1,n,r2);
				newComponent.reinterpret_dimensions({components[i].dimensions[0]*2, components[i].dimensions[1], 1});
				result.set_component(i, newComponent);
			}
		}
		return result;
	}
	
	TTTensor TTTangentVector::change_direction(const TTTensor &_base) const {
		TTTensor result(_base);// = change_direction_not_orthogonalized(_base);
		result.set_component(0, components[0]);
		for (size_t i=1; i<components.size(); ++i) {
			TTTensor tmp(_base);
			tmp.set_component(i, components[i]);
			result += tmp;
		}
		result.move_core(0);
		return result;
	}
	
	TTTensor TTTangentVector::added_to_base(const TTTensor &_base) const {
		TTTensor result(components.size());
		Index i1,i2,n,r1,r2;
		if (components.size() == 1) {
			FullTensor newComponent(3);
			newComponent(r1,n,r2) = _base.get_component(0)(r1,n,r2) + components[0](r1,n,r2);
			result.set_component(0, newComponent);
		} else {
			result = change_direction_not_orthogonalized(_base);
			size_t pos = components.size()-1;
			const Tensor &baseComponent = _base.get_component(pos);
			FullTensor newComponent(4);
			newComponent(i1,r1,n,r2) = Tensor::dirac({2},{0})(i1) * (components[pos](r1,n,r2) + baseComponent(r1,n,r2))
										+Tensor::dirac({2},{1})(i1) * baseComponent(r1,n,r2);
			newComponent.reinterpret_dimensions({components[pos].dimensions[0]*2, components[pos].dimensions[1], 1});
			result.set_component(pos, newComponent);
		}
		result.move_core(0, true);
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
		
		_tangentVector = TTTangentVector(_newBase, _tangentVector.change_direction(_oldBase));
	}
}
