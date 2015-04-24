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

#include <vector>
#include <limits>
#include "exceptions.h"
#include "missingFunctions.h"
#include "namedLogger.h"
#include "test.h"

START_MISC_NAMESPACE

/* 


#include "simpleNumerics.h"

int main() {
	LOG(pi, "pi starting");
	LOG(pi, std::setw(14) << "n" << std::setw(14) << "sum" << std::setw(14) << "shanks" << std::setw(14) << "err" << std::setw(14) << "richard" << std::setw(14) << "err");
	double sum = 0;
	ShanksTransformation<double> s;
	RichardsonExtrapolation<double> r;
	s.push_back(sum);
	r.push_back(sum);
	for (size_t n=0; n<25; ++n) {
		sum += (n%2==0?4:-4) / double(2*n+1);
		s.push_back(sum); r.push_back(sum);
		LOG(pi, std::setw(14) << n << std::setw(14) << sum << std::setw(14) << s.best_estimate() << std::setw(14) << s.error_approximate() << std::setw(14) << r.best_estimate() << std::setw(14) << r.error_approximate());
	}
	
	LOG(pi, "pi^2/6 starting");
	LOG(pi, std::setw(14) << "n" << std::setw(14) << "sum" << std::setw(14) << "shanks" << std::setw(14) << "err" << std::setw(14) << "richard" << std::setw(14) << "err");
	sum = 0;
	s.reset(); s.push_back(sum);
	r.reset(); r.push_back(sum);
	for (size_t n=1; n<25; ++n) {
		sum += 1/double(n*n);
		s.push_back(sum); r.push_back(sum);
// 		LOG(asdkj, r.values);
		LOG(pi, std::setw(14) << n << std::setw(14) << sum << std::setw(14) << s.best_estimate() << std::setw(14) << s.error_approximate() << std::setw(14) << r.best_estimate() << std::setw(14) << r.error_approximate());
	}
	
	return 0;
}
 */


/// performs a Romberg Integration (richardson extrapolation of regular riemannian sum) + adaptive refinement
double integrate(const std::function<double(double)> &_f, double _a, double _b, double _eps=std::numeric_limits<double>::epsilon(), 
				 uint _minIter=4, uint _maxIter=6, uint _branchFactor=7, 
				 uint _maxRecursion=10, bool _relativeError=true);

double integrate_segmented(const std::function<double(double)> &_f, double _a, double _b, double _segmentation, 
						   double _eps=1e-8, uint _minIter=4, uint _maxIter=6, uint _branchFactor=8,
						   uint _maxRecursion=10);




/// class to represent a polynomial by its vector of coefficients
struct Polynomial {
	std::vector<double> coefficients;
	Polynomial(const std::vector<double> _coeff) : coefficients(_coeff) {}
	Polynomial() {}
	
	size_t terms() const {
		return coefficients.size();
	}
	Polynomial &operator-=(const Polynomial &_rhs) {
		if (terms() < _rhs.terms()) {
			coefficients.resize(_rhs.terms());
		}
		for (size_t i=0; i<std::min(terms(), _rhs.terms()); ++i) {
			coefficients[i] -= _rhs.coefficients[i];
		}
		return *this;
	}
	Polynomial operator*(const Polynomial &_rhs) const {
		Polynomial result;
		result.coefficients.resize(_rhs.terms()+terms()-1);
		for (size_t i=0; i<terms(); ++i) {
			for (size_t j=0; j<_rhs.terms(); ++j) {
				result.coefficients[i+j] += coefficients[i]*_rhs.coefficients[j];
			}
		}
		return result;
	}
	Polynomial &operator/=(double _rhs) {
		for (size_t i=0; i<terms(); ++i) {
			coefficients[i] /= _rhs;
		}
		return *this;
	}
	Polynomial &operator*=(double _rhs) {
		for (size_t i=0; i<terms(); ++i) {
			coefficients[i] *= _rhs;
		}
		return *this;
	}
	Polynomial operator*(double _rhs) const {
		Polynomial result;
		result.coefficients.resize(terms());
		for (size_t i=0; i<terms(); ++i) {
			result.coefficients[i] = coefficients[i]*_rhs;
		}
		return result;
	}
	
	double operator()(double x) const {
		double result = 0;
		for (size_t i=terms(); i>0; --i) {
			result *= x;
			result += coefficients[i-1];
		}
		return result;
	}
	double scalar_product(const Polynomial &_rhs, const std::function<double (double)> &_weight, double _minX, double _maxX) const {
		return integrate([&](double x){
			return (*this)(x) * _rhs(x) * _weight(x);
		}, _minX, _maxX, 1e-10);
	}
	double norm(const std::function<double (double)> &_weight, double _minX, double _maxX) const {
		return std::sqrt(scalar_product(*this, _weight, _minX, _maxX));
	}
	
	/// orthogonalizes this polynomial with respect to the provided (@note orthogonal!) basis
	Polynomial &orthogonolize(const std::vector<Polynomial> &_orthoBase, const std::function<double (double)> &_weight, double _minX, double _maxX) {
		for (size_t i=0; i<_orthoBase.size(); ++i) {
			(*this)-=_orthoBase[i]*scalar_product(_orthoBase[i], _weight, _minX, _maxX);
			REQUIRE(std::abs(scalar_product(_orthoBase[i], _weight, _minX, _maxX)) < 1e-12, i << " " << std::abs(scalar_product(_orthoBase[i], _weight, _minX, _maxX)));
		}
		(*this)/=norm(_weight, _minX, _maxX);
		return *this;
	}
	
	/// returns @a _N pairwise orthogonal polynomials w.r.t. a scalar product defined by the @a _weight
	static std::vector<Polynomial> build_orthogonal_base(uint _N, const std::function<double (double)> &_weight, double _minX, double _maxX) {
		std::vector<Polynomial> base;
		while (base.size() < _N) {
			Polynomial next;
			next.coefficients.resize(base.size()+1);
			next.coefficients.back() = 1; // next = x^(base.size)
			next.orthogonolize(base, _weight, _minX, _maxX);
			base.push_back(next);
			LOG(debug, "next basis polynomial " << next.coefficients);
		}
		return base;
	}
};


/// classes that can extract an estimate of the limit of a sequence
template<class ft_type>
class LimitExtractor {
public:
	virtual void push_back(ft_type _val) = 0;
	virtual ft_type best_estimate() const = 0;
	virtual ft_type error_approximate() const = 0;
	virtual void reset() = 0;
};

/** limit extraction using the shanks transformation aka Aitken process
 * derivation by assuming the sequence to go as x_n = x_inf + alpha * q^n for large n
 */
template<class ft_type>
class ShanksTransformation : public LimitExtractor<ft_type> {
public:
	static constexpr ft_type epsilon = std::numeric_limits<ft_type>::epsilon();
public:
	std::vector<ft_type> values;
	
	static ft_type shanks(ft_type x1, ft_type x2, ft_type x3) {
		ft_type denominator = x1 - 2*x2 + x3;
		if (std::abs(denominator) < epsilon * std::max(x1, std::max(x2, x3))) {
			return x2;
		}
		return (x1*x3 - x2*x2) / denominator;
	}
	
public:
	void push_back(ft_type _val) override {
		values.push_back(_val);
		for (size_t i=values.size()-1; i>=2; i-=2) {
			values[i-2] = shanks(values[i-2], values[i-1], values[i]);
		}
	}
	
	ft_type best_estimate() const override {
		if (values.size() == 0) {
			XERUS_THROW(generic_error() << "tried to extract limit of empty sequence"); 
		} else {
			return values[(values.size()-1) % 2];
		}
	}
	
	ft_type error_approximate() const override {
		size_t i = (values.size()-1) % 2;
		if (i+1 >= values.size()) {
			return 1;
		} else {
			return std::abs(values[i] - values[i+1]);
		}
	}
	
	void reset() {
		values.clear();
	}
};


//TODO the following is crap. implement Levin-t, Levin-u Levin-v
/** limit extraction using the richardson extrapolation
 * derivation by assuming that x_inf - x_n = alpha * n^(-1)
 */
template<class ft_type>
class RichardsonExtrapolation : public LimitExtractor<ft_type> {
public:
	static constexpr ft_type epsilon = std::numeric_limits<ft_type>::epsilon();
public:
	std::vector<ft_type> values;
	
	static ft_type richard(size_t n, ft_type x1, ft_type x2) {
		return double(n+1)*x2 - double(n)*x1;
	}
	
public:
	void push_back(ft_type _val) override {
		values.push_back(_val);
		for (size_t i=values.size()-1; i>=1; i-=1) {
			values[i-1] = richard(i-1, values[i-1], values[i]);
		}
	}
	
	ft_type best_estimate() const override {
		if (values.size() == 0) {
			XERUS_THROW(generic_error() << "tried to extract limit of empty sequence"); 
		} else {
			return values.front();
		}
	}
	
	ft_type error_approximate() const override {
		if (values.size() < 2) {
			return 1;
		} else {
			return std::abs(values[0] - values[1]);
		}
	}
	
	void reset() {
		values.clear();
	}
};

END_MISC_NAMESPACE
