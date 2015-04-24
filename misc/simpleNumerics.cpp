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

#include "simpleNumerics.h"
#include "missingFunctions.h"

START_MISC_NAMESPACE

double integrate(const std::function<double(double)> &_f, double _a, double _b, double _eps, 
				 uint _minIter, uint _maxIter, uint _branchFactor, 
				 uint _maxRecursion, bool _relativeError) {
	
	REQUIRE(_minIter > 1, "");
	REQUIRE(_branchFactor > 1, "");
	double min = std::min(_a,_b);
	double max = std::max(_a,_b);
	if (_relativeError) _eps = std::max(_eps, std::numeric_limits<double>::epsilon());
	std::vector<double> iterants;
	iterants.push_back((max - min)*(_f(min) + _f(max))/2);
	double h = (max-min);
	double sum;
	double error = 1;
	double maxVal = std::abs(iterants[0]);
	for (uint iter=0; iter<_maxIter; ++iter) {
		sum = 0;
		for (double x = min+h/2; x<max; x += h) {
			double f = _f(x);
			sum += f;
			maxVal = std::max(maxVal, std::abs(f));
		}
		
		h /= 2;
		sum *= h;
		sum += iterants.back()/2;
		iterants.push_back(sum);
		double oldIt0 = iterants[0];
		for (size_t k=0; k<iterants.size()-1; ++k) {
			size_t i=iterants.size()-1-k;
			iterants[i-1] = iterants[i] + (iterants[i]-iterants[i-1])/(pow(2.0, 2*(k+1))-1);
		}
		if (_relativeError) {
			error = std::abs((iterants[0]-oldIt0)/oldIt0);
			if (std::isnan(error)) error = std::abs(iterants[0]-oldIt0);
		} else {
			error = std::abs(iterants[0]-oldIt0);
		}
		if (iter >= _minIter && error <= _eps) {
			return (_a>_b?-1:1)*iterants[0];
		}
	}
	if (_maxRecursion == 0) {
		return (_a>_b?-1:1)*iterants[0];
	}
	// divide and conquer
	// this dynamicaly divides those parts that did not converge easily
	LOG(debug, (_relativeError?"relative":"absolute") << " error " << error << " is larger than eps = " << _eps);
	sum = 0;
	h = (max-min) / _branchFactor;
	double newEps = (_relativeError
						? std::max(std::abs(iterants[0]), maxVal)*_eps 
						: std::max(
							_eps, // /std::sqrt(_branchFactor), 
							std::sqrt(_branchFactor)*std::numeric_limits<double>::epsilon()*std::max(std::abs(iterants[0]), maxVal)
						  )
					);
	for (uint i=0; i<_branchFactor; ++i) {
		sum += integrate(_f, min+i*h, min+(i+1)*h, newEps, _minIter, _maxIter, _branchFactor, _maxRecursion-1, false);
	}
	return (_a>_b?-1:1)*sum;
}


double integrate_segmented(const std::function<double(double)> &_f, double _a, double _b, double _segmentation, 
						   double _eps, uint _minIter, uint _maxIter, uint _branchFactor,
						   uint _maxRecursion) {
	double min = std::min(_a,_b);
	double max = std::max(_a,_b);
	double res = 0;
	for (double x=min; x<max; x+=_segmentation) {
		res += integrate(_f, x, std::min(x+_segmentation, max), _eps, _minIter, _maxIter, _branchFactor, _maxRecursion);
	}
	return (_a>_b?-1:1)*res;
}


END_MISC_NAMESPACE
