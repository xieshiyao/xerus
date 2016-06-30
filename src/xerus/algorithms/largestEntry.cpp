
#include <xerus/algorithms/largestEntry.h>

namespace xerus {
	template<bool isOperator>
	size_t find_largest_entry(const TTNetwork<isOperator> &_T, const double _accuracy, const value_t _lowerBound) {
		_T.require_correct_format();
		
		// There is actual work to be done
		if(misc::sum(_T.ranks()) >= _T.degree()) {
			const double alpha = _accuracy;
			
			TTNetwork<isOperator> X = _T;
			X.round(size_t(1));
			double Xn = std::max(_T[find_largest_entry(X, 0.0, 0.0)], _lowerBound);
			double tau = (1-alpha)*alpha*Xn*Xn/(2.0*double(_T.degree()-1));
			
			X = _T;
			while(misc::sum(X.ranks()) >= _T.degree()) {
				X.entrywise_square();
				
				X.soft_threshold(tau, true);
				
				TTNetwork<isOperator> Y = X;
				Y.round(1);
				const size_t yMaxPos = find_largest_entry(Y, 0.0, 0.0);
				
				Xn = std::max(X[yMaxPos], (1-(1-alpha)*alpha/2.0)*Xn*Xn);
				
				const double fNorm = X.frob_norm();
				Xn /= fNorm;
				X /= fNorm;
				tau = (1-alpha)*alpha*Xn*Xn/(2.0*double(_T.degree()-1));
			}
			return find_largest_entry(X, 0.0, 0.0);
			
		// We are already rank one
		} else {
			const size_t numComponents = _T.degree()/(isOperator?2:1);
			size_t position = 0;
			size_t factor = misc::product(_T.dimensions);
			for(size_t c = 0; c < numComponents; ++c) {
				const size_t localSize = isOperator ? _T.dimensions[c]*_T.dimensions[numComponents+c] : _T.dimensions[c];
				factor /= localSize;
				
				size_t maxPos = 0;
				for(size_t i = 1; i < localSize; ++i) {
					if(std::abs(_T.get_component(c)[i]) > std::abs(_T.get_component(c)[maxPos])) {
						maxPos = i;
					}
				}
				position += maxPos*factor;
			}
			return position;
		}
	}
	
	template size_t find_largest_entry(const TTNetwork<true> &, double, value_t);
	template size_t find_largest_entry(const TTNetwork<false> &, double, value_t);
}
