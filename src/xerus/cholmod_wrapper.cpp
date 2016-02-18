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
* @brief Implementation of the suitesparse wrapper functions.
*/

#include <xerus/cholmod_wrapper.h>
#include <xerus/index.h>
#include <xerus/tensor.h>
#include <xerus/misc/performanceAnalysis.h>

#include <xerus/misc/check.h>

namespace xerus { namespace internal {
	
	CholmodCommon::RestrictedAccess::RestrictedAccess(cholmod_common* _c, std::mutex& _lock) 
		: c(_c), lock(_lock)
	{
		lock.lock();
	}
	
	CholmodCommon::RestrictedAccess::operator cholmod_common*() const {
		return c;
	}
	
	CholmodCommon::RestrictedAccess::~RestrictedAccess() {
		lock.unlock();
	}

	static void error_handler(int status, const char *file, int line, const char *message) {
		if (status < 0) {
			LOG(fatal, "CHOLMOD had a fatal error in " << file << ":" << line << " (status: " << status << ") msg: " << message);
		} else {
			LOG(cholmod_warning, "CHOLMOD warns in " << file << ":" << line << " (status: " << status << ") msg: " << message);
		}
	}
	
	CholmodCommon::CholmodCommon() : c(new cholmod_common()) {
		cholmod_start(c.get());
		c->itype = CHOLMOD_LONG;
		c->dtype = CHOLMOD_DOUBLE;
		c->error_handler = &error_handler;
// 		c->print = 1;
		REQUIRE(c->status == 0, "unable to initialize CHOLMOD");
	}
	
	CholmodCommon::~CholmodCommon() {
		cholmod_finish(c.get());
	}

	CholmodCommon::RestrictedAccess CholmodCommon::operator&() {
		return RestrictedAccess(c.get(), lock);
	}

	std::function<void(cholmod_sparse*)> CholmodCommon::get_deleter() {
		return [&](cholmod_sparse* _toDelete) {
			cholmod_free_sparse(&_toDelete, &(*this));
		};
	}
	
	static thread_local CholmodCommon cholmodObject;

	
	
	CholmodSparse::CholmodSparse(cholmod_sparse* _matrix) : matrix(_matrix, cholmodObject.get_deleter()) { }

	
	CholmodSparse::CholmodSparse(const size_t _m, const size_t _n, const size_t _N) 
		 : matrix(cholmod_allocate_sparse(_m, _n, _N, 1, 1, 0, CHOLMOD_REAL, &cholmodObject), cholmodObject.get_deleter())
	{
		REQUIRE(matrix, "cholmod_allocate_sparse did not allocate anything... status: " << cholmodObject.c->status << " call: " << _m << " " << _n << " " << _N << " alloc: " << cholmodObject.c->malloc_count);
	}

	CholmodSparse::CholmodSparse(const std::map<size_t, double>& _input, const size_t _m, const size_t _n, const bool _transpose) 
		: CholmodSparse(_m, _n, _input.size())
	{
		size_t entryPos = 0;
		LOG(test, uintptr_t(matrix->i));
		long* i = static_cast<long*>(matrix->i);
		long* p = static_cast<long*>(matrix->p);
		double* x = static_cast<double*>(matrix->x);
		i[0] = 0;
		
		// create compressed column storage of A^T aka compressed row storage of A
		
		long currRow = -1;
		
		for(const std::pair<size_t, value_t>& entry : _input) {
			x[entryPos] = entry.second;
			i[entryPos] = static_cast<long>(entry.first%_m);
			while(currRow < static_cast<long>(entry.first/_m)) {
				p[++currRow] = long(entryPos);
			}
			entryPos++;
		}
		
		REQUIRE(size_t(currRow) < _n && entryPos == _input.size(), "Internal Error " << currRow << ", " << _n << " | " << entryPos << ", " <<  _input.size());
		
		while(currRow < static_cast<long>(_n)) {
			p[++currRow] = long(entryPos);
		}
		
		if(!_transpose) {
			// we didn't want A^T, so transpose the data to get compressed column storage of A
			matrix = ptr_type(cholmod_transpose(matrix.get(), 1, &cholmodObject), cholmodObject.get_deleter());
		}
	}

	std::map<size_t, double> CholmodSparse::to_map(double _alpha) const {
		std::map<size_t, double> result;
		long* mi = reinterpret_cast<long*>(matrix->i);
		long* p = reinterpret_cast<long*>(matrix->p);
		double* x = reinterpret_cast<double*>(matrix->x);
		
		for(size_t i = 0; i < matrix->ncol; ++i) {
			for(long j = p[i]; j < p[i+1]; ++j) {
				IF_CHECK( auto ret = ) result.emplace(mi[size_t(j)]*matrix->ncol+i, _alpha*x[j]);
				REQUIRE(ret.second, "Internal Error");
			}
		}
		
		return result;
	}

	CholmodSparse CholmodSparse::operator*(const CholmodSparse& _rhs) const {
		return CholmodSparse(cholmod_ssmult(matrix.get(), _rhs.matrix.get(), 0, 1, 1, &cholmodObject));
	}

	
	void CholmodSparse::matrix_matrix_product( std::map<size_t, double>& _C,
								const size_t _leftDim,
								const size_t _rightDim,
								const double _alpha,
								const std::map<size_t, double>& _A,
								const bool _transposeA,
								const size_t _midDim,
								const std::map<size_t, double>& _B,
								const bool _transposeB ) 
	{
		const CholmodSparse lhsCs(_A, _leftDim, _midDim, _transposeA);
		const CholmodSparse rhsCs(_B, _midDim, _rightDim, _transposeB);
		const CholmodSparse resultCs = lhsCs * rhsCs;
		_C = resultCs.to_map(_alpha);
	}
}}
