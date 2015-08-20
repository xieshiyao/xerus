#include <xerus/misc/allocator.h>

using xma = xerus::misc::AllocatorStorage;

unsigned long xma::allocCount[xma::SMALLEST_NOT_CACHED_SIZE];
long xma::maxAlloc[xma::SMALLEST_NOT_CACHED_SIZE];
long xma::currAlloc[xma::SMALLEST_NOT_CACHED_SIZE];

__attribute__((constructor)) void initAllocatorStorage() {
	void *newPool = std::malloc(xma::POOL_SIZE);
	xma::pools.emplace_back(newPool, newPool); // TODO threads
	for (unsigned long i=0; i<xma::SMALLEST_NOT_CACHED_SIZE; ++i) {
		xma::allocCount[i]=0;
		xma::maxAlloc[i]=0;
		xma::currAlloc[i]=0;
	}
}


void* operator new(std::size_t n) {
    if (n>=xma::SMALLEST_NOT_CACHED_SIZE) {
		void *res = std::malloc(n+sizeof(size_t));
		*(static_cast<size_t*>(res)) = n;
		return static_cast<void*>(static_cast<size_t*>(res)+1);
	} else {
		void *res;
		if (xma::buckets[n].empty()) {
			if (xma::pools.back().second + n >= xma::pools.back().first + xma::POOL_SIZE) {
				void *newPool = std::malloc(xma::POOL_SIZE);
				xma::pools.emplace_back(newPool, newPool);
			}
			res = xma::pools.back().second;
			xma::pools.back().second += n+sizeof(size_t);
			*(static_cast<size_t*>(res)) = n;
			res = static_cast<void*>(static_cast<size_t*>(res)+1);
		} else {
			res = xma::buckets[n].back();
			xma::buckets[n].pop_back();
		}
		xma::allocCount[n] += 1;
		xma::currAlloc[n] += 1;
		if (xma::currAlloc[n] > xma::maxAlloc[n]) {
			xma::maxAlloc[n] = xma::currAlloc[n];
		}
		return res;
	}
}
void operator delete(void* ptr) noexcept {
    size_t n = *(static_cast<size_t*>(ptr)-1);
	if (n<xma::SMALLEST_NOT_CACHED_SIZE) {
		xma::currAlloc[n] -= 1;
		xma::buckets[n].push_back(ptr);
	} else {
		std::free(ptr);
	}
}
