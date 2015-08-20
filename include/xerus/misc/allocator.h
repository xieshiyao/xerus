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
 * @brief Header file for xerus::misc::allocators
 */


#pragma once
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>

namespace xerus { namespace misc {

template <class T>
struct Mallocator {
	typedef T value_type;
	Mallocator();
	template <class S> Mallocator(const Mallocator<S>& other) {}
	T* allocate(std::size_t n) {
		if (~0ul / sizeof(T) < n) {
			throw new std::bad_alloc();
		}
		return std::malloc(n*sizeof(T));
	}
	void deallocate(T* p, std::size_t ) {
		std::free(p);
	}
};
template <class T, class U>
bool operator==(const Mallocator<T>&, const Mallocator<U>&) {
	return true;
}
template <class T, class U>
bool operator!=(const Mallocator<T>&, const Mallocator<U>&) {
	return false;
}

struct AllocatorStorage {
	static constexpr size_t SMALLEST_NOT_CACHED_SIZE = 2048;
	static constexpr size_t POOL_SIZE = 9*1024*1024;
	
	static thread_local std::array<std::vector<void*, Mallocator<void*>>, SMALLEST_NOT_CACHED_SIZE> buckets;
	static thread_local std::vector<std::pair<char*, char*>, Mallocator<void*>> pools;
	
	static unsigned long allocCount[SMALLEST_NOT_CACHED_SIZE];
	static long maxAlloc[SMALLEST_NOT_CACHED_SIZE];
	static long currAlloc[SMALLEST_NOT_CACHED_SIZE];
};

}}

void* operator new(std::size_t n);
void operator delete(void* ptr) noexcept;





