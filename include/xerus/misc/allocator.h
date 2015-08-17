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
#include <ext/new_allocator.h>

namespace xerus { namespace misc {

struct allocatorStorage {
	static unsigned long allocCount[1000];
	static long maxAlloc[1000];
	static long currAlloc[1000];
};

template <class Tp>
struct DebugAllocator : __gnu_cxx::new_allocator<Tp> {
	typedef Tp value_type;
	
	DebugAllocator() {};
	template <class T> DebugAllocator(const DebugAllocator<T>& other) {};
	Tp* allocate(unsigned long n) {
		n = n * sizeof(Tp);
		if (n<1000) {
			allocatorStorage::allocCount[n] += 1;
			allocatorStorage::currAlloc[n] += 1;
			if (allocatorStorage::currAlloc[n] > allocatorStorage::maxAlloc[n]) {
				allocatorStorage::maxAlloc[n] = allocatorStorage::currAlloc[n];
			}
		}
		// TODO check whether n*sizeof overflows
		return static_cast<Tp*>(::operator new(n));
	}
	void deallocate(Tp* p, unsigned long n) {
		if (p == nullptr) return;
		n = n * sizeof(Tp);
		if (n<1000) {
			allocatorStorage::currAlloc[n] -= 1;
// 			if (allocatorStorage::currAlloc[n] < 0) {
// 				allocatorStorage::currAlloc[n] = 0;
// 			}
		}
		::operator delete(p);
	}
};
template <class T, class U>
bool operator==(const DebugAllocator<T>&, const DebugAllocator<U>&) {
	return true;
}
template <class T, class U>
bool operator!=(const DebugAllocator<T>&, const DebugAllocator<U>&) {
	return false;
}

}}



