#include <xerus/misc/allocator.h>

unsigned long xerus::misc::allocatorStorage::allocCount[1000];
long xerus::misc::allocatorStorage::maxAlloc[1000];
long xerus::misc::allocatorStorage::currAlloc[1000];

__attribute__((constructor)) void bla() {
	for (unsigned long i=0; i<1000; ++i) {
		xerus::misc::allocatorStorage::allocCount[i]=0;
		xerus::misc::allocatorStorage::maxAlloc[i]=0;
		xerus::misc::allocatorStorage::currAlloc[i]=0;
	}
}

