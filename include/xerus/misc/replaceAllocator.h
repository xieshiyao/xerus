#ifndef _GLIBCXX_CXX_ALLOCATOR_H
#define _GLIBCXX_CXX_ALLOCATOR_H 1

#include "allocator.h"

namespace std
{
  /// sets the default allocator for all stl containers (needs to be set before any stl library is included!)
  template<typename _Tp>
    using __allocator_base = xerus::misc::DebugAllocator<_Tp>;
}

#endif
