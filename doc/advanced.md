# Notes for Advanced Users

## Multi-Threading

Please note, that `xerus` is only thread-safe up to 1024 threads at the time of this writing. This number is due to the internal handling of indices. To ensure that indices can 
be compared and identified uniquely, they store a unique id of which the first 10 bits are reserved to denote the number of the current thread. With more than 1024 threads (when these
10 bits overflow) it can thus lead to collisions and indices that were shared between threads and were meant to be seperate suddenly evaluate to be equal to all algorithms.


## R-Value References

Any movable object of the type `xerus::IndexedTensorMovable<T>` is assumed to be manipulatable be the receiving functions without influence to any outside storage.

Storing an r-value reference to such an object (`xerus::IndexedTensorMovable<T> &&a = ...`) will thus break the library. Badly. Just don't do it...

This is because these object typically arise as the returned objects of operators on indexed tensors. Eg. `A(i,j)*B(j,k)` will create a temporary of exactly this type. That it is
of the type `xerus::IndexedTensorMovable<T>` and an r-value reference denotes to further functions (like further multiplications) that this is a temporary and may thus be recycled to 
reduce the total number of temporary objects created...

