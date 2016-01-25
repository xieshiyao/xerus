# Notes for Advanced Users

## Multi-Threading

Please note, that `xerus` is only thread-safe up to 1024 threads at the time of this writing. This number is due to the internal handling of indices. To ensure that indices can 
be compared and identified uniquely, they store a unique id of which the first 10 bits are reserved to denote the number of the current thread. With more than 1024 threads (when these
10 bits overflow) it can thus lead to collisions and indices that were shared between threads and were meant to be seperate suddenly evaluate to be equal to all algorithms.


