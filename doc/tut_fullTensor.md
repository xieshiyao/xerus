# Getting started

In this part of the tutorial the basic functionalities of xerus for FullTensors are introduced. 

In xerus the class FullTensor is used to represent all dense tensors, i.e. tensors that are not sparse and not decomposed into a tensor network. 
Lets start of by creating a simple 3x3 matrix, i.e. a tensor of order two and dimensions `{3,3}`.
\code{.cpp}
int main() {
    xerus::FullTensor A({3,3});
\endcode
By default xerus zero-initializes every entry of a created tensor. 
The elements of a tensor object can be accesed via the `[]` operator. For example 
\code{.cpp}
    A[{0,0}] = 1.0;
    A[{1,1}] = 1.0;
    A[{2,2}] = 1.0;
\endcode
makes `A` the 3x3 identity matrix. Xerus also allows to create random constructed tensors via `C++11` random facilities. For example
\code{.cpp}
    std::mt19937_64 rnd;
    std::normal_distribution<value_t> dist (0.0, 10.0);
    xerus::FullTensor B = xerus::FullTensor::construct_random({10,10,5}, rnd, dist);
\endcode
creates a `10x10x5` tensor with independend gaussian distributed random entries.

Xerus supports all basic arithmetic for FullTensors, e.g.



A different approach to create a new Tensor is to use a function based constructor, e.g.
\code{.cpp}
    xerus::FullTensor B({3,3}, [](const std::vector<size_t>& _position) {
        if (_position[0] == _position[1]) { 
            return 1.0; 
        } else { 
            return 0.0; 
        }
    }));
\endcode
which const


to set the entry `2,2` one can use

alternativelty on can use the single integer overload of the `[]` operator to directly operate on the data array. Xerus uses row-major order in all representations. In our example this means
\code{.cpp}
    B[0] = 1.0; // Sets the entry (0,0)
    B[1] = 0.0; // Sets the entry (0,1)
    B[2] = 0.0; // Sets the entry (0,2)
    B[3] = 0.0; // Sets the entry (1,0)
    B[4] = 1.0; // Sets the entry (1,1)
    B[5] = 0.0; // Sets the entry (1,2)
    B[6] = 0.0; // Sets the entry (2,0)
    B[7] = 0.0; // Sets the entry (2,1)
\endcode
Now `B` is the 3x3 identity matrix. 



constructing a random fulltensor
some contractions with different index order, `i^d` and `i&1`
