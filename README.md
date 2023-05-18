## Linear algebra C++ static library prototype

The idea of this code was to provide an easy way to operate with matrices in C++, with the aim of comparing the numerical stability of two variants of the QR matrix factorization algorithm, more precisely:

1. the one based on Householder reflectors,

2. the classic version based on Gram-Schmidt.

There are several modifications that must still be undertaken to make this an actual functional C++ library, to name a few:

- Fixing the "rigidity" of the structs that serve as containers for matrices. This can probably be done by the use of dynamimc vectors instead of arrays of doubles.

- Adding a function to calculate the determinant of matrices.

- Arranging the code and files in the appropiate way for a static c++ library, which implies moving the functions definitions from the header file into a .cpp file.

On another note, there is still the task to change the random matrix generator function to have the appropiate quantity of decimals in each entry of the matrices so as to explicitly be able to note the differences in the numerical stability of the two variants of the QR algorithm.
