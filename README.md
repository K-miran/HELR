# HELR

HELR is a software project for performing a logistic regression training on encrypted data (Secure Logistic Regression based on Homomorphic Encryption: Design and Evaluation (https://medinform.jmir.org/2018/2/e19/))

## Setting up HEMat library 

### Dependencies

Our library requires a c++ compiler and the following libraries:

* `GMP`(GNU Multi-Precision library), which is available at https://gmplib.org,

* `NTL`, which is available at http://www.shoup.net/ntl/,  (with pThread)

* `HEAAN`,  which is an implementation of the paper "Homomorphic Encryption for Arithmetic of Approximate Numbers" (https://eprint.iacr.org/2016/421.pdf). We refered to the underlying HE library in the "src" folder. You can build the libarary “libheaan.a" by typing "$make all" in the "/src" directory.

```sh
make clean
make all
```

### Installing HELR library

HELR is easy to configure and build in Linux and macOS. You can then install our library by the following work with command line tools:

```sh
make new
```

## Running a test source code

For example, you run a test program main.cpp with a filename and a degree of approximating polynomial of the sigmoid function.
For example, you can write:

```sh
g++ -std=c++11 -O2 -I/usr/local/include -pthread main.cpp libHELR.a  -o foo -L/usr/local/lib -lntl -lgmp -lm
./foo data/edint.txt 3 
```

In particular, our program supports the evaluation of the gradient descent algorithm based on polynomial of degree 3 or degree 7 with several optimizations.
