# HELR
Secure Logistic Regression based on Homomorphic Encryption

Our code requires the NTL library which is available at http://www.shoup.net/ntl/, and a c++ compiler. 
Our underlying homomorphic encryption scheme is “Homomorphic Encryption for Arithmetic of Approximate Numbers” 
and we referred to its library in the "src" folder.

  $make clean
  
  $make all

This will build our library “libHELR.a”. Then the following should work:

  $g++ -std=c++11 -O2 -I/usr/local/include -pthread main.cpp libHELR.a  -o foo -L/usr/local/lib -lntl -lgmp -lm

for compiling the code and making your program main.cpp
The program will run with a filename and a degree of approximating polynomial of the sigmoid function.
For example, you can write:

  $./foo data/edint.txt 3 

In particular, our program supports the evaluation of the gradient descent algorithm based on polynomial of degree 3 or degree 7 with several optimizations.
