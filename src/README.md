This is our underlying homomorphic encryption scheme, 
which is “Homomorphic Encryption for Arithmetic of Approximate Numbers” (HEAAN).

We referred to its library in gitHub (https://github.com/kimandrik/HEAAN) and 
here is the version at 2017/09.

$ gcc -c -I/usr/local/include -I. -std=c++11 *.cpp

$ ar rc libheaan.a ../src/*.o

This will build the HEAAN library “libheaan.a”. 
