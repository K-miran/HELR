
PREFIX = $(HOME)
CC = g++
AR = ar
CFLAGS= -g -O2 -std=c++11 -pthread
#-std=c++11

# path to the library files
LDLIBS  = -L/usr/local/lib -lntl -lgmp -lm  

# link in a library       
LDFLAGS = -I/usr/local/include      
LIBS    = -lntl -lm



HEADER = CZZ.h  CZZX.h EvaluatorUtils.h KsiPows.h NumUtils.h Ring2Utils.h StringUtils.h   TimeUtils.h   Cipher.h Message.h Params.h PubKey.h SchemeAux.h Scheme.h  SecKey.h Database.h LRtest.h LogisticReg.h
#LinearReg.h  
SRC =   CZZ.cpp  CZZX.cpp  EvaluatorUtils.cpp KsiPows.cpp NumUtils.cpp Ring2Utils.cpp StringUtils.cpp  TimeUtils.cpp    Cipher.cpp Message.cpp Params.cpp PubKey.cpp SchemeAux.cpp  Scheme.cpp  SecKey.cpp Database.cpp LRtest.cpp LogisticReg.cpp
#LinearReg.cpp 
OBJ= $(SRC:.cpp=.o)



all: libHELR.a

clean:
	rm *.o libHELR.a  || echo nothing to clean

new:
	make clean
	make all
	g++ -std=c++11 -O2 -I/usr/local/include -pthread main.cpp libHELR.a  -o foo -L/usr/local/lib -lntl -lgmp -lm 


runedin3:
	./foo data/edin.txt 3

runedin7:
	./foo data/edin.txt 7


obj: $(OBJ)


%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHELR.a: $(OBJ)
	$(AR) -q libHELR.a  $(OBJ)
