
PREFIX = $(HOME)
CC = g++
AR = ar
CFLAGS= -g -O2 -std=c++11 -pthread

# path to the library files
LDLIBS  = -L/usr/local/lib -lntl -lgmp -lm  

# link in a library       
LDFLAGS = -I/usr/local/include      
LIBS    = -lntl -lm



HEADER =  Database.h  LRtest.h     HELR.h

SRC =    Database.cpp LRtest.cpp  HELR.cpp

OBJ= $(SRC:.cpp=.o)




lib: ar rc libheaan.a ../src/*.o

all: libHELR.a

clean:
	rm *.o libHELR.a  || echo nothing to clean
	
new:
	make clean
	make all


test:
	g++ -std=c++11 -O2 -I/usr/local/include -pthread Test.cpp libheaan.a libHELR.a -o foo -L/usr/local/lib -lntl -lgmp -lm
#./foo
#data/edin.txt



obj: $(OBJ)


%.o: %.cpp $(HEADER)
	 $(CC) $(CFLAGS) -c $(LDFLAGS)  $<  
  
libHELR.a: $(OBJ)
	$(AR) -q libHELR.a  $(OBJ)
