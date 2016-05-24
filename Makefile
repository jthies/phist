include ./make.defines

CXX = mpiCC
CF = mpif90
CFFLAGS += -fPIC -openmp -O3 -warn all 
CPPFLAGS +=  -fPIC -Wall -openmp -openmp-link=static -pthread -O3 -pedantic $(INCLUDES) $(DEFINES)
LIBS +=
INCLUDES +=
TARGET = libcpaft.so
LFLAGS =

OBJS = cp.o aft.o


$(TARGET): $(OBJS)
	$(CXX) -shared -o $(TARGET) $(OBJS) 

cp.o: cp.cpp cp.h cp_array.h cp_ghostdensemat.h
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LIBS) -c cp.cpp

aft.o: aft.cpp aft.h aft_macros.h
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LIBS) -c aft.cpp

clean:
	rm -f *.o *.so
