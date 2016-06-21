include ./make.defines

CXX = mpiCC
CF = mpif90
CFFLAGS += -fPIC -openmp -O3 -warn all 
CPPFLAGS +=  -fPIC -Wall -openmp -pthread -O3 -DSCR -pedantic 
LIBS +=
INCLUDES +=
TARGET = libcpaft.so
LFLAGS =

OBJS = cp.o aft.o cp_options.o 


$(TARGET): $(OBJS)
	$(CXX)  -shared -o $(TARGET) $(OBJS) $(INCLUDES) $(LIB_PATH) $(LIBS)  

cp.o: cp.cpp cp.h cp_array.h cp_ghostdensemat.h cpPOD.h
	$(CXX) $(CPPFLAGS)  $(INCLUDES) $(LIB_PATH) $(LIBS) -c cp.cpp

cp_options.o: cp_options.cpp cp_options.h 
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LIB_PATH) $(LIBS) -c cp_options.cpp

aft.o: aft.cpp aft.h aft_macros.h
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(LIB_PATH) $(LIBS) -c aft.cpp

clean:
	rm -f *.o *.so
