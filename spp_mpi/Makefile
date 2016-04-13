# Makefile for the test programs

.SUFFIXES: .o .cpp

# Compiler and loader definitions
LD = ld
LDFLAGS =

CXX = mpicxx -std=c++0x -O4 -fopenmp -I ./include
CXXFLAGS = -g

MAKEFILE = Makefile

# list of all object and source files

OBJS =  LinkedList.o HashTable.o json_value.o json_reader.o json_writer.o Point.o Window.o EdgeCorrection.o SpatialPointPattern.o spp_mpi.o 
SRCS =	LinkedList.cpp HashTable.cpp json_value.cpp json_reader.cpp json_writer.cpp Point.cpp Window.cpp EdgeCorrection.cpp SpatialPointPattern.cpp spp_mpi.cpp 

all:		spp_mpi 

spp_mpi:	$(OBJS) 
		$(CXX) -o $@ $(OBJS) $(LDFLAGS)
		$(PURIFY) $(CXX) -o $@ $(OBJS) $(LDFLAGS)

.cpp.o:
		$(CXX) $(CXXFLAGS) -c $<

test:
	mpiexec -n 1 ./spp_mpi data/pnt/random_5000_sim50k.geojson data/us/us_sim50k.geojson 500000 2000 1 1 0 1 4
testcloud:
	mpiexec -hosts 10.128.96.246,10.128.96.250 -n 2 ./spp_mpi data/pnt/random_5000_sim50k.geojson data/us/us_sim50k.geojson 500000 2000 1 1 0 1 4
clean:
		rm -f core \#* *.bak *~ *.o spp_mpi vgcore.*
