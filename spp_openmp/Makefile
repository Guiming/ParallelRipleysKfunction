# Makefile for the test programs

.SUFFIXES: .o .cpp

# Compiler and loader definitions
LD = ld
LDFLAGS =

CXX = g++ -std=c++0x -O4 -fopenmp -I ./include -m64
CXXFLAGS = -g

MAKEFILE = Makefile

# list of all object and source files

OBJS =  LinkedList.o HashTable.o json_value.o json_reader.o json_writer.o Point.o Window.o EdgeCorrection.o SpatialPointPattern.o spp.o 
SRCS =	LinkedList.cpp HashTable.cpp json_value.cpp json_reader.cpp json_writer.cpp Point.cpp Window.cpp EdgeCorrection.cpp SpatialPointPattern.cpp spp.cpp 

all:		spp 

spp:	$(OBJS) 
		$(CXX) -o $@ $(OBJS) $(LDFLAGS)
		$(PURIFY) $(CXX) -o $@ $(OBJS) $(LDFLAGS)
test:
	./spp data/pnt/random_5000_sim50k.geojson data/us/us_sim50k.geojson output/result_sim50k.csv 700000 2000 1 1 1 100000 4 1 1 1 1 16 > output/log.log

.cpp.o:
		$(CXX) $(CXXFLAGS) -c $<

clean:
		rm -f core \#* *.bak *~ *.o spp vgcore.*

