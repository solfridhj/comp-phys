CXX = g++
CXXF = -g -std=c++11 -Wall
CXXFLAGS =  -O3 -larmadillo -llapack -lblas

main: main.o lattice.o 
	$(CXX) -o main main.o lattice.o  $(CXXFLAGS)

main.o: main.cpp lattice.hpp 
	$(CXX) -c main.cpp $(CXXFLAGS)

lattice.o: lattice.hpp 


