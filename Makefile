#MAKEFILE FOR WIT PROJECT

.PHONY: clean all info

CC = gcc-8
CXX = g++-8

TARGETS := stiff tests main
SOURCES := $(TARGETS:=.cpp)
OBJS    := $(TARGETS:=.o)

OFLAGS := -O2 -O3
#-ffast-math
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread -DUSEOMP
#DEBUGFLAGS:= -g -Wno-stack-protector
CXXFLAGS := -std=c++14
LDFLAGS  :=
LIBS :=  -lstdc++ -I../Eigen
LIPOPT :=
OPENMPLIB :=

EXAMPLE_DEPS = Makefile

all: main

clean:
	rm -f $(OBJS) $(TARGETS)

info:
	@echo Compiler:	CC = $(CC) and CXX = $(CXX)
	@echo Compile command: COMPILE.cc = $(COMPILE.cc)
	@echo Link command: LINK.cc = $(LINK.cc)

stiff.o: stiff.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o stiff.o stiff.cpp
main.o: main.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o main.o main.cpp
main: main.o stiff.o
	@$(CXX) $(LDFLAGS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o main main.o stiff.o $(LIBS)
