#MAKEFILE FOR WIT PROJECT
#THE WINNING PEAR

.PHONY: clean all main tests info play

CC = gcc-8
CXX = g++-8

SRCDIR := src/
BINDIR := bin/

TARGETS := stiff boundary function tests main play
SOURCES := $(TARGETS:=.cpp)
OBJS    := $(TARGETS:=.o)

OFLAGS := -O2 -O3
#-ffast-math
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread -DUSEOMP
#DEBUGFLAGS:= -g -Wno-stack-protector
CXXFLAGS := -std=c++14
LDFLAGS  :=
LIBS :=  -lstdc++ -I Eigen
LIPOPT :=
OPENMPLIB :=

EXAMPLE_DEPS = Makefile

all: tests $(BINDIR)main

play: $(BINDIR)play
	./$(BINDIR)play

tests: $(BINDIR)tests
	./$(BINDIR)tests

clean:
	rm -rf $(BINDIR)*

info:
	@echo Compiler:	CC = $(CC) and CXX = $(CXX)
	@echo Compile command: COMPILE.cc = $(COMPILE.cc)
	@echo Link command: LINK.cc = $(LINK.cc)

$(BINDIR)stiff.o: $(SRCDIR)stiff.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)stiff.o $(SRCDIR)stiff.cpp
$(BINDIR)boundary.o: $(SRCDIR)boundary.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)boundary.o $(SRCDIR)boundary.cpp
$(BINDIR)function.o: $(SRCDIR)function.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)function.o $(SRCDIR)function.cpp

$(BINDIR)tests.o: $(SRCDIR)tests.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)tests.o $(SRCDIR)tests.cpp
$(BINDIR)tests: $(BINDIR)tests.o $(BINDIR)stiff.o $(BINDIR)boundary.o $(BINDIR)function.o
	@$(CXX) $(LDFLAGS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)tests $(BINDIR)tests.o $(BINDIR)stiff.o $(BINDIR)boundary.o $(BINDIR)function.o $(LIBS)

$(BINDIR)main.o: $(SRCDIR)main.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)main.o $(SRCDIR)main.cpp
$(BINDIR)main: $(BINDIR)main.o $(BINDIR)stiff.o $(BINDIR)boundary.o $(BINDIR)function.o
	@$(CXX) $(LDFLAGS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)main $(BINDIR)main.o $(BINDIR)stiff.o $(BINDIR)boundary.o $(BINDIR)function.o $(LIBS)

$(BINDIR)play.o: $(SRCDIR)play.cpp $(EXAMPLE_DEPS)
	@$(CXX) -c $(CXXFLAGS) $(OFLAGS) $(LIBS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)play.o $(SRCDIR)play.cpp
$(BINDIR)play: $(BINDIR)play.o
	@$(CXX) $(LDFLAGS) $(PARALLELFLAGS) $(DEBUGFLAGS) -o $(BINDIR)play $(BINDIR)play.o $(LIBS)
