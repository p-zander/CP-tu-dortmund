CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -Wfatal-errors -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')
CXXFLAGS+=-Weverything -Wno-missing-prototypes

all: poisson

poisson: poisson.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) poisson.cpp 
	$(CXX) poisson.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy $(CXXFLAGS) -o poisson

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f poisson
	rm -f poisson.o
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	rm -f poisson
	rm -f poisson.o
	tar -czf abgabe5.tar.gz *

.PHONY: all clean cleanlatex pack
