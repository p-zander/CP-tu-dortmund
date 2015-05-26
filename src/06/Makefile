CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -Wfatal-errors -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')
CXXFLAGS+=-Wno-missing-prototypes -Wno-c++98-compat -Wno-unused-variable -Weverything

all: plots plots2

poisson: poisson.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) poisson.cpp 
	$(CXX) poisson.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy $(CXXFLAGS) -o poisson

schroedinger: schroedinger.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) schroedinger.cpp 
	$(CXX) schroedinger.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy $(CXXFLAGS) -o schroedinger

plots: poisson
	./poisson

plots2: schroedinger
	./schroedinger

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f poisson
	rm -f poisson.o
	rm -f schroedinger
	rm -f schroedinger.o
	rm -f *.pdf
	rm -f *.tar.gz
	rm -f Welle.gif

pack: all
	rm -f poisson
	rm -f poisson.o
	rm -f schroedinger
	rm -f schroedinger.o
	tar -czf abgabe5.tar.gz *

.PHONY: all clean cleanlatex pack plots
