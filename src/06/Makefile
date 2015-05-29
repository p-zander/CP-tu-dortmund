CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -Wfatal-errors -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')
CXXFLAGS+=-Wno-missing-prototypes -Wno-c++98-compat -Wno-unused-variable -Weverything

all: 1_plots 2_plots abgabe6.pdf

poisson: poisson.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) poisson.cpp 
	$(CXX) poisson.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy -o poisson

schroedinger: schroedinger.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) schroedinger.cpp 
	$(CXX) schroedinger.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy -o schroedinger

2_plots: poisson
	./poisson b
	./poisson c
	./poisson e 2

1_plots: schroedinger
	./schroedinger

abgabe6.pdf: abgabe6.tex
	latexmk -lualatex -f -silent abgabe6.tex
	$(MAKE) cleanlatex

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
	rm -f 1_Welle.gif

pack: all
	rm -f poisson
	rm -f poisson.o
	rm -f schroedinger
	rm -f schroedinger.o
	tar -czf abgabe6.tar.gz *

.PHONY: all clean cleanlatex pack 1_plots 2_plots
