CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -Wfatal-errors -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')
CXXFLAGS+=-Wno-missing-prototypes -Wno-unused-variable -Weverything -Wno-c++98-compat 

all: plots 

bifurk: bifurk.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) bifurk.cpp 
	$(CXX) bifurk.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy $(CXXFLAGS) -o bifurk

plots: plots.py 
	./bifurk

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f bifurk
	rm -f bifurk.o
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	rm -f bifurk
	rm -f bifurk.o
	tar -czf abgabe7.tar.gz *

.PHONY: all clean cleanlatex pack plots
