CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -Wfatal-errors -pedantic -std=c++11 
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g') $(shell python2.7-config --cflags) -Wno-missing-prototypes
CXXFLAGS+=-Wno-unused-variable -Wno-c++98-compat
LDFLAGS+=$(shell python2.7-config --ldflags) -lboost_python -lboost_numpy

all: plots 

bifurk: bifurk.cpp
	$(CXX) -c $(CXXFLAGS) bifurk.cpp 
	$(CXX) bifurk.o $(LDFLAGS) -o bifurk

feigenbaum: feigenbaum.cpp
	$(CXX) -c $(CXXFLAGS) feigenbaum.cpp 
	$(CXX) feigenbaum.o $(LDFLAGS) -o feigenbaum

plots: plots.py bifurk feigenbaum
	./bifurk
	./feigenbaum > 2_nullstellen.txt

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f bifurk
	rm -f bifurk.o
	rm -f feigenbaum
	rm -f feigenbaum.o
	rm -f *.pdf
	rm -f *.txt
	rm -f *.tar.gz

pack: all
	rm -f bifurk
	rm -f bifurk.o
	rm -f feigenbaum
	rm -f feigenbaum.o
	tar -czf abgabe7.tar.gz *

.PHONY: all clean cleanlatex pack plots
