CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: md

# md: md.cpp
md: md.cpp plots.py
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) md.cpp 
	$(CXX) md.o $(shell python2.7-config --ldflags) -lboost_python $(CXXFLAGS) -o md

cleanlatex:
	rm -f *.au
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f md
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	rm -f md
	rm -f md*.pdf
	tar -czf abgabe5.tar.gz *

.PHONY: all clean cleanlatex pack
