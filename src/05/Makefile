CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -Wfatal-errors -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')
CXXFLAGS+=-Wconversion -Wno-unused-parameter -Wno-unused-variable

all: plots

md: md.cpp
	$(CXX) -c $(shell python2.7-config --cflags) $(CXXFLAGS) md.cpp 
	$(CXX) md.o $(shell python2.7-config --ldflags) -lboost_python -lboost_numpy $(CXXFLAGS) -o md

plots: md plots.py
	./md

cleanlatex:
	rm -f *.au
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f md
	rm -f md.o
	rm -f *.pdf
	rm -f *.mp4
	rm -f *.gif
	rm -f *.png
	rm -f *.tar.gz

pack: all
	rm -f md
	rm -f md.o
	rm -f md*.pdf
	tar -czf abgabe5.tar.gz *

.PHONY: all clean cleanlatex pack plots
