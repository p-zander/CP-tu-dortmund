CXX=g++
CXXFLAGS=-O2 -Wall -Wextra -Werror -Wfatal-errors -pedantic -std=c++11 -fopenmp
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: 

rndwalk: rndwalk.cpp

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f rndwalk
	rm -f *.txt
	rm -f *.png
	rm -f *.tar.gz

pack: all
	tar -czf abgabe9.tar.gz *

.PHONY: all clean cleanlatex pack
