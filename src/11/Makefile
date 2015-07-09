CXX=g++
CXXFLAGS=-O2 -Wall -Wextra -Werror -Wfatal-errors -pedantic -std=c++11 -fopenmp

all: plots

mc_1: mc_1.cpp

mc_2: mc_2.cpp

plots: mc_1 mc_2
	./mc_1
	./mc_2
	python2 plots.py

clean:
	rm -f mc_1
	rm -f mc_2
	rm -f *.txt
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	tar -czf abgabe11.tar.gz *

.PHONY: all plots clean pack
