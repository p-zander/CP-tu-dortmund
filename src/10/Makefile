CXX=g++
CXXFLAGS=-O2 -Wall -Wextra -Werror -Wfatal-errors -pedantic -std=c++11 -fopenmp

all: plots 

rndwalk: rndwalk.cpp

rnd.txt: rndwalk
	./rndwalk

mc_int: mc_int.cpp

mc_int.txt: mc_int
	./mc_int

mc_int_hist.txt: mc_int
	./mc_int

mc_int_rndwalk.txt: mc_int
	./mc_int

plots: rnd.txt mc_int.txt mc_int_hist.txt mc_int_rndwalk.txt plots.py
	python2 plots.py

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f rndwalk
	rm -f mc_int
	rm -f *.txt
	rm -f *.png
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	tar -czf abgabe10.tar.gz *

.PHONY: all plots clean cleanlatex pack
