CXX=g++
CXXFLAGS=-O2 -Wall -Wextra -Werror -Wfatal-errors -pedantic -std=c++11 -fopenmp

all: abgabe 

abgabe: plots
	pdflatex abgabe.tex
	make cleanlatex

rndwalk: rndwalk.cpp

rnd_a_1.txt: rndwalk
	./rndwalk

mc_int: mc_int.cpp

mc_int.txt: mc_int
	./mc_int >> Integrale.txt

mc_int_hist.txt: mc_int
	./mc_int

mc_int_rndwalk.txt: mc_int
	./mc_int

plots: rnd_a_1.txt mc_int.txt mc_int_hist.txt mc_int_rndwalk.txt plots.py
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
	rm -f rndwalk
	rm -f mc_int
	rm -f *.png
	rm -f mc*.pdf
	rm -f mc*.txt rnd*.txt
	tar -czf abgabe10.tar.gz *

.PHONY: all plots clean cleanlatex pack
