CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: Abgabe

harmoszi: harmoszi.cpp

kepler: kepler.cpp

plots: harmoszi kepler
	./harmoszi
	./kepler
	python2 plots.py

Abgabe: plots
	latexmk -pdf -pdflatex="lualatex" -use-make abgabe.tex
	make cleanlatex

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f harmoszi
	rm -f kepler
	rm -f *.txt
	rm -f *.tar.gz
	rm -f *pdf
pack: all
	rm -f harmoszi
	rm -f kepler
	tar -czf abgabe.tar.gz *
