CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3)

all: harmoszi

harmoszi: harmoszi.cpp

#plots: elektrostatik
#	./elektrostatik
#	python2 plots.py
#
#latex: plots
#	pdflatex abgabe.tex
#	make cleanlatex
#
#cleanlatex:
#	rm -f *.aux
#	rm -f *.fdb_latexmk
#	rm -f *.fls
#	rm -f *.log
#	rm -f *.synctex.gz

clean: 
	rm -f harmoszi
	rm -f *.tar.gz

pack: all
	rm -f harmoszi
	tar -czf abgabe.tar.gz *
