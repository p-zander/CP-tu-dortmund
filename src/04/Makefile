CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: Abgabe4.pdf

zeeman: zeeman.cpp

bifurk: bifurk.cpp

plots: zeeman bifurk
	./zeeman
	./bifurk
	python2 plots.py

Abgabe4.pdf: Abgabe4.tex plots
	latexmk -lualatex -use-make -f -silent Abgabe4.tex
	$(MAKE) cleanlatex

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean:	cleanlatex
	rm -f zeeman
	rm -f bifurk
	rm -f *.txt
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	rm -f zeeman
	rm -f bifurk
	rm -f 1*.pdf
	rm -f 2*.pdf
	tar -czf abgabe4.tar.gz *
