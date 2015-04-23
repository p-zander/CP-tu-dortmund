CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11

all: hauptwert elektrostatik latex

hauptwert: hauptwert.cpp

elektrostatik: elektrostatik.cpp

plots: elektrostatik
	./elektrostatik
	python2 plots.py

latex: plots
	pdflatex abgabe.tex
	make cleanlatex

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f hauptwert
	rm -f elektrostatik
	rm -f *.txt
	rm -f *.pdf
	rm -f *.tar.gz

pack: clean
	tar -czf abgabe.tar.gz *
