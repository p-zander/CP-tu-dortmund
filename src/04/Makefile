CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: plots

zeeman: zeeman.cpp

bifurk: bifurk.cpp

plots: zeeman bifurk
	./zeeman
	./bifurk
	python2 plots.py

# latex: plots
# 	pdflatex abgabe.tex
# 	make cleanlatex

# cleanlatex:
# 	rm -f *.aux
# 	rm -f *.fdb_latexmk
# 	rm -f *.fls
# 	rm -f *.log
# 	rm -f *.synctex.gz

clean:
	rm -f zeeman
	rm -f bifurk
	rm -f *.txt
	rm -f *.pdf
	rm -f *.tar.gz

pack: all
	rm -f zeeman
	rm -f bifurk
	tar -czf abgabe.tar.gz *
