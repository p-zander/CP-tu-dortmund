CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: Abgabe4.pdf

zeeman: zeeman.cpp

bifurk: bifurk.cpp

1_ab.txt 1_b_16_0_4.txt 1_b_16_2_4.txt 2_a.txt 2_c.txt 1_b_16_0_2.txt 1_b_16_2_2.txt 1_c.txt 2_b.txt: zeeman.cpp bifurk.cpp
	./zeeman
	./bifurk

1_ab_F.pdf 1_b.pdf 2_a_0.pdf 2_a_2.pdf 2_a_4.pdf 2_b.pdf 1_ab_V.pdf 1_c.pdf 2_a_1.pdf 2_a_3.pdf 2_a_5.pdf 2_c.pdf: \
zeeman bifurk plots.py 1_ab.txt 1_b_16_0_4.txt 1_b_16_2_4.txt 2_a.txt 2_c.txt 1_b_16_0_2.txt 1_b_16_2_2.txt 1_c.txt 2_b.txt
	python2 plots.py
	ls *.pdf | grep -v 'Abgabe4.pdf' | parallel -N 1 pdfcrop {1} {1}

Abgabe4.pdf: Abgabe4.tex plots.py zeeman.cpp bifurk.cpp
	latexmk -lualatex -use-make -f -silent Abgabe4.tex
	$(MAKE) cleanlatex

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
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

.PHONY: all clean cleanlatex pack
