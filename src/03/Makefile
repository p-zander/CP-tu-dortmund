CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: abgabe.pdf

harmoszi: harmoszi.cpp

kepler: kepler.cpp

2a.txt 2b.txt: harmoszi 
	./harmoszi

a_tn_37_N_100_vy_0.88.txt a_tn_37_N_100_vy_1.3.txt a_tn_37_N_25_vy_1.3.txt bcd.txt e_alpha_0.9.txt e_alpha_1.1.txt: kepler
	./kepler

2a.pdf 2b_2.pdf 2b.pdf a1.pdf a2.pdf a3.pdf b.pdf c1.pdf c2.pdf d.pdf e09EL.pdf e09LR.pdf e11EL.pdf e11LR.pdf: \
2a.txt 2b.txt a_tn_37_N_100_vy_0.88.txt a_tn_37_N_100_vy_1.3.txt a_tn_37_N_25_vy_1.3.txt bcd.txt e_alpha_0.9.txt e_alpha_1.1.txt plots.py
	python2 plots.py

abgabe.pdf: 2a.pdf 2b_2.pdf 2b.pdf a1.pdf a2.pdf a3.pdf b.pdf c1.pdf c2.pdf d.pdf e09EL.pdf e09LR.pdf e11EL.pdf e11LR.pdf
	latexmk -lualatex -silent abgabe.tex
	$(MAKE) cleanlatex

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
	rm -f *.pdf

pack: all
	rm -f harmoszi
	rm -f kepler
	tar -czf abgabe.tar.gz *

.PHONY : all pack clean cleanlatex
.SILENT: pack clean cleanlatex 
