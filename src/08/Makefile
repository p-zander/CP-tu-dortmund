CXX=g++
CXXFLAGS=-O2 -Wall -Wextra -Werror -pedantic -std=c++11

all: lcg_corr_3.png lcg_corr.png lcg_hist.png dists.png

lcg: lcg.cpp

lcg.txt: lcg
	./lcg

lcg_corr_3.png: lcg.txt plots.py 
	python2 plots.py

lcg_corr.png: lcg.txt plots.py 
	python2 plots.py

lcg_hist.png: lcg.txt plots.py 
	python2 plots.py

rand: rand.cpp

rand.txt: rand
	./rand

dists.png: rand.txt plots2.py
	python2 plots2.py

cleanlatex:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.synctex.gz

clean: cleanlatex
	rm -f lcg
	rm -f rand
	rm -f *.txt
	rm -f *.png
	rm -f *.tar.gz

pack: all
	rm -f lcg
	rm -f rand
	tar -czf abgabe8.tar.gz *

.PHONY: all clean cleanlatex pack
