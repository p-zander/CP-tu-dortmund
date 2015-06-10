CXX=g++
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3 | sed 's/I/isystem/g')

all: plots plots2

lcg:lcg.cpp

lcg.txt: lcg
	./lcg

plots:lcg.txt plots.py 
	python2 plots.py

rand:rand.cpp

rand.txt:rand
	./rand

plots2:rand.txt plots2.py
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
	rm -f *.pdf
	rm -f *.png
	rm -f *.tar.gz

pack: all
	rm -f lcg
	rm -f rand
	tar -czf abgabe8.tar.gz *

.PHONY: all clean cleanlatex pack
