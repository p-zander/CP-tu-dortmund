CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
#CXXFLAGS+=$(shell pkg-config --cflags eigen3)

all: hauptwert elektrostatik plots

hauptwert: hauptwert.cpp

elektrostatik: elektrostatik.cpp

plots: elektrostatik
	./elektrostatik 11 80 > Aufgabe2a.txt
	./elektrostatik  0  9 > Aufgabe2b.txt
	python2 plot.py

clean:
	rm -f *.txt
	rm -f *.pdf
	rm -f hauptwert
	rm -f elektrostatik