CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11

all: hauptwert elektrostatik plots

hauptwert: hauptwert.cpp

elektrostatik: elektrostatik.cpp

plots: elektrostatik
	./elektrostatik
	python2 plots.py

clean:
	rm -f *.txt
	rm -f *.pdf
	rm -f hauptwert
	rm -f elektrostatik