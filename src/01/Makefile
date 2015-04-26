CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS+=$(shell pkg-config --cflags eigen3)

drehmoment: drehmoment.cpp

data: drehmoment
	./drehmoment

integration: integration.cpp

plots: data
	./plots.py

clean:
	rm -f drehmoment
	rm -f *.txt
	rm -f integration
	rm -f *.pdf