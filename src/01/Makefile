CXX=clang++
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
CXXFLAGS=$(shell pkg-config --cflags eigen3)

drehmoment: drehmoment.cpp

clean:
	rm -f drehmoment
	rm -f *.txt
