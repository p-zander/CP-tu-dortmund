# CXX=g++
# CXXFLAGS=-O2 -Wall -std=c++11
# CXXFLAGS+=$(shell python2.7-config --cflags --ldflags) -lboost_python 

default:
	g++ -c $(shell python2.7-config --cflags) --std=c++11 utils_test.cpp 
	g++ utils_test.o $(shell python2.7-config --ldflags) -lboost_python --std=c++11 -o utils_test

clean:
	rm -f utils_test
	rm -f utils_test.o
	rm -f *.pdf
