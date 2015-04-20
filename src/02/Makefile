CXX=clang++-3.5
CXXFLAGS=-O2 -Wall -Werror -pedantic -std=c++11
#CXXFLAGS+=$(shell pkg-config --cflags eigen3)

all: hauptwert elektrostatik

hauptwert: hauptwert.cpp

elektrostatik: elektrostatik.cpp

clean:
	rm -f hauptwert
	rm -f elektrostatik