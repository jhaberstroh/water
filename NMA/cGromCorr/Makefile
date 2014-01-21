CPP = g++
CPPFLAGS:=-g -O2 -Wall -std=c++11

parse: main.cpp grom2atom.o
	$(CPP) $(CPPFLAGS) -o parse $^

clean: .
	rm *.o
