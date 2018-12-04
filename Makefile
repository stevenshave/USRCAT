CC = clang++

ldflags =  
libs = -lm -lopenbabel
debug = 
ccopt = -O3 -I /usr/include/openbabel-2.0/ -Wall -std=c++1y
all: usrcat

usrcat: Usrcat.cpp
	$(CC) -o $@ $^ $(libs) $(ccopt)

.PHONY: clean
clean:
	/bin/rm -f usrcat

