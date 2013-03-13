# This is a makefile for the DupHMM Viterbi algorithm


#############
# VARIABLES #
#############

APP = viterbi
SRC = viterbi.cpp
OBJ = viterbi.o
CC = g++
CFLAGS = -Wall -O3

###########
# TARGETS #
###########

default:
	make g++
	make clean_o

g++:
	$(CC) -c $(SRC) -o $(OBJ)
	$(CC) $(OBJ) -o $(APP) $(CFLAGS)

clean:
	rm -f *.o $(APP)

clean_o:
	rm -f *.o

test:
	echo "No test yet"
