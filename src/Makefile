# For Linux or any machines with gcc compiler
CC = gcc 
CFLAGS = -std=c99 -Wall -pedantic
BIN = ..\bin

all: Sec3_1

clean:
	del *.o *.exe $(BIN)\*.exe lib\*.o 

OBJ = lib\tiff.o lib\allocate.o lib\randlib.o lib\qGGMRF.o lib\solve.o

Sec3_1: Sec3_1.o $(OBJ)
	$(CC) $(CFLAGS) -o Sec3_1 Sec3_1.o $(OBJ) -lm
