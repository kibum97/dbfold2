CC = mpicc
CFLAGS = -Iinclude -Wall -Wextra -pedantic
LDFLAGS = 

SRC = src/mpi_util.c main.c
OBJ = $(SRC:.c=.o)
TARGET = project

all: $(TARGET)

$(TARGET): $(OBJ)
    $(CC) $(LDFLAGS) -o $@ $^

clean:
    rm -f $(OBJ) $(TARGET)

.PHONY: all clean