CC = mpicxx
CFLAGS = -Iinclude -Wall -Wextra -pedantic
LDFLAGS = 

SRC = src/mpi_util.cpp src/align.cpp src/fold.cpp main.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = project

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

# Pattern rule to compile .cpp files to .o files
src/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean