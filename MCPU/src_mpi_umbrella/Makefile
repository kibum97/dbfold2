# 1. Compiler and Flags
CC      = mpicc
CFLAGS  = -O3 -Wall -DDEBUG # -DDEBUG activates your #if DEBUG blocks
LIBS    = -lm              # Links the math library for sqrt, sin, cos

# 2. Identify Source Files
# This automatically finds all .c files in the current directory
SRCS    = $(wildcard *.c)

# 3. Derive Object Files
# This converts the list "vector.c init.c" into "vector.o init.o"
OBJS    = $(SRCS:.c=.o)

# 4. Name of your final executable
TARGET  = mcpu

# 5. Default Rule (The one that runs when you just type 'make')
all: $(TARGET)

# 6. Linking Rule
# This joins all the .o files into the final program
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# 7. Compilation Rule (Pattern Rule)
# This says: to create a .o file, look for a .c file of the same name.
# It also ensures that if a .h file changes, the .c file is recompiled.
%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

# 8. Clean Rule (Run 'make clean' to delete old builds)
clean:
	rm -f $(OBJS) $(TARGET)

# 9. Special rule to prevent confusion with files named 'clean' or 'all'
.PHONY: all clean