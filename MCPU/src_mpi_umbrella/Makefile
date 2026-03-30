# 1. Compiler and Flags
CC      = mpicc
# -Iinclude tells the compiler to look in the include directory for .h files
CFLAGS  = -O3 -Wall -DDEBUG -Iinclude 
LIBS    = -lm

# 2. Directories
SRCDIR = src
INCDIR = include
OBJDIR = build
BINDIR = bin

# 3. Identify Source and Object Files
SRCS = $(wildcard $(SRCDIR)/*.c)
# This maps src/file.c to build/file.o
OBJS = $(SRCS:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

TARGET = $(BINDIR)/mcpu

# 4. Rules
all: $(OBJDIR) $(BINDIR) $(TARGET)

# Create directories if they don't exist
$(OBJDIR) $(BINDIR):
	mkdir -p $@

# Linking
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Compilation: Notice we no longer depend on %.h in the rule 
# because headers are now managed via -Iinclude and include guards
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean