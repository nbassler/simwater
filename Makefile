TARGET1   = h2ocalc
TARGET2   = time

CC       = gcc
CFLAGS   = -Wall -I.
LFLAGS   = -lgsl -lgslcblas -lm

SRCDIR   = src
OBJDIR   = obj

TARGETS  := $(TARGET1) $(TARGET2)
SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)


#$(TARGET): $(OBJECTS)
#	$(CC) $(OBJECTS) $(LFLAGS) -o $@

$(TARGET1): $(OBJDIR)/$(TARGET1).o
	$(CC) $(OBJDIR)/$(TARGET1).o $(LFLAGS) -o $@

$(TARGET2): $(OBJDIR)/$(TARGET2).o
	$(CC) $(OBJDIR)/$(TARGET2).o $(LFLAGS) -o $@

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -rf $(TARGETS) $(OBJECTS)


.PHONY: all
all: $(TARGETS)
