# Debug
DEBUG ?= 1
# Target executable
TARGET := main

# C files
CFILES := 
#    gpujpeg_reader.c \
#    main.c

#if you want to add the cpp files , please put here.
# CPP files
CPPFILES := \
    main.cpp \
    ClContext.cpp \
    ColorConv.cpp \

# Compilers
CC := g++ -fPIC

# Debug
ifeq ($(DEBUG),1)
    COMMONFLAGS += -g -D_DEBUG -O0
else
    COMMONFLAGS += -O2
endif

# Common flags
COMMONFLAGS += -I. -I./include
# C flags
CFLAGS += $(COMMONFLAGS) -std=c99
# CPP flags
CPPFLAGS += `pkg-config --cflags opencv`
CPPFLAGS += $(COMMONFLAGS)
# Linker flags
LDFLAGS += `pkg-config --libs opencv`
LDFLAGS += -Wl
all: $(TARGET)

# Lists of object files
COBJS=$(CFILES:.c=.c.o)
CPPOBJS=$(CPPFILES:.cpp=.cpp.o)

$(TARGET): $(COBJS) $(CPPOBJS)
	$(CC) -g $(COBJS) $(CPPOBJS) $(LDFLAGS) -o $(TARGET)

# Clean
clean:
	$(RM) $(TARGET) *.o core.*

run:
	./$(TARGET) ../images/sea256x256.png out.jpg 90
	#./$(TARGET) ../images/sea16x16.png out.jpg 90
gdb:
	gdb --args ./$(TARGET) ../images/sea16x16.png out.jpg 90
pyc:
	~/software/bin/pyclewn -w left --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ../images/sea16x16.png out.jpg 90'
valg:
	valgrind --db-attach=yes --leak-check=yes --tool=memcheck --num-callers=16 --leak-resolution=high ./$(TARGET) ../images/sea16x16.png out.jpg 90

# Pattern rule for compiling C files
%.c.o: %.c 
	$(CC) $(CFLAGS) -c $< -o $@

# Pattern rule for compiling CPP files
%.cpp.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@;

#if you want to add the cpp files , please write the file dependencies here.
# Set file dependencies
ClContext.cpp.o: ClContext.cpp
ColorConv.cpp.o: ColorConv.cpp
main.cpp.o: main.cpp
