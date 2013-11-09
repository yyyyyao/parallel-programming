DEBUG ?= 1

CXX = g++

OPENCV_LIB += `pkg-config --libs opencv`
OPENCV_INCLUDE += `pkg-config --cflags opencv`

TARGET = jpeg-enc

OPENCL_INCLUDE = -I /usr/local/cuda/include
OPENCL_LIB = -lOpenCL

CFLAGS := $(OPENCL_INCLUDE_PATH) $(OPENCV_INCLUDE)
LDFLAGS := $(OPENCL_LIB) $(OPENCV_LIB)

SRCDIR = .
OBJDIR = ./obj

CFILES := 
CXXFILES := main.cpp

all: $(TARGET)

COBJS = $(CFILES:.c=.c.o)
CXXOBJS = $(CXXFILES:.cpp=.cpp.o)

$(TARGET) : $(COBJS) $(CXXOBJS)
	$(CXX) -g $(COBJS) $(CXXOBJS) $(CFLAGS) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f *.o $(TARGET)

run:
	./$(TARGET) 256 256 3 ./darvish_256x256.raw
gdb:
	~/software/bin/gdb --args ./$(TARGET) 256 256 3 ./darvish_256x256.raw

%.c.o: %.c
	$(CC) $(CFLAGS) -g -c $< -o $@

%.cpp.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS) -g -c $< -o $@

# Set file dependencies
main.cpp.o: main.cpp
