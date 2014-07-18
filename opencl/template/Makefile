CC = g++
RM = rm
#LDFLAGS += -L/usr/local/cuda/lib64
LDFLAGS += -L/opt/intel/opencl-1.2-4.4.0.117/lib64
LDFLAGS += -lOpenCL
#LDFLAGS += `pkg-config --libs opencv`
#CFLAGS += -I/usr/local/cuda/include
CFLAGS += -I/opt/intel/opencl-1.2-4.4.0.117/include
CFLAGS += -I/usr/local/include
#CFLAGS += `pkg-config --cflags opencv`
EXTENSION = cpp

TARGET = main
$(TARGET): $(TARGET).$(EXTENSION)
	$(CC) $(TARGET).$(EXTENSION) -o $(TARGET) $(CFLAGS) $(LDFLAGS)

clean:
	$(RM) $(TARGET)

run:
	./$(TARGET)
