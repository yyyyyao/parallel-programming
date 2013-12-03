CC = g++
RM = rm
LDFLAGS += -L/usr/local/cuda/lib64
LDFLAGS += -L/usr/lib64/nvidia
LDFLAGS += -lOpenCL
LDFLAGS += `pkg-config --libs opencv`
CFLAGS += -I/usr/local/cuda/include
CFLAGS += `pkg-config --cflags opencv`
EXTENSION = cpp

TARGET = main
$(TARGET): $(TARGET).$(EXTENSION)
	$(CC) -o $(TARGET) $(LDFLAGS) $(CFLAGS) $(TARGET).$(EXTENSION)

clean:
	$(RM) $(TARGET)

run:
	./$(TARGET)
