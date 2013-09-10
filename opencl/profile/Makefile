CC = g++
RM = rm
LDFLAGS += -L/usr/local/cuda/lib64 -lOpenCL
CFLAGS += -I/usr/local/cuda/cuda/include

EXTENSION = cpp

TARGET = main
$(TARGET): $(TARGET).$(EXTENSION)
	$(CC) -g -o $(TARGET) $(LDFLAGS) $(CFLAGS) $(TARGET).$(EXTENSION)

clean:
	$(RM) $(TARGET)

run:
	./$(TARGET) 4096
