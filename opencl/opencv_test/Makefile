CC = g++
RM = rm
LDFLAGS += `pkg-config --libs opencv`
CFLAGS += `pkg-config --cflags opencv`
DEBUG += -g
EXTENSION = cpp

TARGET = blur
$(TARGET): $(TARGET).$(EXTENSION)
	$(CC) $(TARGET).$(EXTENSION) $(DEBUG) -o $(TARGET) $(LDFLAGS) $(CFLAGS)

clean:
	$(RM) $(TARGET)

run:
	./$(TARGET) ./dat.jpg 3
