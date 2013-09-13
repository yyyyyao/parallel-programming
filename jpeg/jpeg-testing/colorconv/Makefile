CC = g++
RM = rm
LDFLAGS += `pkg-config --libs opencv`
CFLAGS += `pkg-config --cflags opencv`
DEBUG += -g
EXTENSION = cpp

TARGET = colorConv
$(TARGET): $(TARGET).$(EXTENSION)
	$(CC) $(DEBUG) -o $(TARGET) $(LDFLAGS) $(CFLAGS) $(TARGET).$(EXTENSION)

clean:
	$(RM) $(TARGET)

run:
	./$(TARGET) ../images/sea256x256.png 3
gdb:
	gdb --args ./$(TARGET) ../images/sea256x256.png 3
pyc:
	~/software/bin/pyclewn --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ../images/sea256x256.png'
