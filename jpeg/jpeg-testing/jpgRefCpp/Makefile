CC = g++
RM = rm
LDFLAGS += `pkg-config --libs opencv`
CFLAGS += `pkg-config --cflags opencv`
DEBUG += -g
EXTENSION = cpp

TARGET = jpgRef
$(TARGET): $(TARGET).$(EXTENSION)
	$(CC) $(DEBUG) -o $(TARGET) $(LDFLAGS) $(CFLAGS) $(TARGET).$(EXTENSION)

clean:
	$(RM) $(TARGET)

run:
	./$(TARGET) ./sea16x16.png res.jpg 90
	#./$(TARGET) ../images/sea256x256.png
gdb:
	gdb --args ./$(TARGET) ./sea16x16.png 3
pyc:
	~/software/bin/pyclewn -w left --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ./sea16x16.png res.jpg 90'
	#~/software/bin/pyclewn -w left --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ../images/sea256x256.png'
valg:
	valgrind --db-attach=yes --leak-check=yes --tool=memcheck --num-callers=16 --leak-resolution=high ./$(TARGET) ../images/sea256x256.png res.jpg 90
