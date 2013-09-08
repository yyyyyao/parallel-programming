CC = g++
LDFLAGS +=
CFLAGS += -fdump-rtl-expand
TARGET = tga2jpg
EXTENSION = cpp
OBJS = jpgd.o jpge.o tga2jpg.o timer.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -g $(OBJS) $(CFLAGS) $(LDFLAGS) -o $(TARGET)

.cpp.o:
	$(CC) -g -c $(CFLAGS) $(LDFLAGS) $<

.c.o:
	$(CC) -g -c $(CFLAGS) $(LDFLAGS) $<

clean:
	rm -rf $(TARGET) *.o

run:
	./$(TARGET) ../jpeg-testing/images/sea16x16.png ./sea16x16.jpg 90

pyc:
	~/software/bin/pyclewn -w left --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ../jpeg-testing/images/sea16x16.png ./sea16x16_.jpg 90'
