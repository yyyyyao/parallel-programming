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
	#./tga2jpg ./opencv-logo-150.png ./logo.jpg 80
	./tga2jpg ../jpeg-testing/jpgRefUpdate/sea16x16.png ./ref16x16.jpg 90

pyc:
	~/software/bin/pyclewn --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ./opencv-logo-150.png ./logo.jpg 80'
