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
	./tga2jpg ../images/sea256x256.png ./ref256x256.jpg 90

pyc:
	~/software/bin/pyclewn --cargs "-S ~/.pyclewn_mapping.vim" -a '--args ./$(TARGET) ./opencv-logo-150.png ./logo.jpg 80'
