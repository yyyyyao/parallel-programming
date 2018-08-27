CC      = g++
GDB	= gdb
CFLAGS  = -g -MMD -MP -O0
LDFLAGS = 
LIBS    = 
INCLUDE = -I ./include
SRC_DIR = ./src
OBJ_DIR = ./build
SOURCES = $(shell ls $(SRC_DIR)/*.c) 
OBJS    = $(subst $(SRC_DIR),$(OBJ_DIR), $(SOURCES:.c=.o))
TARGET  = matmul
DEPENDS = $(OBJS:.o=.d)

all: $(TARGET)

run:
	./$(TARGET)

$(TARGET): $(OBJS) $(LIBS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c 
	@if [ ! -d $(OBJ_DIR) ]; \
	    then echo "mkdir -p $(OBJ_DIR)"; mkdir -p $(OBJ_DIR); \
	    fi
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $< 

clean:
	$(RM) $(OBJS) $(TARGET) $(DEPENDS)

gdb:
	$(GDB) -tui ./$(TARGET)

-include $(DEPENDS)

.PHONY: all clean
