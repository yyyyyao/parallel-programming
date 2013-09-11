# Use constant memory for huffman coder tables?
HUFFMAN_CODER_TABLES_IN_CONSTANT ?= 1
# Debug
DEBUG ?= 1
# CUDA install path
CUDA_INSTALL_PATH ?= /usr/local/cuda

# Target executable
TARGET := gpujpeg

# C files
CFILES := \
    gpujpeg_common.c \
    gpujpeg_encoder.c \
    gpujpeg_decoder.c \
    gpujpeg_table.c \
    gpujpeg_dct_cpu.c \
    gpujpeg_huffman_cpu_encoder.c \
    gpujpeg_huffman_cpu_decoder.c \
    gpujpeg_writer.c \
    gpujpeg_reader.c \
    main.c
# CUDA files
CUFILES := \
    gpujpeg_dct_gpu.cu \
    gpujpeg_preprocessor.cu \
    gpujpeg_huffman_gpu_encoder.cu \
    gpujpeg_huffman_gpu_decoder.cu

# Compilers
CC := gcc -fPIC
LINK := g++ -fPIC
NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc -Xcompiler -fPIC

# Debug
ifeq ($(DEBUG),1)
    COMMONFLAGS += -g -D_DEBUG -O0
    NVCCFLAGS += -G
else
    COMMONFLAGS += -O2
endif

# Common flags
COMMONFLAGS += -I. -I$(CUDA_INSTALL_PATH)/include
# C flags
CFLAGS += $(COMMONFLAGS) -std=c99
# CUDA flags
NVCCFLAGS += $(COMMONFLAGS) \
	--ptxas-options="-v" \
	-gencode arch=compute_20,code=sm_20 #\
	-gencode arch=compute_11,code=sm_11 \
	-gencode arch=compute_10,code=sm_10
# Linker flags
LDFLAGS += -shared -Wl,-soname,$(SONAME)

# Other Flags
ifeq ($(HUFFMAN_CODER_TABLES_IN_CONSTANT),1)
    COMMONFLAGS += -DGPUJPEG_HUFFMAN_CODER_TABLES_IN_CONSTANT
endif

# Do 32bit vs. 64bit setup
LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
    # 64bit
    LDFLAGS += -L$(CUDA_INSTALL_PATH)/lib64
else
    # 32bit
    LDFLAGS += -L$(CUDA_INSTALL_PATH)/lib
endif
LDFLAGS += -lcudart

all: $(TARGET)

# Lists of object files
COBJS=$(CFILES:.c=.c.o)
CUOBJS=$(CUFILES:.cu=.cu.o)

$(TARGET): $(COBJS) $(CUOBJS)
	$(NVCC) -g $(COBJS) $(CUOBJS) -o $(TARGET)

# Clean
clean:
	rm -f *.o $(TARGET) $(SONAME) $(LIBRARY_NAME)
	rm -f *.i *.ii 
	rm -f *.cudafe1.c *.cudafe1.cpp *.cudafe1.gpu *.cudafe1.stub.c
	rm -f *.cudafe2.c *.cudafe2.gpu *.cudafe2.stub.c
	rm -f *.fatbin *.fatbin.c *.ptx *.hash *.cubin *.cu.cpp

run:
	#./gpujpeg -v --encode --size=3456x5184 --quality=80 ./female.raw female.jpg
	./gpujpeg -v --encode --size=600x338 --quality=80 ./gakuryoku600x338.raw gakuryoku.jpg

pyc:
	./gpujpeg -v --encode --size=5184x3456 --quality=80 ./female3456x5184.raw female.jpg
# Set suffix for CUDA files
.SUFFIXES: .cu

# Pattern rule for compiling C files
%.c.o: %.c 
	$(CC) $(CFLAGS) -c $< -o $@

# Pattern rule for compiling CUDA files
%.cu.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@;

# Set file dependencies
gpujpeg_common.c.o: gpujpeg_common.c libgpujpeg/gpujpeg_common.h
gpujpeg_encoder.c.o: gpujpeg_encoder.c libgpujpeg/gpujpeg_encoder.h
gpujpeg_decoder.c.o: gpujpeg_decoder.c libgpujpeg/gpujpeg_decoder.h
gpujpeg_table.c.o: gpujpeg_table.c libgpujpeg/gpujpeg_table.h
gpujpeg_preprocessor.cu.o: gpujpeg_preprocessor.cu gpujpeg_preprocessor.h gpujpeg_colorspace.h
gpujpeg_dct_cpu.c.o: gpujpeg_dct_cpu.c gpujpeg_dct_cpu.h
gpujpeg_dct_gpu.cu.o: gpujpeg_dct_gpu.cu gpujpeg_dct_gpu.h
gpujpeg_huffman_cpu_encoder.c.o: gpujpeg_huffman_cpu_encoder.c gpujpeg_huffman_cpu_encoder.h
gpujpeg_huffman_gpu_encoder.cu.o: gpujpeg_huffman_gpu_encoder.cu gpujpeg_huffman_gpu_encoder.h
gpujpeg_huffman_cpu_decoder.c.o: gpujpeg_huffman_cpu_decoder.c gpujpeg_huffman_cpu_decoder.h
gpujpeg_huffman_gpu_decoder.cu.o: gpujpeg_huffman_gpu_decoder.cu gpujpeg_huffman_gpu_decoder.h
gpujpeg_writer.c.o: gpujpeg_writer.c libgpujpeg/gpujpeg_writer.h
gpujpeg_reader.c.o: gpujpeg_reader.c libgpujpeg/gpujpeg_reader.h
main.c.o: main.c
