# This acts as a makefile for the ptx files

# This only needs to be run once. ptx and cu files are saved to repository. 
# Need to have a CUDA toolkit and VS2013 installed.

# Adding system path for nvcc to compile with nvcc
#setenv('PATH', [getenv('PATH') ';C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\bin']);
# Adding system path for VS2013 to compile with cl
#setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin']);
#setenv('PATH', [getenv('PATH') ';/usr/local/cuda-6.5/bin']);
#
#[a, b] = system('nvcc -ptx cuda_gaussMLEv2\cuda_gaussMLEv2.cu -o  ../../ptx/cuda_gaussMLEv2.ptx')
#copyfile('cuda_gaussMLEv2\cuda_gaussMLEv2.cu','../../ptx/cuda_gaussMLEv2.cu')

BASE = cuda_gaussBlobROIStack

SRCS = $(BASE).cu *.cuh *.h
OBJS = $(BASE).ptx

NVCC = nvcc

all: $(OBJS)
	mv -f $(BASE).ptx ../../../ptx
	cp -f $(BASE).cu  ../../../ptx

clean:
	rm -f $(OBJS) ../../../ptx/$(BASE).cu ../../../ptx/$(BASE).ptx

.SUFFIXES: $(SUFFIXES) .cu .ptx

.cu.ptx:
	$(NVCC) -ptx $<
