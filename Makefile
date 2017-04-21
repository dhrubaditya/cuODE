# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
compile=nvcc -dc -DOUBLE -arch=sm_20
#MODEL=VB_dynamo
MODEL=harmonic
#MODEL=particle_in_mag
CUDA=CUDA
# ---------------------------------------------------------------------
default:
	(cd src/;ln -sf models/${MODEL}.h model.h; ln -sf models/${MODEL}.cu model.cu; ln -sf CUDA/${CUDA}.h CUDA.h; ln -sf CUDA/cuda.cu .; ln -sf CUDA/mycomplex.h .; ln -sf toolbox/Random.h .; ln -sf toolbox/SpecialFunctions.h .; ln -sf toolbox/Random.cu .; ln -sf toolbox/SpecialFunctions.cu .; ln -sf ../input.h .; make; mv ode.exe ..)
clean:
	cd src/; rm -f *.o *.mod *.exe
