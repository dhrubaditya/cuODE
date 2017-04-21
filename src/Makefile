compile=nvcc 
copts=-DOUBLE -arch=sm_20 
objects = ode.o solve.o model.o Random.o SpecialFunctions.o cuda.o 
headers = input.h ode.h model.h Random.h SpecialFunctions.h CUDA.h 
#default: ode.cu solve.cu model.cu Random.cu SpecialFunctions.cu
#	${compile} ${copts} ode.cu solve.cu model.cu Random.cu SpecialFunctions.cu  -o ode.exe 
all: ${objects} ${headers}
	${compile} ${copts}  ${objects} -o ode.exe 
rand: rand.o Random.o 
	${compile} ${copts}  rand.o Random.o -o ode.exe 
rand.o: rand.cu
	${compile}  -arch=sm_20 -I. -dc rand.cu  
Random.o: Random.cu
	${compile}  -arch=sm_20 -I. -dc Random.cu  
%.o: %.cu
	${compile}  -arch=sm_20 -I. -dc $< -o $@
test: test.o Random.o SpecialFunctions.o cuda.o
	${compile} ${copts}  test.o Random.o SpecialFunctions.o cuda.o 
clean:
	rm -f *.o *.exe 
