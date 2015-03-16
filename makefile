#make file for brownian_dynamics package

CC = g++ -std=c++0x
CC1 = icpc -std=c++0x
OMP = -fopenmp
NVC = nvcc -arch=sm_20 -c
NVL = nvcc -lcudart -o

reg: test.cpp
	${CC} -c test.cpp
	${CC} -c proteins.cpp
	${CC} -c brownsys.cpp
	${CC} -o run test.o proteins.o brownsys.o 
	rm -f *.o

par: test.cpp
	${CC} -c test.cpp ${OMP}
	${CC} -c proteins.cpp ${OMP} 
	${CC} -c brownsys.cpp ${OMP}
	${NVC} neighbors.cu
	${NVL} run test.o proteins.o brownsys.o neighbors.o
