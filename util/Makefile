
all: \
../build/DataTran.o \
../build/H_Border.o \
../build/h_Coord.o \
../build/Inside.o \
../build/Partition.o \
../build/RWsgy.o \
../build/testHFWI3D.o 

../build/DataTran.o: DataTran.cpp ../include/DataTran.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h ../include/Partition.h ../include/RWsgy.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
../build/H_Border.o: H_Border.cpp ../include/Partition.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
../build/h_Coord.o: h_Coord.cpp ../include/Partition.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
../build/Inside.o: Inside.cpp ../include/Partition.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
../build/Partition.o: Partition.cpp ../include/Partition.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
../build/RWsgy.o: RWsgy.cpp ../include/RWsgy.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h ../include/Partition.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
../build/testHFWI3D.o: testHFWI3D.cpp ../include/testHFWI3D.h \
 ../include/DataTran.h ../include/common.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpio.h \
 /usr/include/mpich-x86_64/mpi.h /usr/include/mpich-x86_64/mpicxx.h \
 ../include/const_value.h ../include/Partition.h ../include/RWsgy.h
	mpicxx -Wall -lfftw3 -g -I /home/xtay/3D_FDTD_dwx/include/ -c $< -o $@
