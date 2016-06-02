./3D_FDTD_DWX: ./build/H_Border.o ./build/h_Coord.o ./build/Inside.o ./build/Partition.o ./build/RWsgy.o ./build/DataTran.o ./build/testHFWI3D.o ./build/main.o
	mpicxx -o 3D_FDTD_DWX ./build/H_Border.o ./build/h_Coord.o ./build/Inside.o ./build/Partition.o  ./build/RWsgy.o ./build/DataTran.o ./build/testHFWI3D.o ./build/main.o -lfftw3 -Wall

./build/H_Border.o: ./include/Partition.h ./util/H_Border.cpp
	mpicxx -c -I ./include ./util/H_Border.cpp -o $@ -lfftw3 -Wall -g

./build/h_Coord.o: ./include/Partition.h ./util/h_Coord.cpp
	mpicxx -c -I ./include  ./util/h_Coord.cpp -o $@ -lfftw3 -Wall -g

./build/Inside.o: ./include/Partition.h ./util/Inside.cpp
	mpicxx -c -I ./include  ./util/Inside.cpp -o $@ -lfftw3 -Wall -g

./build/Partition.o: ./util/Partition.cpp ./include/Partition.h ./util/H_Border.cpp ./util/h_Coord.cpp ./util/Inside.cpp
	mpicxx -c -I ./include ./util/Partition.cpp -o $@ -lfftw3 -Wall -g

./build/RWsgy.o: ./util/RWsgy.cpp ./include/RWsgy.h ./include/Partition.h
	mpicxx -c -I ./include ./util/RWsgy.cpp -o $@ -lfftw3 -Wall -g

./build/DataTran.o: ./util/DataTran.cpp ./include/DataTran.h ./include/Partition.h 
	mpicxx -c -I ./include ./util/DataTran.cpp -o $@ -lfftw3 -Wall -g

./build/testHFWI3D.o: ./util/testHFWI3D.cpp ./include/testHFWI3D.h ./include/Partition.h
	mpicxx -c -I ./include ./util/testHFWI3D.cpp -o $@ -lfftw3 -Wall -g

./build/main.o: main.cpp ./include/testHFWI3D.h
	mpicxx -c -I ./include main.cpp -o $@ -lfftw3 -Wall -g

clean:
	rm -rf ./build/*.o
	rm -rf 3D_FDTD_DWX
	ls *.txt | grep -v information.txt | xargs rm -rf 
	rm -rf ./result/CalGrad/*
	rm -rf ./result/CalGrad_comu_/*
	rm -rf ./result/CalStepLength/*
	rm -rf ./result/CalStepLength_comu_/*
	rm -rf ./result/CalTrueWF/*
	rm -rf ./result/CalTrueWF_comu_/*
	rm -rf ./result/PreProcess_comu_/*
	rm -rf ./result/read_time/*
	rm -rf ./result/Total_time/*
