sanwei-bing: H_Border.o h_Coord.o Inside.o Partition.o RWsgy.o DataTran.o testHFWI3D.o main.o
	mpicxx -o sanwei-bing H_Border.o h_Coord.o Inside.o Partition.o  RWsgy.o DataTran.o testHFWI3D.o main.o -lfftw3 -Wall

H_Border.o: Partition.h H_Border.cpp
	mpicxx -c H_Border.cpp -lfftw3 -Wall

h_Coord.o: Partition.h h_Coord.cpp
	mpicxx -c h_Coord.cpp -lfftw3 -Wall

Inside.o: Partition.h Inside.cpp
	mpicxx -c Inside.cpp -lfftw3 -Wall

Partition.o: Partition.cpp Partition.h H_Border.cpp h_Coord.cpp Inside.cpp
	mpicxx -c Partition.cpp -lfftw3 -Wall

RWsgy.o: RWsgy.cpp RWsgy.h Partition.h
	mpicxx -c RWsgy.cpp -lfftw3 -Wall

DataTran.o: DataTran.cpp DataTran.h Partition.h 
	mpicxx -c DataTran.cpp -lfftw3 -Wall

testHFWI3D.o: testHFWI3D.cpp testHFWI3D.h Partition.h
	mpicxx -c testHFWI3D.cpp -lfftw3 -Wall

main.o: main.cpp 
	mpicxx -c main.cpp -lfftw3 -Wall

clean:
	rm -rf *.o
	rm -rf sanwei-bing
	ls *.txt | grep -v information.txt | xargs rm -rf 
	rm -rf ./CalGrad/*
	rm -rf ./CalGrad_comu_/*
	rm -rf ./CalStepLength/*
	rm -rf ./CalStepLength_comu_/*
	rm -rf ./CalTrueWF/*
	rm -rf ./CalTrueWF_comu_/*
	rm -rf ./PreProcess_comu_/*
	rm -rf ./read_time/*
	rm -rf ./Total_time/*
