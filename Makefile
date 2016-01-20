sanwei-bing: H_Border.o h_Coord.o Inside.o Partition.o RWsgy.o DataTran.o testHFWI3D.o main.o
	mpicxx -o sanwei-bing H_Border.o h_Coord.o Inside.o Partition.o  RWsgy.o DataTran.o testHFWI3D.o main.o -lfftw3

H_Border.o: Partition.h H_Border.cpp
	mpicxx -c -std=c++11 H_Border.cpp -lfftw3

h_Coord.o: Partition.h h_Coord.cpp
	mpicxx -c -std=c++11 h_Coord.cpp -lfftw3

Inside.o: Partition.h Inside.cpp
	mpicxx -c -std=c++11 Inside.cpp -lfftw3

Partition.o: Partition.cpp Partition.h H_Border.cpp h_Coord.cpp Inside.cpp
	mpicxx -c -std=c++11 Partition.cpp -lfftw3

RWsgy.o: RWsgy.cpp RWsgy.h Partition.h
	mpicxx -c -std=c++11 RWsgy.cpp -lfftw3

DataTran.o: DataTran.cpp DataTran.h Partition.h 
	mpicxx -c -std=c++11 DataTran.cpp -lfftw3

testHFWI3D.o: testHFWI3D.cpp testHFWI3D.h Partition.h
	mpicxx -c -std=c++11 testHFWI3D.cpp -lfftw3

main.o: main.cpp 
	mpicxx -c -std=c++11 main.cpp -lfftw3

clean:
	rm -rf *.o
	rm -rf sanwei-bing
	rm -rf *.txt
