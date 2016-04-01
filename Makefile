sanwei-bing: H_Border.o h_Coord.o Inside.o Partition.o RWsgy.o DataTran.o testHFWI3D.o main.o
	mpicxx -o sanwei-bing H_Border.o h_Coord.o Inside.o Partition.o  RWsgy.o DataTran.o testHFWI3D.o main.o -lfftw3

H_Border.o: Partition.h H_Border.cpp
	mpicxx -c H_Border.cpp -lfftw3

h_Coord.o: Partition.h h_Coord.cpp
	mpicxx -c h_Coord.cpp -lfftw3

Inside.o: Partition.h Inside.cpp
	mpicxx -c Inside.cpp -lfftw3

Partition.o: Partition.cpp Partition.h H_Border.cpp h_Coord.cpp Inside.cpp
	mpicxx -c Partition.cpp -lfftw3

RWsgy.o: RWsgy.cpp RWsgy.h Partition.h
	mpicxx -c RWsgy.cpp -lfftw3

DataTran.o: DataTran.cpp DataTran.h Partition.h 
	mpicxx -c DataTran.cpp -lfftw3

testHFWI3D.o: testHFWI3D.cpp testHFWI3D.h Partition.h
	mpicxx -c testHFWI3D.cpp -lfftw3

main.o: main.cpp 
	mpicxx -c main.cpp -lfftw3

clean:
	rm -rf *.o
	rm -rf sanwei-bing
	ls *.txt | grep -v information.txt | xargs rm -rf 
