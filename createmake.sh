#CC=mpicxx
#C_FLAGS="-Wall -lfftw3 -g -I $(pwd)/include/"
#export CC
#export C_FLAGS

currentdir=$(pwd)

#echo $C_FLAGS
#echo $currentdir;
#echo $1

cd $1
rm Makefile
echo '' >> Makefile

echo -n "all: " >> Makefile
for file in $(ls *.cpp)
do
	obj=$(echo ${file} | cut -d '.' -f 1 | sed 's/^/..\/build\//g')
	echo '\' >> Makefile
	echo -n "${obj}.o " >> Makefile
done

echo '' >> Makefile
echo '' >> Makefile

for file in $(ls *.cpp)
do
	name=$(echo ${file} | cut -d '.' -f 1)
	${CC} -MM -I ../include/ ${file} | sed "s/$name.o/..\/build\/$name.o/g" >> Makefile
	echo -e "\t${CC} ${C_FLAGS} -c \$< -o \$@" >> Makefile
#obj=$
done
cd $currentdir

