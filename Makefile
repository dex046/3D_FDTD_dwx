CC=mpicxx
C_FLAGS=-Wall -lfftw3 -g -I $(shell pwd)/include/###shell执行shell调用,shell脚本用bash
export CC
export C_FLAGS

CURRENT_DIR=.
UTIL_DIR=util
INCLUDE_DIR=include
BUILD_DIR=build
RESULT_DIR=result

RESULT=3D_FDTD_DWX

SRC:=$(wildcard ${UTIL_DIR}/*.cpp)
SRC+=$(wildcard ${CURRENT_DIR}/*.cpp)

OBJ=$(SRC:.cpp=.o)
OBJ:=$(notdir $(OBJ))
OBJ:=$(addprefix ${BUILD_DIR}/, ${OBJ})

HEADERS=$(wildcard $(INCLUDE_DIR)/*.h)

#vpath %.h ${currentdir}/include/
#vpath %.cpp ${currentdir}/util/

B=$(A)
A=later
all: dir $(RESULT)

dir:
	@if [ ! -d $(BUILD_DIR) ]; then\
		mkdir $(BUILD_DIR);\
	fi
	@if [ ! -d ${RESULT_DIR} ]; then\
		mkdir ${RESULT_DIR};\
	fi

$(RESULT): makeutil $(BUILD_DIR)/main.o
	$(CC) $(C_FLAGS) $(OBJ) -o $(RESULT)

#excute Makefile in the util dir
makeutil: $(UTIL_DIR)/Makefile
	$(MAKE) -C $(UTIL_DIR)

$(BUILD_DIR)/main.o: $(CURRENT_DIR)/main.cpp
	$(CC) ${C_FLAGS} -c $< -o $@

#create or recreate the Makefile while new files added into the util dir
$(UTIL_DIR)/Makefile: $(UTIL_DIR)
	bash createmake.sh $(UTIL_DIR);

.PHONY:clean
clean:
	-rm -r ./build/*.o
	-rm -r 3D_FDTD_DWX
	-rm ls *.txt | grep -v information.txt | xargs rm -rf 
	-rm -r ./result/CalGrad/*
	-rm -r ./result/CalGrad_comu_/*
	-rm -r ./result/CalStepLength/*
	-rm -r ./result/CalStepLength_comu_/*
	-rm -r ./result/CalTrueWF/*
	-rm -r ./result/CalTrueWF_comu_/*
	-rm -r ./result/PreProcess_comu_/*
	-rm -r ./result/read_time/*
	-rm -r ./result/Total_time/*
veryclean:
	-rm -rf ./build/*.o
	-rm -rf 3D_FDTD_DWX
	-rm -rf ./util/Makefile
	ls *.txt | grep -v information.txt | xargs rm -rf 
	-rm -rf ./result/CalGrad/*
	-rm -rf ./result/CalGrad_comu_/*
	-rm -rf ./result/CalStepLength/*
	-rm -rf ./result/CalStepLength_comu_/*
	-rm -rf ./result/CalTrueWF/*
	-rm -rf./result/CalTrueWF_comu_/*
	-rm -rf ./result/PreProcess_comu_/*
	-rm -rf ./result/read_time/*
	-rm -rf ./result/Total_time/*
