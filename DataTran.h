/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/
#ifndef DATATRAN_H
#define DATATRAN_H
#include "Partition.h"
#include <mpi.h>
#include "RWsgy.h"

#define TOP_TO_BOTTOM 0
#define LEFT_TO_RIGHT 1
#define BOTTOM_TO_TOP 2
#define RIGHT_TO_LEFT 3
#define FRONT_TO_BACK 4
#define BACK_TO_FRONT 5

#define STEP_U 0
#define STEP_V 6
#define STEP_K 12
#define STEP_W 18

#define STEP_VP 24

#define STEP_FEN 30
#define STEP_MAX 31
#define STEP_OBJ 32
#define STEP_MAX_RL 33
#define STEP_MIN_RL 34

#define STEP_SUMRESTRIAL 35
#define STEP_SUMRESCURR 36

#define STEP_GRAD 37

//#define UpdateVpPML 3
void dataTransport(float *data, const Partition& pt, int tag, int it, const MPI_Comm& mycomm);
void dataTransport_Vp(float *data, const Partition& pt, int tag, const AFDP3D &Pa, const MPI_Comm &mycomm);

void dataGather(float *data, const Partition& pt, int tag, const MPI_Comm &mycomm);

void copydatatobuf(float *data, float *buf, const Partition& pt, uint transportlen_side, int tag, int flag);
void copydatatobuf_Vp(float *data, float *buf, const Partition& pt, uint transportlen_side, int flag, const AFDP3D &Pa);
void copybuftodata(float *buf, float *data, const Partition& pt, uint transportlen_side, int tag, int flag);
#endif // DATATRAN_H
