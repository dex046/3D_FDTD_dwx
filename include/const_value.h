/******************************************
 * author:dwx
 ******************************************/
#ifndef CONST_VALUE_H
#define CONST_VALUE_H

typedef     unsigned int    uint;
const float e = 12.0f;

const char * const TrueSg = "./result/TrueSG.sgy";
const char * const GradVp = "./result/GradientVp.sgy";
const char * const InvertedVp = "./result/InvertedVp.sgy";

const char * const TrueVp = "./data/TrueVp-3D.sgy";
const char * const InitVp= "./data/InitVp-3D.sgy";

const uint cpu_x = 1, cpu_y = 8, cpu_z = 1;
const uint shot_x = 4, shot_y = 1, shot_z = 2;
#endif // CONST_VALUE_H
