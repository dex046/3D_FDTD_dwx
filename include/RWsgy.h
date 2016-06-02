/************************************************************************/
/*
* 文件信息：本程序包含使用C语言来读写Sgy文件的相关变量声明

* 文件名: RWsgy.h

* 版本：0.1

* 作者：高照奇

* 日期：2013-12-08
*/
/************************************************************************/
#ifndef __RWSGY_H__
#define __RWSGY_H__

#pragma once
#include "common.h"
#include "Partition.h"

#define     ROOT_ID     0
#define     WRITE_ALL   0
#define     WRITE_INTER 1
/*
SEG-Y总格式：
    |卷头	|		|道头1 数据1|道头2 数据2
    3200	400		240			240
*/

// 卷头前3200个字节是EBCDIC码，其转换成ASCII码可以通过查表实现
const unsigned char E2A[256]={
    0, 1, 2, 3,156, 9,134,127,151,141,142, 11, 12, 13, 14, 15, 16, 17, 18, 19,157,133,
    8,135, 24, 25,146,143, 28, 29, 30, 31, 128,129,130,131,132,10,23,
    27,136,137,138,139,140, 5, 6, 7, 144,145, 22,147,148,149,150, 4,152,153,154,155,
    20, 21,158, 26, 32,160,161,162,163,164,165,166,167,168, 91, 46, 60, 40, 43, 33,
    38,169,170,171,172,173,174,175,176,177, 93, 36, 42, 41, 59, 94, 45,
    47,178,179,180,181,182,183,184,185,124, 44, 37, 95, 62, 63,
    186,187,188,189,190,191,192,193,194, 96, 58, 35, 64, 39, 61, 34, 195, 97, 98,
    99,100,101,102,103,104,105,196,197,198,199,200,201,
    202,106,107,108,109,110,111,112,113,114,203,204,205,206,207,
    208,209,126,115,116,117,118,119,120,121,122,210,211,212,213,214,215,
    216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231, 123, 65, 66, 67,
    68, 69, 70,71, 72, 73,232,233,234,235,236,237, 125, 74, 75, 76,
    77,78,79,80,81,82,238,239,240,241,242,243, 92,159, 83, 84, 85, 86, 87, 88, 89,
    90,244,245,246,247,248,249, 48,49, 50, 51, 52, 53, 54, 55, 56,
    57,250,251,252,253,254,255
};
// SEGY卷头文件400字节内容

typedef struct  Struct_reelb400
{    /* bhed - binary header */
    int jobid;  /* job identification number */
                // 任务识别码
    int lino;   /* line number (only one line per reel) */
                // 测线号
    int reno;   /* reel number */
                // 卷号
    short ntrpr;    /* number of data traces per record */
                // 每个记录的道数
    short nart; /* number of auxiliary traces per record */
                // 每个记录的辅助道数
    unsigned short hdt; /* sample interval in micro secs for this reel */
                        // 该卷采样间隔
    unsigned short dto; /* same for original field recording */
                        // 原始记录的采样间隔
    unsigned short hns; /* number of samples per trace for this reel */
                        // 该卷每道的样点数
    unsigned short nso; /* same for original field recording */
                        // 原始记录每道的采样点
    short format;   /* data sample format code:
                    // 数据采样的记录格式
                    1 = floating point (4 bytes)
                    2 = fixed point (4 bytes)
                    3 = fixed point (2 bytes)
                    4 = fixed point w/gain code (4 bytes) */
    short fold; /* CDP fold expected per CDP ensemble */
                // CDP覆盖次数
    short tsort;    /* trace sorting code:
                    // 道分选码
                    1 = as recorded (no sorting)
                    2 = CDP ensemble
                    3 = single fold continuous profile
                    4 = horizontally stacked */
    short vscode;   /* vertical sum code:
                    // 垂直和码
                    1 = no sum
                    2 = two sum ...
                    N = N sum (N = 32,767) */
    short hsfs; /* sweep frequency at start */
                // 起始扫描频率
    short hsfe; /* sweep frequency at end */
                // 终止扫描频率
    short hslen;    /* sweep length (ms) */
                    // 扫描时窗
    short hstyp;    /* sweep type code:
                    // 扫描类型
                    1 = linear
                    2 = parabolic
                    3 = exponential
                    4 = other */
    short schn; /* trace number of sweep channel */
                // 扫描频带道数
    short hstas;    /* sweep trace taper length at start if
                    tapered (the taper starts at zero time
                    and is effective for this length) */
                    // 斜坡扫描起始时间
    short hstae;    /* sweep trace taper length at end (the ending
                    taper starts at sweep length minus the taper
                    length at end) */
                    // 斜坡扫描终止时间
    short htatyp;   /* sweep trace taper type code:
                    1 = linear                 2 = cos-squared
                    3 = other */
                    // 斜坡扫描类型
    short hcorr;    /* correlated data traces code:
                    1 = no
                    2 = yes */
                    // 相关数据道
    short bgrcv;    /* binary gain recovered code:
                    1 = yes
                    2 = no */
                    // 增益恢复
    short rcvm; /* amplitude recovery method code:
                1 = none
                2 = spherical divergence
                3 = AGC
                4 = other */
                // 振幅恢复方法
    short mfeet;    /* measurement system code:
                    1 = meters
                    2 = feet */
                    // 测量系统
    short polyt;    /* impulse signal polarity code:
                    1 = increase in pressure or upward
                    geophone case movement gives
                    negative number on tape
                    2 = increase in pressure or upward
                    geophone case movement gives
                    positive number on tape */
                    // 激发信号极性
    short vpol; /* vibratory polarity code:
                code    seismic signal lags pilot by
                1   337.5 to  22.5 degrees
                2    22.5 to  67.5 degrees
                3    67.5 to 112.5 degrees
                4   112.5 to 157.5 degrees
                5   157.5 to 202.5 degrees
                6   202.5 to 247.5 degrees
                7   247.5 to 292.5 degrees
                8   293.5 to 337.5 degrees */
                // 振动极性
    short hunass[170];  /* unassigned */
}Struct_bheader;


struct Struct_Su240
{
    int tracl; /* Trace sequence number within line
               --numbers continue to increase if the
               same line continues across multiple
               SEG Y files.*/
    int tracr; /* Trace sequence number within SEG Y file
               ---each file starts with trace sequence
               one*/
    int fldr; /* Original field record number */
    int tracf; /* Trace number within original field record */
    int ep;  /* energy source point number
             ---Used when more than one record occurs
             at the same effective surface location.*/
    int cdp; /* Ensemble number (i.e. CDP, CMP, CRP,...) */
    int cdpt; /* trace number within the ensemble
              ---each ensemble starts with trace number one.*/
    short trid; /* trace identification code:
                -1 = Other
                0 = Unknown
                1 = Seismic data
                2 = Dead
                3 = Dummy
                4 = Time break
                5 = Uphole
                6 = Sweep
                7 = Timing
                8 = Water break
                9 = Near-field gun signature
                10 = Far-field gun signature
                11 = Seismic pressure sensor
                12 = Multicomponent seismic sensor
                - Vertical component
                13 = Multicomponent seismic sensor
                - Cross-line component
                14 = Multicomponent seismic sensor
                - in-line component
                15 = Rotated multicomponent seismic sensor
                - Vertical component
                16 = Rotated multicomponent seismic sensor
                - Transverse component
                17 = Rotated multicomponent seismic sensor     - Radial component
                18 = Vibrator reaction mass
                19 = Vibrator baseplate
                20 = Vibrator estimated ground force
                21 = Vibrator reference
                22 = Time-velocity pairs
                23 ... N = optional use
                (maximum N = 32,767)

                Following are CWP id flags:

                109 = autocorrelation
                110 = Fourier transformed - no packing
                xr[0],xi[0], ..., xr[N-1],xi[N-1]
                111 = Fourier transformed - unpacked Nyquist
                xr[0],xi[0],...,xr[N/2],xi[N/2]
                112 = Fourier transformed - packed Nyquist
                even N:
                xr[0],xr[N/2],xr[1],xi[1], ...,
                xr[N/2 -1],xi[N/2 -1]
                (note the exceptional second entry)
                odd N:
                xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
                xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
                (note the exceptional second & last entries)
                113 = Complex signal in the time domain
                xr[0],xi[0], ..., xr[N-1],xi[N-1]
                114 = Fourier transformed - amplitude/phase
                a[0],p[0], ..., a[N-1],p[N-1]
                115 = Complex time signal - amplitude/phase
                a[0],p[0], ..., a[N-1],p[N-1]
                116 = Real part of complex trace from 0 to Nyquist
                117 = Imag part of complex trace from 0 to Nyquist
                118 = Amplitude of complex trace from 0 to Nyquist
                119 = Phase of complex trace from 0 to Nyquist
                121 = Wavenumber time domain (k-t)
                122 = Wavenumber frequency (k-omega)
                123 = Envelope of the complex time trace
                124 = Phase of the complex time trace
                125 = Frequency of the complex time trace
                130 = Depth-Range (z-x) traces
                143 = Seismic Data, Vertical Component
                144 = Seismic Data, Horizontal Component 1
                145 = Seismic Data, Horizontal Component 2     146 = Seismic Data, Radial Component
                147 = Seismic Data, Transverse Component
                201 = Seismic data packed to bytes (by supack1)
                202 = Seismic data packed to 2 bytes (by supack2)
                */
    short nvs; /* Number of vertically summed traces yielding
               this trace. (1 is one trace,
               2 is two summed traces, etc.)*/
    short nhs; /* Number of horizontally summed traces yielding
               this trace. (1 is one trace
               2 is two summed traces, etc.)*/
    short duse; /* Data use:
                1 = Production
                2 = Test*/
    int offset; /* Distance from the center of the source point
                to the center of the receiver group
                (negative if opposite to direction in which
                the line was shot).*/
    int gelev; /* Receiver group elevation from sea level
               (all elevations above the Vertical datum are
               positive and below are negative).*/
    int selev; /* Surface elevation at source. */
    int sdepth; /* Source depth below surface (a positive number). */
    int gdel; /* Datum elevation at receiver group. */
    int sdel; /* Datum elevation at source. */
    int swdep; /* Water depth at source. */
    int gwdep; /* Water depth at receiver group. */
    short scalel; /* Scalar to be applied to the previous 7 entries
                  to give the real value.
                  Scalar = 1, +10, +100, +1000, +10000.
                  If positive, scalar is used as a multiplier,
                  if negative, scalar is used as a divisor.*/
    short scalco; /* Scalar to be applied to the next 4 entries
                  to give the real value.
                  Scalar = 1, +10, +100, +1000, +10000.
                  If positive, scalar is used as a multiplier,
                  if negative, scalar is used as a divisor.*/
    int  sx; /* Source coordinate - X */
    int  sy; /* Source coordinate - Y */
    int  gx; /* Group coordinate - X */
    int  gy; /* Group coordinate - Y */
    short counit; /* Coordinate units: (for previous 4 entries and
                  for the 7 entries before scalel)
                  1 = Length (meters or feet)       2 = Seconds of arc
                  3 = Decimal degrees
                  4 = Degrees, minutes, seconds (DMS)

                  In case 2, the X values are longitude and
                  the Y values are latitude, a positive value designates
                  the number of seconds east of Greenwich
                  or north of the equator

                  In case 4, to encode +-DDDMMSS
                  counit = +-DDD*10^4 + MM*10^2 + SS,
                  with scalco = 1. To encode +-DDDMMSS.ss
                  counit = +-DDD*10^6 + MM*10^4 + SS*10^2
                  with scalco = -100.*/
    short wevel; /* Weathering velocity. */
    short swevel; /* Subweathering velocity. */
    short sut; /* Uphole time at source in milliseconds. */
    short gut; /* Uphole time at receiver group in milliseconds. */
    short sstat; /* Source static correction in milliseconds. */
    short gstat; /* Group static correction  in milliseconds.*/
    short tstat; /* Total static applied  in milliseconds.
                 (Zero if no static has been applied.) */
    short laga; /* Lag time A, time in ms between end of 240-
                byte trace identification header and time
                break, positive if time break occurs after
                end of header, time break is defined as
                the initiation pulse which maybe recorded
                on an auxiliary trace or as otherwise
                specified by the recording system */
    short lagb; /* lag time B, time in ms between the time break
                and the initiation time of the energy source,
                may be positive or negative */
    short delrt; /* delay recording time, time in ms between
                 initiation time of energy source and time
                 when recording of data samples begins
                 (for deep water work if recording does not
                 start at zero time) */
    short muts; /* mute time--start */
    short mute; /* mute time--end */
    unsigned short ns; /* number of samples in this trace */
    unsigned short dt; /* sample interval; in micro-seconds */
    short gain; /* gain type of field instruments code:
                1 = fixed
                2 = binary     3 = floating point
                4 ---- N = optional use */
    short igc; /* instrument gain constant */
    short igi; /* instrument early or initial gain */
    short corr; /* correlated:
                1 = no
                2 = yes */
    short sfs; /* sweep frequency at start */
    short sfe; /* sweep frequency at end */
    short slen; /* sweep length in ms */
    short styp; /* sweep type code:
                1 = linear
                2 = cos-squared
                3 = other */
    short stas; /* sweep trace length at start in ms */
    short stae; /* sweep trace length at end in ms */
    short tatyp; /* taper type: 1=linear, 2=cos^2, 3=other */
    short afilf; /* alias filter frequency if used */
    short afils; /* alias filter slope */
    short nofilf; /* notch filter frequency if used */
    short nofils; /* notch filter slope */
    short lcf; /* low cut frequency if used */
    short hcf; /* high cut frequncy if used */
    short lcs; /* low cut slope */
    short hcs; /* high cut slope */
    short year; /* year data recorded */
    short day; /* day of year */
    short hour; /* hour of day (24 hour clock) */
    short minute; /* minute of hour */
    short sec; /* second of minute */
    short timbas; /* time basis code:
                  1 = local
                  2 = GMT
                  3 = other */
    short trwf; /* trace weighting factor, defined as 1/2^N
                volts for the least sigificant bit */
    short grnors; /* geophone group number of roll switch
                  position one */
    short grnofr; /* geophone group number of trace one within
                  original field record */
    short grnlof; /* geophone group number of last trace within
                  original field record */
    short gaps; /* gap size (total number of groups dropped) */
    short otrav; /* overtravel taper code:     1 = down (or behind)
                 2 = up (or ahead) */

    /* cwp local assignments */
    float d1; /* sample spacing for non-seismic data */
    float f1; /* first sample location for non-seismic data */
    float d2; /* sample spacing between traces */
    float f2; /* first trace location */
    float ungpow; /* negative of power used for dynamic range compression */
    float unscale; /* reciprocal of scaling factor to normalize range */
    int ntr;  /* number of traces */
    short mark; /* mark selected traces */
    short shortpad; /* alignment padding */

    short unass[14]; /* unassigned--NOTE: last entry causes
                     a break in the word alignment, if we REALLY
                     want to maintain 240 bytes, the following
                     entry should be an odd number of short/UINT2
                     OR do the insertion above the "mark" keyword
                     entry */
};

// 卷头400个字节
union REEL
{
    Struct_reelb400 reelstruct;
    short int binary[200];
};

// IBM float
union IBMFLOAT4
{
    float a;
    int b;
    char ch[4];
};

// IEEE float
union IEEEFLOAT4
{
    float a;
};

// INT4
union INT4
{
    int a;
    char ch[4];
};

// INT2
union INT2
{
    short a;
    char ch[2];
};

union Head240
{
    Struct_Su240 headstruct;
    int h4[60];
    short h2[120];
};

// Sgy数据中一道的相关信息：包含240字节的道头和数据
struct Trace{
    Head240 head;
    float *data;
};

/* 字节交换(swapbytes)函数 */
// swap_short_2
void swap(short *tni2);
// swap u_short_2
void swap(unsigned short *tni2);
// swap_int_4
void swap(int *tni4);
// swap_u_int_4
void swap(unsigned int *tni4);
// swap_long_4
void swap(long *tni4);
// swap_u_long_4
void swap(unsigned long *tni4);
// swap_float_4
void swap(float *tnf4);
// swap_double_8
void swap(double *tndd8);
// IBMSwarp
float IBMF4Swap(float x);


/* 读取Sgy文件的相关信息*/
bool InfoOfSgy(const char *FileName, REEL reel, unsigned short *TraceNum, unsigned short *SampleNum,
               unsigned short *SampleInt, short *DFormat, bool *BReel, bool *BIBM);

/* 读Sgy中的数据*/
bool ReadSgyData(const char *FileName, Trace *trace, REEL reel,
                 unsigned short *SampleNum, short *DFormat,
                 bool *BReel, bool *BIBM, const Partition &pt);

/* 将数据写到Sgy文件中:写成微机格式，IEEE的浮点类型 */
bool WriteSgy(const char * const FileName, unsigned char *f3200, Trace *trace, unsigned short TraceNum, unsigned short SampleNum,
              unsigned short SampleInt, const Partition &pt, const AFDP3D &Pa, MPI_Offset filesize, usht tag);


#endif //__RWSGY_H__

