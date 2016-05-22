/**************************************************************************
* 三维常密度声波方程有限差分正演

* 使用非分裂完全匹配曾（NPML）技术处理吸收边界

* 使用2阶位移运动方程

* 正演使用空间8阶时间2阶精度的交错网格有限差分技术

* 为了降低内存的存储量，使用混合反演算法，在时间域做正演，频率域求取梯度和步长

* 作者：高照奇

* 日期：2015-12-22
**************************************************************************/
#ifndef		_HybridFullWaveformInversion_Acoustic_3D_
#define		_HybridFullWaveformInversion_Acoustic_3D_

#include    "DataTran.h"
#include    "Partition.h"

#include	"RWsgy.h"
#include	"fftw3.h"

#define     TOP     0
#define     LEFT    1
#define     BOTTOM  2
#define     RIGHT   3
#define     FRONT   4
#define     BACK    5

#define     ADD         1
#define     SUBSTRACT   0

//#define		uint			unsigned int
//#define		usht			unsigned short
//#define		PI				3.141592653589793f

//// 空间8阶交错网格有限差分系数
//#define		C1_4			1.196289062506247f
//#define		C2_4			-0.079752604166839f
//#define		C3_4			0.009570312500021f
//#define		C4_4			-0.000697544642859f

//#define		PowerGrad		2.0f
//#define		StartFreq		2.0f * PI * (0.1f)
//#define		InterFreq		2.0f * PI * (0.4f)

//// 复数定义
//struct Complex{
//    float x;
//    float y;
//};

//// 震源位置
//typedef struct{
//    uint Sx;
//    uint Sy;
//    uint Sz;
//} SL;

//// 检波器位置
//typedef struct{
//    uint Rx;
//    uint Ry;
//    uint Rz;
//} RL;

//// 每一炮的信息
//struct Shot
//{
//    SL s;				// 震源位置
//    RL *re;				// 检波器位置数组
//    float *wf_t;		// 观测波场
//    float *wf_c;		// 正演模拟波场
//    float *wf_r;		// 残差波场
//    uint rn;			// 检波器个数
//};

//// 有限差分正演参数
//struct AFDP3D			// acoustic forward modeling parameters in 3D
//{
//    uint Nx;			// x方向网格数
//    uint Ny;			// y方向网格数
//    uint Nz;			// z方向网格数
//    uint Nt;			// 正演的时间步数
//    float dx;			// x方向网格间距
//    float dy;			// y方向网格间距
//    float dz;			// z方向网格间距
//    float dt;			// 时间采样间隔
//    float f0;			// ricker子波的主频
//    uint PMLx;			// x方向PML边界的网格数
//    uint PMLy;			// y方向PML边界的网格数
//    uint PMLz;			// z方向PML边界的网格数
//};

// 反演过程中使用的全局变量
struct CPUVs
{
    // 反演中使用的子波
    float *Wavedata;

    // NPML相关参数
    float *h_dx;		// x方向的NPML参数
    float *h_dy;		// y方向的NPML参数
    float *h_dz;		// z方向的NPML参数
    float *h_Bx;		// x方向的NPML参数
    float *h_By;		// y方向的NPML参数
    float *h_Bz;		// z方向的NPML参数
    float *h_PHIx_V_x;	// 波场变量V在x方向的卷积系数
    float *h_PHIy_K_y;	// 波场变量K在y方向的卷积系数
    float *h_PHIz_W_z;	// 波场变量W在z方向的卷积系数
    float *h_PHIx_U_x;	// 波场变量U在x方向的卷积系数
    float *h_PHIy_U_y;	// 波场变量U在y方向的卷积系数
    float *h_PHIz_U_z;	// 波场变量U在z方向的卷积系数
    float *h_PHIx_V_x_r;	// 残差反传时波场变量V在x方向的卷积系数
    float *h_PHIy_K_y_r;	// 残差反传时波场变量K在y方向的卷积系数
    float *h_PHIz_W_z_r;	// 残差反传时波场变量W在z方向的卷积系数
    float *h_PHIx_U_x_r;	// 残差反传时波场变量U在x方向的卷积系数
    float *h_PHIy_U_y_r;	// 残差反传时波场变量U在y方向的卷积系数
    float *h_PHIz_U_z_r;	// 残差反传时波场变量U在z方向的卷积系数

    // 波场
    float *h_U_past;	// (n-1)时刻的U波场快照
    float *h_U_now;		// (n+1)时刻的U波场快照
    float *h_U_next;	// (n)时刻的U波场快照
    float *h_V;			// (n)时刻的V波场快照
    float *h_K;			// (n)时刻的K波场快照
    float *h_W;			// (n)时刻的W波场快照
    float *h_U_past_r;		// (n-1)时刻的U波场快照(残差反传)
    float *h_U_now_r;		// (n)时刻的U波场快照(残差反传)
    float *h_U_next_r;		// (n+1)时刻的U波场快照(残差反传)
    float *h_V_r;			// (n)时刻的V波场快照(残差反传)
    float *h_K_r;			// (n)时刻的K波场快照(残差反传)
    float *h_W_r;			// (n)时刻的W波场快照(残差反传)

    float *h_Vp;		// 正演中使用的速度
    RL *h_re;			// 检波器位置

    // 波场记录的炮集
    float *h_TrueWF;	// 观测炮集
    float *h_CurrWF;	// 正演模拟炮集
    float *h_ResWF;		// 残差炮集

    // 频率域波场
    Complex *h_U_freq;		// 正传波场的频率域波场
    Complex *h_U_freq_r;	// 反传波场的频率域波场

    float *h_Grad;		// 梯度
    float *h_Obj;		// 目标函数值

    // 计算迭代的步长(目前暂定在频率域求取)
    float *h_TrailWF;
    Complex *h_ResTrail_freq;
    Complex *h_ResCurr_freq;
    float *h_Temp;
    float *h_SumFenzi;
    float *h_SumFenmu;

    // 自相关(拟Hessian)
    float *h_Grad_fz;
    float *h_Grad_fm;

    // 求取残差反传震源做FFTW的参数
    fftw_complex *fft_in;
    fftw_complex *fft_out;
};

//// 反演中要使用的参数
//struct IP{
//    uint IterN;			// 迭代的步数
//    uint ShotN;			// 反演中使用的炮数
//    Shot *St;			// 每一炮的相关信息
//    float *TrueVp;		// 真实的速度模型
//    float *CurrVp;		// 当前迭代的速度模型
//    float *GradVp;		// 目标函数对速度参数的梯度
//    float Alpha;		// 步长
//    float *ObjIter;		// 每一次迭代的目标函数值
//    uint FreqN;			// 反演中使用的频率个数
//};

/*------------------------ 反演程序中的函数 --------------------------*/
// Ricker子波
float Ricker(float f0,
             float t);

// 从Sgy文件中读取数据
void ReadData(const char *FileName,
              float *Data,
              const Partition &pt,
              usht flag);

void write_sgs_t_Data(const char * const FileName,
                      usht SampleNum,
                      usht TraceNum,
                      usht SampleInt,
                      float *data,
                      const Partition &pt,
                      const AFDP3D &Pa,
                      usht flag);

// 将数据写入Sgy文件
void WriteData(const char * const FileName,
            usht SampleNum, //nnz 234
            usht TraceNum, //nnx 610
            usht SampleInt,
            float *data,
            usht flag,
            const Partition& pt,
            const AFDP3D Pa,
            usht tag);

// 在内存上为变量申请空间
void MallocVariables(AFDP3D Pa,
                     IP *ip,
                     CPUVs *plan, const Partition &pt);

// 生成NPML系数
void GenerateNPML(AFDP3D Pa,
                  CPUVs *plan, const Partition &pt);

// 找出数组中的最大值
float MaxValue(float *Array,
               uint len);

// 正演加震源
void AddSource(AFDP3D Pa,
               float *h_U,
               SL s,
               float Wave,
               float *h_Vp);

// 一步更新波场U的卷积项
void StepPHIU(AFDP3D Pa,
              float *h_U,
              float *h_PHIx_U_x,
              float *h_PHIy_U_y,
              float *h_PHIz_U_z,
              float *h_Bx,
              float *h_By,
              float *h_Bz, const Partition &pt, int it);

// 一步更新波场V,K,W的卷积项
void StepPHIVKW(AFDP3D Pa,
                float *h_V,
                float *h_K,
                float *h_W,
                float *h_PHIx_V_x,
                float *h_PHIy_K_y,
                float *h_PHIz_W_z,
                float *h_Bx,
                float *h_By,
                float *h_Bz, const Partition &pt, int it);

// 一步更新波场U
void StepU(AFDP3D Pa,
           float *h_U_next,
           float *h_U_now,
           float *h_U_past,
           float *h_V,
           float *h_K,
           float *h_W,
           float *h_PHIx_V_x,
           float *h_PHIy_K_y,
           float *h_PHIz_W_z,
           float *h_Bx,
           float *h_By,
           float *h_Bz,
           float *h_Vp, const Partition &pt, int it);

// 一步更新波场V,K,W
void StepVKW(AFDP3D Pa,
             float *h_U,
             float *h_V,
             float *h_K,
             float *h_W,
             float *h_PHIx_U_x,
             float *h_PHIy_U_y,
             float *h_PHIz_U_z,
             float *h_Bx,
             float *h_By,
             float *h_Bz, const Partition &pt, int it);

// 一步记录检波器位置的波场值
void StepShotGather(AFDP3D Pa,
                    float *h_U,
                    float *h_SG,
                    RL *h_re,
                    uint nstep, const Partition &pt);

// 求取残差波场
void StepResidual(AFDP3D Pa,
                  float *h_TrueWF,
                  float *h_CurrWF,
                  float *h_ResWF,
                  uint Rn);

// 求取残差反传波场震源
void StepResidualConj(AFDP3D Pa,
                   CPUVs *plan,
                   uint Rn, const Partition &pt);

// 残差反传加震源
void AddResidual(AFDP3D Pa,
                 float *h_ResWF,
                 float *h_U_r,
                 RL *h_re,
                 uint nstep,
                 uint Rn,
                 float *h_Vp, const Partition &pt, int it);

// 梯度后处理
void PostProcessGrad(AFDP3D Pa,
                     float *h_Grad,
                     float *h_Vp, const Partition &pt);

// 更新内部网格的速度值
void UpdateVpInnerGrid(AFDP3D Pa,
                       float *h_Vp,
                       float *h_Grad,
                       float e, const Partition &pt);

// 更新PML边界内网格的速度值
void UpdateVpPML(AFDP3D Pa,
                 float *h_Vp,
                 float *h_Grad,
                 float e, const Partition &pt);

// 求取复数向量中每一维的共轭
void ComplexConjugate(Complex *h_data1,
                      Complex *h_data2,
                      uint len);

// 求取复数向量中每一维的模的平方
void CalSquareAbs(float *h_abs,
                  Complex *h_data,
                  uint len);


// 复数向量乘以常数
void ComplexMulConst(Complex *h_data1,
                     Complex *h_data2,
                     Complex a,
                     uint len);

// 一步在频率域计算梯度
void StepCalGradFreq(AFDP3D Pa,
                     float *h_Grad,
                     Complex *h_WF_f,
                     Complex *h_WF_b, const Partition &pt);

// 一步计算频率域的波场值
void StepCalWFFreq(AFDP3D Pa,
                   float *h_U,
                   float freq,
                   uint nt,
                   Complex *h_U_freq, const Partition &pt);

// 求取观测波场
void CalTrueWF(AFDP3D Pa,
               IP *ip,
               CPUVs *plan,
               float *sgs_t, const Partition &pt, const MPI_Comm &mycomm);

// 求取梯度
void CalGrad(AFDP3D Pa,
             IP *ip,
             CPUVs *plan,
             float *sgs_t,
             float *sgs_c,
             float *sgs_r,
             uint It,
             const Partition& pt,
             const MPI_Comm &mycomm);

// 计算频率域炮集
void CalSGFreq(AFDP3D Pa,
               uint Nt,
               uint Rn,
               float Freq,
               float *h_SG,
               Complex *h_SG_freq);

// 计算两个复向量点乘的实部
void ComplexDotMulReal(float *h_data,
                       Complex *h_data_f1,
                       Complex *h_data_f2,
                       uint len);


// 计算迭代步长
void CalStepLength(AFDP3D Pa,
                   IP *ip,
                   CPUVs *plan,
                   float *sgs_t,
                   float *sgs_c,
                   float e, const Partition &pt, const MPI_Comm &mycomm);

// 下次迭代前的预处理
void PreProcess(AFDP3D Pa,
                IP *ip,
                CPUVs *plan, const Partition &pt, const MPI_Comm &mycomm);






#endif
