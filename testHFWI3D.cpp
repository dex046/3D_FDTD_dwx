/**************************************************************************
* 三维常密度声波方程有限差分正演

* 使用非分裂完全匹配曾（NPML）技术处理吸收边界

* 使用2阶位移运动方程

* 正演使用空间8阶时间2阶精度的交错网格有限差分技术

* 为了降低内存的存储量，使用混合反演算法，在时间域做正演，频率域求取梯度和步长

* 作者：高照奇

* 日期：2015-12-23
**************************************************************************/

#include	"testHFWI3D.h"

using namespace std;

/*------------------------------------------------------------------------
function: Ricker

Information:
返回t时刻Ricker子波的值

Inputs:
f0: Ricker子波的主频
t: 时刻

Outputs:
t时刻Ricker子波的值
------------------------------------------------------------------------*/
float Ricker(float f0,
             float t)
{
    float RV = 0.0f, temp = 0.0f;

    temp = PI * f0 * (t - 1.0f / f0);
    RV = (1 - 2 * temp * temp) * exp(-1.0f * temp * temp);

    return RV;
}

/*------------------------------------------------------------------------
function: ReadData

Information:
read data from Sgy file

Inputs:
FileName: name of the Sgy file
Data: output data
------------------------------------------------------------------------*/
void ReadData(char FileName[],
            float *Data,
            const Partition &pt,
            usht flag)//flag ?
{
    FILE *f1;
    unsigned char f3200[3200];
    //Struct_reelb400 FileHeader;
    REEL *Reel;
    unsigned short *TraceNum, *SampleNum, *SampleInt;
    short *DFormat;
    bool BReel, BIBM;
    uint n,m;

    int rank = pt.getrank();

    TraceNum = new unsigned short[1];
    memset((void *)TraceNum, 0, sizeof(unsigned short));
    SampleNum = new unsigned short[1];
    memset((void *)SampleNum, 0, sizeof(unsigned short));
    SampleInt = new unsigned short[1];
    memset((void *)SampleInt, 0, sizeof(unsigned short));
    DFormat = new short[1];
    memset((void *)DFormat, 0, sizeof(short));
    Reel = new REEL[1];
    memset((void *)Reel, 0, sizeof(REEL));
//cout << "ww" << endl;
    f1 = fopen(FileName, "rb");
    if (f1 == NULL)
    {
        cout << "File Open error!" << "\n" << endl;
    }
    if(rank == ROOT_ID)
    {
        // read the first 3200 bytes
        fread(f3200, 3200, 1, f1);
    }

    // read the information of Sgy file
    InfoOfSgy(FileName, *Reel, TraceNum,
        SampleNum, SampleInt, DFormat,
        &BReel, &BIBM);
//coiut << 1231421 << endl;
//    cout << BReel << endl;
//    cout << *TraceNum << endl;
//    cout << *SampleNum << endl;
//    cout << *SampleInt << endl;
//    cout << *DFormat << endl;
//    cout << BIBM << endl;

    int length_x = pt.getblockLength_x();
    int length_z = pt.getblockLength_z();

//cout << "rank" << rank << "rank" << length_x << " " << length_z << " " << endl;
    Trace *trace;
    trace = new Trace[length_x];
    memset((void *)trace, 0, sizeof(Trace) * length_x);
    for (n = 0; n < length_x; n++)
    {
        trace[n].data = new float[length_z];
        memset((void *)trace[n].data, 0, sizeof(float) * length_z);
    }

    // read the data
    ReadSgyData(FileName, trace, *Reel, SampleNum,
        DFormat, &BReel, &BIBM, pt);


//    if(rank == 2)
//    {
//        for(n = 0; n < length_x; ++n)
//        {
//            for(m = 0; m < length_z; ++m)
//            {
//                cout << trace[n].data[m];
//            }
//        }
//    }
    // write the trace data to Data
    for(n = 0; n < length_x; ++n)
    {
        for(m = 0; m < length_z; ++m)
        {
            if(flag == 0)
            {
                Data[m * length_x + n] = trace[n].data[m];
                //cout << trace[n].data[m];
            }
            else
            {
                Data[n * length_z + m] = trace[n].data[m];
                //cout << trace[n].data[m];
            }
        }
    }

    // free memory
    for (n = 0; n < length_x; n++)
    {
        delete []trace[n].data;
    }
    delete []trace;
    delete []Reel;
    delete []DFormat;
    delete []SampleInt;
    delete []TraceNum;
    delete []SampleNum;

}
void write_sgs_t_Data(const char * const FileName, usht SampleNum, usht TraceNum, usht SampleInt, float *data, const Partition &pt, const AFDP3D &Pa, usht flag)
{
    unsigned char f3200[3200];
    memset((void*)f3200, 0, 3200);
    uint m, n;

    uint RL_num = pt.getRL_num();
    Trace *trace;
    trace = new Trace[RL_num];
    memset((void *)trace, 0, sizeof(Trace) * RL_num);

    int begin_num = pt.getRL_beginnum();
    int end_num = pt.getRL_endnum();

    int rank = pt.getrank();

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, FileName, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_Offset offset = 0, filesize = 0;
    MPI_Status status;


    filesize = 3600 + (SampleNum * sizeof(float) + 240) * TraceNum;

    MPI_File_set_size(fh, filesize);

    for (int n = 0; n < RL_num; n++)
    {
        trace[n].data = new float[SampleNum];
        memset((void *)trace[n].data, 0, sizeof(float) * SampleNum);
        for (m = 0; m < SampleNum; m++)
        {
            if (flag == 0)
            {
                trace[n].data[m] = *(data + m * RL_num + n);
            }
            else
            {
                trace[n].data[m] = *(data + n * SampleNum + m);
            }
        }
    }

    // trace header
    for (n = 0; n < RL_num; n++)
    {
        for (int m = 0; m < 60; m++)
        {
            trace[n].head.h4[m] = 0;
        }
    }

    if(rank == ROOT_ID)
    {
        // reel header
        for (n = 0; n < 3200; n++)
        {
            f3200[n] = 1;
        }
    }

    //cout << rank << endl;

    if(rank == ROOT_ID)
    {
        MPI_File_write_at(fh, offset, &f3200[0], 3200, MPI_BYTE, &status);

        // 写卷头中400个字节
        REEL reel;
        reel.reelstruct.hns = SampleNum;
        reel.reelstruct.hdt = SampleInt;
        reel.reelstruct.format = 5; // IEEE float
        reel.reelstruct.mfeet = 1;
        //fwrite(&reel.reelstruct, 400, 1, fp);
        MPI_File_write_at(fh, offset + 3200, &reel.reelstruct, 400, MPI_BYTE, &status);
    }

    offset += 3600 + begin_num * (240 + SampleNum * sizeof(float));

    for (int i = 0; i < RL_num; i++)
    {
        // 写道头
        trace[i].head.headstruct.cdp = i + begin_num;//
        trace[i].head.headstruct.ns = SampleNum;
        trace[i].head.headstruct.dt = SampleInt;
        trace[i].head.headstruct.sx = 100000;
        trace[i].head.headstruct.sy = 1000000 + (i + begin_num) * 4;//

        MPI_File_write_at(fh, offset, &trace[i].head.headstruct, 240, MPI_BYTE, &status);

        offset += 240;
        MPI_File_write_at(fh, offset, &trace[i].data[0], SampleNum, MPI_FLOAT, &status);
        offset += SampleNum * sizeof(float);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&fh);
    MPI_Barrier(MPI_COMM_WORLD);

    // free memory
    for (n = 0; n < RL_num; n++)
    {
        delete []trace[n].data;
    }
    delete []trace;
}

/*------------------------------------------------------------------------
function: WriteData

Information:
write data to Sgy file

Inputs:
FileName: name of the Sgy file
SampleNum: number of samples in each trace
TraceNum: number of traces
SampleInt: sample interval
data: output data
------------------------------------------------------------------------*/
void WriteData(const char * const FileName,
            usht SampleNum, //nnz 234
            usht TraceNum, //nnx 610
            usht SampleInt,
            float *data,
            usht flag,
            const Partition& pt,
            const AFDP3D Pa,
            usht tag)
{
    unsigned char f3200[3200];
    memset((void*)f3200, 0, 3200);
    uint m, n;

    MPI_Offset filesize = (SampleNum * sizeof(float) + 240) * TraceNum + 3600;

    int indexmin_x = pt.getindexmin_x();
    int indexmin_z = pt.getindexmin_z();

    uint block_x = 0;
    uint block_z = 0;

    int rank = pt.getrank();

    if(tag == WRITE_INTER)
    {
        block_x = pt.getinteriorLength_x();
        block_z = pt.getinteriorLength_z();

//        if(!block_x || !block_z)
//            return;
    }
    else if(tag == WRITE_ALL)
    {
        block_x = pt.getblockLength_x();
        block_z = pt.getblockLength_z();
    }
    else
    {

    }

    Trace *trace;
    trace = new Trace[block_x];
    memset((void *)trace, 0, sizeof(Trace) * block_x);

    for (int n = 0; n < block_x; n++)
    {
        trace[n].data = new float[block_z];
        memset((void *)trace[n].data, 0, sizeof(float) * block_z);
        for (m = 0; m < block_z; m++)
        {
            if (flag == 0)
            {
                trace[n].data[m] = *(data + m * block_x + n);
            }
            else
            {
                trace[n].data[m] = *(data + n * block_z + m);
            }
        }
    }

    // trace header
    for (n = 0; n < block_x; n++)
    {
        for (int m = 0; m < 60; m++)
        {
            trace[n].head.h4[m] = 0;
        }
    }

    if(rank == ROOT_ID)
    {
        // reel header
        for (n = 0; n < 3200; n++)
        {
            f3200[n] = 1;
        }
    }

    WriteSgy(FileName, &f3200[0], trace,
        TraceNum, SampleNum, SampleInt, pt, Pa, filesize, tag);///hai

    // free memory
    for (n = 0; n < block_x; n++)
    {
        delete []trace[n].data;
    }
    delete []trace;

}

/*------------------------------------------------------------------------
function: MallocVariables

Information:
在内存上为变量申请空间

Inputs:
Pa: 正演参数
ip: 反演参数
plan: 全局变量
------------------------------------------------------------------------*/
void MallocVariables(AFDP3D Pa,
                     IP *ip,
                     CPUVs *plan,
                     const Partition& pt)
{
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

//    if(pt.getrank() == 2)
//    {
//        cout << block_x << " " << block_y << " " << block_z << endl;
//    }

    uint RL_num = pt.getRL_num();

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_y = pt.getinteriorLength_y();
    uint interiorlength_z = pt.getinteriorLength_z();

    H_Border h_U = pt.geth_U();
    H_Border h_Vp = pt.geth_Vp();
    H_Border h_VW = pt.geth_VW();

    uint sum_h_Coord = pt.geth_Coord_length();





//    if(pt.getrank() == ROOT_ID)
//    {
//        cout << block_x << " " << block_y << " " << block_z << endl;
//        cout << interiorlength_x << " " << interiorlength_y << " " << interiorlength_z << endl;
//    }
    try
    {
        plan->h_dx = new float[block_x];
        plan->h_dy = new float[block_y];
        plan->h_dz = new float[block_z];
        plan->h_Bx = new float[block_x];
        plan->h_By = new float[block_y];
        plan->h_Bz = new float[block_z];

        plan->h_PHIx_V_x = new float[block_z * block_y * block_x];
        plan->h_PHIy_K_y = new float[block_z * block_y * block_x];
        plan->h_PHIz_W_z = new float[block_z * block_y * block_x];
        plan->h_PHIx_U_x = new float[block_z * block_y * block_x];
        plan->h_PHIy_U_y = new float[block_z * block_y * block_x];
        plan->h_PHIz_U_z = new float[block_z * block_y * block_x];

        plan->h_PHIx_V_x_r = new float[block_z * block_y * block_x];
        plan->h_PHIy_K_y_r = new float[block_z * block_y * block_x];
        plan->h_PHIz_W_z_r = new float[block_z * block_y * block_x];
        plan->h_PHIx_U_x_r = new float[block_z * block_y * block_x];
        plan->h_PHIy_U_y_r = new float[block_z * block_y * block_x];
        plan->h_PHIz_U_z_r = new float[block_z * block_y * block_x];

        plan->h_U_past = new float[block_z * block_y * block_x];
        plan->h_U_now = new float[h_U.length_z * h_U.length_y * h_U.length_x];
        plan->h_U_next = new float[block_z * block_y * block_x];
        plan->h_V = new float[h_VW.length_z * h_VW.length_y * h_VW.length_x];
        plan->h_K = new float[h_VW.length_z * h_VW.length_y * h_VW.length_x];
        plan->h_W = new float[h_VW.length_z * h_VW.length_y * h_VW.length_x];

        plan->h_U_past_r = new float[block_z * block_y * block_x];
        plan->h_U_now_r = new float[h_U.length_z * h_U.length_y * h_U.length_x];
        plan->h_U_next_r = new float[block_z * block_y * block_x];
        plan->h_V_r = new float[h_VW.length_z * h_VW.length_y * h_VW.length_x];
        plan->h_K_r = new float[h_VW.length_z * h_VW.length_y * h_VW.length_x];
        plan->h_W_r = new float[h_VW.length_z * h_VW.length_y * h_VW.length_x];

        plan->h_Vp = new float[h_Vp.length_z * h_Vp.length_y * h_Vp.length_x];
        plan->h_re = new RL[ip->St[0].rn];

        if(RL_num)
        {
            plan->h_TrueWF = new float[RL_num * Pa.Nt];
            plan->h_CurrWF = new float[RL_num * Pa.Nt];
            plan->h_ResWF = new float[RL_num * Pa.Nt];

            plan->h_Obj = new float[RL_num * Pa.Nt];

            plan->h_TrailWF = new float[RL_num * Pa.Nt];
            plan->h_ResTrail_freq = new Complex[RL_num];
            plan->h_ResCurr_freq = new Complex[RL_num];

            plan->h_Temp = new float[RL_num];
        }

        if(interiorlength_x && interiorlength_y && interiorlength_z)
        {
            plan->h_U_freq = new Complex[interiorlength_z * interiorlength_y * interiorlength_x * ip->FreqN];
            plan->h_U_freq_r = new Complex[interiorlength_z * interiorlength_y * interiorlength_x * ip->FreqN];

            plan->h_Grad = new float[interiorlength_z * interiorlength_y * interiorlength_x];

            plan->h_Grad_fz = new float[interiorlength_z * interiorlength_y * interiorlength_x];
            plan->h_Grad_fm = new float[interiorlength_z * interiorlength_y * interiorlength_x];
        }

        plan->h_SumFenzi = new float[1];
        plan->h_SumFenmu = new float[1];

        plan->fft_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Pa.Nt);
        plan->fft_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Pa.Nt);

//        plan->fft_in = new fftw_complex[Pa.Nt];
//        plan->fft_out = new fftw_complex[Pa.Nt];
    }
    catch(const std::bad_alloc& e)
    {
        cout << "Allocation failed: www   " << pt.getrank() << endl;
        cout << block_x << " " << block_y << " " << block_z << endl;
        cout << interiorlength_x << " " << interiorlength_y << " " << interiorlength_z << endl;
        cout << RL_num << endl;
        cout << h_VW.length_z * h_VW.length_y * h_VW.length_x << endl;
        cout << h_Vp.length_z * h_Vp.length_y * h_Vp.length_x << endl;
    }



    // 给全局变量申请内存空间
    plan->Wavedata = new float[Pa.Nt];





    // 初始化上述内存空间
    memset((void *)plan->Wavedata,		0,	sizeof(float) * Pa.Nt);

    memset((void *)plan->h_dx,			0,	sizeof(float) * block_x);
    memset((void *)plan->h_dy,			0,	sizeof(float) * block_y);
    memset((void *)plan->h_dz,			0,	sizeof(float) * block_z);
    memset((void *)plan->h_Bx,			0,	sizeof(float) * block_x);
    memset((void *)plan->h_By,			0,	sizeof(float) * block_y);
    memset((void *)plan->h_Bz,			0,	sizeof(float) * block_z);

    memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIy_U_y,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIy_K_y,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_y * block_x);

    memset((void *)plan->h_PHIx_U_x_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIy_U_y_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIz_U_z_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIx_V_x_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIy_K_y_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_PHIz_W_z_r,	0,	sizeof(float) * block_z * block_y * block_x);

    memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_U_now,		0,	sizeof(float) * h_U.length_z * h_U.length_y * h_U.length_x);
    memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
    memset((void *)plan->h_K,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
    memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);

    memset((void *)plan->h_U_past_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_U_now_r,		0,	sizeof(float) * h_U.length_z * h_U.length_y * h_U.length_x);
    memset((void *)plan->h_U_next_r,	0,	sizeof(float) * block_z * block_y * block_x);
    memset((void *)plan->h_V_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
    memset((void *)plan->h_K_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
    memset((void *)plan->h_W_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);

    memset((void *)plan->h_Vp,			0,	sizeof(float) * h_Vp.length_z * h_Vp.length_y * h_Vp.length_x);
    memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);

    if(RL_num)
    {
        memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_CurrWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_ResWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_TrailWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_ResTrail_freq,	0,	sizeof(Complex) * RL_num);
        memset((void *)plan->h_ResCurr_freq,	0,	sizeof(Complex) * RL_num);
        memset((void *)plan->h_Temp,		0,	sizeof(float) * RL_num);
        memset((void *)plan->h_Obj,			0,	sizeof(float) * Pa.Nt * RL_num);
    }

    if(interiorlength_x && interiorlength_y && interiorlength_z)
    {
        memset((void *)plan->h_U_freq,		0,	sizeof(Complex) * interiorlength_z * interiorlength_y * interiorlength_x * ip->FreqN);
        memset((void *)plan->h_U_freq_r,	0,	sizeof(Complex) * interiorlength_z * interiorlength_y * interiorlength_x * ip->FreqN);

        memset((void *)plan->h_Grad,		0,	sizeof(float) * interiorlength_z * interiorlength_y * interiorlength_x);

        memset((void *)plan->h_Grad_fz,		0,	sizeof(float) * interiorlength_z * interiorlength_y * interiorlength_x);
        memset((void *)plan->h_Grad_fm,		0,	sizeof(float) * interiorlength_z * interiorlength_y * interiorlength_x);
    }

    memset((void *)plan->h_SumFenzi,	0,	sizeof(float));
    memset((void *)plan->h_SumFenmu,	0,	sizeof(float));

    memset((void *)plan->fft_in, 0, sizeof(fftw_complex) * Pa.Nt);
    memset((void *)plan->fft_out, 0, sizeof(fftw_complex) * Pa.Nt);
}

/*------------------------------------------------------------------------
function: GenerateNPML

Information:
生成NPML系数

Inputs:
Pa: 正演参数
plan: 全局变量
------------------------------------------------------------------------*/
void GenerateNPML(AFDP3D Pa,
                  CPUVs *plan,
                  const Partition& pt)
{
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    uint block_z = pt.getblockLength_z();
    uint block_y = pt.getblockLength_y();
    uint block_x = pt.getblockLength_x();

    uint min_x = pt.getindexmin_x();
    uint min_y = pt.getindexmin_y();
    uint min_z = pt.getindexmin_z();

    // 对NPML参数进行赋值
    // z方向
    float d0 = -3.0f * 2500.0f * logf(1.0e-6f)
        / (2.0f * powf((float)(Pa.PMLz * Pa.dz), 3.0f));

    for (uint iz = 0; iz < block_z; iz++)
    {
        int actual_iz = iz + min_z;
        if (actual_iz < Pa.PMLz)
        {
            plan->h_dz[iz] = d0 * (Pa.PMLz - 1 - actual_iz + 0.5f)
                * (Pa.PMLz - 1 - actual_iz + 0.5f) * Pa.dz * Pa.dz;
        }
        else if (actual_iz > Pa.PMLz + Pa.Nz - 1)
        {
            plan->h_dz[iz] = d0 * (actual_iz - Pa.PMLz - Pa.Nz + 0.5f)
                * (actual_iz - Pa.PMLz - Pa.Nz + 0.5f) * Pa.dz * Pa.dz;
        }
        else
        {
            plan->h_dz[iz] = 0.0f;
        }
        plan->h_Bz[iz] = expf(-1.0f * plan->h_dz[iz] * Pa.dt);



    }

    // x方向
    d0 = -3.0f * 2500.0f * logf(1.0e-6f) /
        (2.0f * powf((float)(Pa.PMLx * Pa.dx), 3.0f));

    for (uint ix = 0; ix < block_x; ix++)
    {
        int actual_ix = ix + min_x;
        if (actual_ix < Pa.PMLx)
        {
            plan->h_dx[ix] = d0 * (Pa.PMLx - 1 - actual_ix + 0.5f) *
                (Pa.PMLx - 1 - actual_ix + 0.5f) * Pa.dx * Pa.dx;
        }
        else if (actual_ix > Pa.PMLx + Pa.Nx - 1)
        {
            plan->h_dx[ix] = d0 * (actual_ix - Pa.PMLx - Pa.Nx + 0.5f) *
                (actual_ix - Pa.PMLx - Pa.Nx + 0.5f) * Pa.dx * Pa.dx;
        }
        else
        {
            plan->h_dx[ix] = 0.0f;
        }
        plan->h_Bx[ix] = expf(-1.0f * plan->h_dx[ix] * Pa.dt);


    }

    // y方向
    d0 = -3.0f * 2500.0f * logf(1.0e-6f) /
        (2.0f * powf((float)(Pa.PMLy * Pa.dy), 3.0f));

    for (uint iy = 0; iy < block_y; iy++)
    {
        int actual_iy = iy + min_y;
        if (actual_iy < Pa.PMLy)
        {
            plan->h_dy[iy] = d0 * (Pa.PMLy - 1 - actual_iy + 0.5f) *
                (Pa.PMLy - 1 - actual_iy + 0.5f) * Pa.dy * Pa.dy;
        }
        else if (actual_iy > Pa.PMLy + Pa.Ny - 1)
        {
            plan->h_dy[iy] = d0 * (actual_iy - Pa.PMLy - Pa.Ny + 0.5f) *
                (actual_iy - Pa.PMLy - Pa.Ny + 0.5f) * Pa.dy * Pa.dy;
        }
        else
        {
            plan->h_dy[iy] = 0.0f;
        }
        plan->h_By[iy] = expf(-1.0f * plan->h_dy[iy] * Pa.dt);


    }
}


/*------------------------------------------------------------------------
function: MaxValue

Information:
找出数组中的最大值

Inputs:
Array: 数组
len：数组的维数
------------------------------------------------------------------------*/
float MaxValue(float *Array,
               uint len)
{
    float MV = fabs(Array[0]);

    for (uint n = 0; n < len; n++)
    {
        if (MV < fabs(Array[n]))
        {
            MV = fabs(Array[n]);
        }
    }

    return MV;
}

/*------------------------------------------------------------------------
function: AddSource

Information:
正演加震源

Inputs:
Pa：正演参数
h_U：当前时刻的波场
s: 震源位置
Wave: 当前时刻波场的值
h_Vp: 当前模型速度值的平方
------------------------------------------------------------------------*/
void AddSource(AFDP3D Pa,
               float *h_U,
               SL s,
               float Wave,
               float *h_Vp,
               const Partition& pt)
{
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;

    H_Border temph_Vp = pt.geth_Vp();
    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_y = pt.getindexmin_y();
    uint indexmin_z = pt.getindexmin_z();

    uint ix = s.Sx - indexmin_x;
    uint iy = s.Sy - indexmin_y;
    uint iz = s.Sz - indexmin_z;

    //uint id = iz * temph_Vp.length_x * temph_Vp.length_y + iy * temph_Vp.length_x + ix;

    //cout << h_U[iz * block_x * block_y + iy * block_x + ix] << endl;
    h_U[iz * block_x * block_y + iy * block_x + ix] += Wave * (Pa.dt * Pa.dt * h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder]);

    //cout << h_U[iz * block_x * block_y + iy * block_x + ix] << endl;
}

/*------------------------------------------------------------------------
function: StepPHIU

Information:
一步更新波场U的卷积项

Inputs:
Pa：正演参数
h_U：当前时刻的波场
h_PHIx_U_x: x方向卷积项
h_PHIy_U_y: y方向卷积项
h_PHIz_U_z: z方向卷积项
h_Bx: x方向NPML参数
h_By: y方向NPML参数
h_Bz：z方向NPML参数
------------------------------------------------------------------------*/
void StepPHIU(AFDP3D Pa,
              float *h_U,
              float *h_PHIx_U_x,
              float *h_PHIy_U_y,
              float *h_PHIz_U_z,
              float *h_Bx,
              float *h_By,
              float *h_Bz,
              const Partition& pt, int it)
{
    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnz = Pa.Nz + 2 * Pa.PMLz;

    H_Border temph_U = pt.geth_U();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_y = pt.gettotallength_y();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_y = pt.getindexmin_y();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_y = pt.getindexmax_y();
    int max_z = pt.getindexmax_z();

    float dUx = 0.0f;
    float dUy = 0.0f;
    float dUz = 0.0f;

    int gap_min_x = 0, gap_min_z = 0, gap_min_y = 0, gap_max_y = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_y < 4)
    {
        gap_min_y = 4 - min_y;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }

    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_y > totallength_y - 5)
    {
        gap_max_y = totallength_y - 5 - max_y;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

//    if(pt.getrank() == 1)
//    {
//        cout << "min_x=" << min_x << endl;
//        cout << temph_U.leftborder << temph_U.backborder << temph_U.topborder << endl;
//        cout << gap_min_z << " " << gap_min_y << " " << gap_min_x << endl;
//        cout << block_z + temph_U.topborder + gap_max_z << " "  << block_y + temph_U.backborder + gap_max_y << " " << block_x + temph_U.leftborder + gap_max_x << endl;
//    }
    // 空间8阶精度的交错网格有限差分
    for (uint iz = 0 + gap_min_z + temph_U.topborder; iz < block_z + temph_U.topborder + gap_max_z; iz++)
    {
        for (uint iy = 0 + gap_min_y + temph_U.backborder; iy < block_y + temph_U.backborder + gap_max_y; iy++)
        {
            for (uint ix = 0 + gap_min_x + temph_U.leftborder; ix < block_x + temph_U.leftborder + gap_max_x; ix++)
            {
                dUx = C1_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 1] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C2_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 2] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix - 1])
                    + C3_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 3] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix - 2])
                    + C4_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 4] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix - 3]);

                dUy = C1_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 1) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C2_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 2) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + (iy - 1) * temph_U.length_x + ix])
                    + C3_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 3) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + (iy - 2) * temph_U.length_x + ix])
                    + C4_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 4) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + (iy - 3) * temph_U.length_x + ix]);

                dUz = C1_4 * (h_U[(iz + 1) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C2_4 * (h_U[(iz + 2) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[(iz - 1) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C3_4 * (h_U[(iz + 3) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[(iz - 2) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C4_4 * (h_U[(iz + 4) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[(iz - 3) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix]);

                //cout << dUx << endl << dUy << endl << dUz << endl;
//                if(pt.getrank() == ROOT_ID)
//                {
//                    //cout << ix << " " << iy << " " << iz << endl;
//                    if(dUx || dUy || dUz)
//                    cout << dUx << " " << dUy << " " << dUz << endl;
//                }

                uint id = (iz - temph_U.topborder) * block_y * block_x + (iy - temph_U.backborder) * block_x + ix - temph_U.leftborder;
                h_PHIx_U_x[id] = h_Bx[ix - temph_U.leftborder] * (h_PHIx_U_x[id] + dUx / Pa.dx) - dUx / Pa.dx;
                h_PHIy_U_y[id] = h_By[iy - temph_U.backborder] * (h_PHIy_U_y[id] + dUy / Pa.dy) - dUy / Pa.dy;
                h_PHIz_U_z[id] = h_Bz[iz - temph_U.topborder] * (h_PHIz_U_z[id] + dUz / Pa.dz) - dUz / Pa.dz;
            }
        }
    }

//    MPI_Barrier(MPI_COMM_WORLD);
//    if(pt.getrank() == ROOT_ID)
//    {
//        cout << "it=" << it << endl;
//    }

}

/*------------------------------------------------------------------------
function: StepPHIVKW

Information:
一步更新波场V,K,W的卷积项

Inputs:
Pa：正演参数
h_V：当前时刻的波场V
h_K: 当前时刻的波场K
h_W：当前时刻的波场W
h_PHIx_V_x: x方向卷积项
h_PHIy_K_y: y方向卷积项
h_PHIz_W_z: z方向卷积项
h_Bx: x向NPML参数
h_By: y向NPML参数
h_Bz：z向NPML参数
------------------------------------------------------------------------*/
void StepPHIVKW(AFDP3D Pa,
                float *h_V,
                float *h_K,
                float *h_W,
                float *h_PHIx_V_x,
                float *h_PHIy_K_y,
                float *h_PHIz_W_z,
                float *h_Bx,
                float *h_By,
                float *h_Bz,
                const Partition& pt,
                int it)
{
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;

    H_Border h_VW = pt.geth_VW();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_y = pt.gettotallength_y();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_y = pt.getindexmin_y();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_y = pt.getindexmax_y();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_z = 0, gap_min_y = 0, gap_max_y = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_y < 4)
    {
        gap_min_y = 4 - min_y;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }

    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_y > totallength_y - 5)
    {
        gap_max_y = totallength_y - 5 - max_y;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

    float dVx = 0.0f;
    float dKy = 0.0f;
    float dWz = 0.0f;

    // 空间8阶精度的交错网格有限差分
    for (uint iz = 0 + gap_min_z + h_VW.topborder; iz < block_z + h_VW.topborder + gap_max_z; iz++)
    {
        for (uint iy = 0 + gap_min_y + h_VW.backborder; iy < block_y + h_VW.backborder + gap_max_y; iy++)
        {
            for (uint ix = 0 + gap_min_x + h_VW.leftborder; ix < block_x + h_VW.leftborder + gap_max_x; ix++)
            {
                dVx = C1_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 1])
                    + C2_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix + 1] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 2])
                    + C3_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix + 2] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 3])
                    + C4_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix + 3] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 4]);

                dKy = C1_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 1) * h_VW.length_x + ix])
                    + C2_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + (iy + 1) * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 2) * h_VW.length_x + ix])
                    + C3_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + (iy + 2) * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 3) * h_VW.length_x + ix])
                    + C4_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + (iy + 3) * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 4) * h_VW.length_x + ix]);

                dWz = C1_4 * (h_W[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 1) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix])
                    + C2_4 * (h_W[(iz + 1) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 2) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix])
                    + C3_4 * (h_W[(iz + 2) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 3) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix])
                    + C4_4 * (h_W[(iz + 3) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 4) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix]);

                uint id = (iz - h_VW.topborder) * block_y * block_x + (iy - h_VW.backborder) * block_x + ix - h_VW.leftborder;
                h_PHIx_V_x[id] = h_Bx[ix - h_VW.leftborder] * (h_PHIx_V_x[id] + dVx / Pa.dx) - dVx / Pa.dx;
                h_PHIy_K_y[id] = h_By[iy - h_VW.backborder] * (h_PHIy_K_y[id] + dKy / Pa.dy) - dKy / Pa.dy;
                h_PHIz_W_z[id] = h_Bz[iz - h_VW.topborder] * (h_PHIz_W_z[id] + dWz / Pa.dz) - dWz / Pa.dz;

            }
        }
    }

}

/*------------------------------------------------------------------------
function: StepU

Information:
一步更新波场U

Inputs:
Pa：正演参数
h_U_next： 下一时刻的波场U
h_U_now： 当前时刻的波场U
h_U_past： 前一时刻的波场U
h_W：当前时刻的波场W
h_K：当前时刻的波场K
h_V：当前时刻的波场V
h_PHIx_V_x: x向卷积项
h_PHIy_K_y: y向卷积项
h_PHIz_W_z: z向卷积项
h_Bx: x向NPML参数
h_By: y向NPML参数
h_Bz：z向NPML参数
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
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
           float *h_Vp,
           const Partition& pt,
           int it)
{
    H_Border h_VW = pt.geth_VW();
    H_Border temph_U = pt.geth_U();
    H_Border temph_Vp = pt.geth_Vp();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_y = pt.gettotallength_y();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_y = pt.getindexmin_y();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_y = pt.getindexmax_y();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_y = 0, gap_min_z = 0, gap_max_x = 0, gap_max_y = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_y < 4)
    {
        gap_min_y = 4 - min_y;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }

    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_y > totallength_y - 5)
    {
        gap_max_y = totallength_y - 5 - max_y;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    float dVx = 0.0f;
    float dKy = 0.0f;
    float dWz = 0.0f;

    // 时间2阶精度，空间8阶精度的交错网格有限差分
    for (uint iz = 0 + gap_min_z + h_VW.topborder; iz < block_z + h_VW.topborder + gap_max_z; iz++)
    {
        for (uint iy = 0 + gap_min_y + h_VW.backborder; iy < block_y + h_VW.backborder + gap_max_y; iy++)
        {
            for (uint ix = 0 + gap_min_x + h_VW.leftborder; ix < block_x + h_VW.leftborder + gap_max_x; ix++)
            {
                dVx = C1_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 1])
                    + C2_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix + 1] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 2])
                    + C3_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix + 2] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 3])
                    + C4_4 * (h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix + 3] - h_V[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix - 4]);

                dKy = C1_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 1) * h_VW.length_x + ix])
                    + C2_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + (iy + 1) * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 2) * h_VW.length_x + ix])
                    + C3_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + (iy + 2) * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 3) * h_VW.length_x + ix])
                    + C4_4 * (h_K[iz * h_VW.length_y * h_VW.length_x + (iy + 3) * h_VW.length_x + ix] - h_K[iz * h_VW.length_y * h_VW.length_x + (iy - 4) * h_VW.length_x + ix]);

                dWz = C1_4 * (h_W[iz * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 1) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix])
                    + C2_4 * (h_W[(iz + 1) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 2) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix])
                    + C3_4 * (h_W[(iz + 2) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 3) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix])
                    + C4_4 * (h_W[(iz + 3) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix] - h_W[(iz - 4) * h_VW.length_y * h_VW.length_x + iy * h_VW.length_x + ix]);

                uint id = (iz - h_VW.topborder) * block_y * block_x + (iy - h_VW.backborder) * block_x + ix - h_VW.leftborder;

                h_U_next[id] = (Pa.dt * Pa.dt * h_Vp[(iz - h_VW.topborder + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy - h_VW.backborder + temph_Vp.backborder) * temph_Vp.length_x + ix - h_VW.leftborder + temph_Vp.leftborder])
                    * (dVx / Pa.dx + h_PHIx_V_x[id] +
                    dKy / Pa.dy + h_PHIy_K_y[id] +
                    dWz / Pa.dz + h_PHIz_W_z[id]) -
                    h_U_past[id] + 2.0f * h_U_now[(iz - h_VW.topborder + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy - h_VW.backborder + temph_U.backborder) * temph_U.length_x + ix - h_VW.leftborder + temph_U.leftborder];

            }
        }
    }
}

/*------------------------------------------------------------------------
function: StepVKW

Information:
一步更新波场V,K和W

Inputs:
Pa：正演参数
h_U： 当前时刻的波场U
h_W：当前时刻的波场W
h_K: 当前时刻的波场K
h_V：当前时刻的波场V
h_PHIx_U_x: x向卷积项
h_PHIy_U_y: y向卷积项
h_PHIz_U_z: x向卷积项
h_Bx: x向NPML参数
h_By: y向NPML参数
h_Bz：z向NPML参数
------------------------------------------------------------------------*/
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
             float *h_Bz,
             const Partition& pt,
             int it)
{
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;

    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_y = pt.gettotallength_y();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_y = pt.getindexmin_y();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_y = pt.getindexmax_y();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_y = 0, gap_min_z = 0, gap_max_x = 0, gap_max_y = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_y < 4)
    {
        gap_min_y = 4 - min_y;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_y > totallength_y - 5)
    {
        gap_max_y = totallength_y - 5 - max_y;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

    float dUx = 0.0f;
    float dUy = 0.0f;
    float dUz = 0.0f;

    // 空间8阶精度的交错网格有限差分
    for (uint iz = 0 + gap_min_z + temph_U.topborder; iz < block_z + temph_U.topborder + gap_max_z; iz++)
    {
        for (uint iy = 0 + gap_min_y + temph_U.backborder; iy < block_y + temph_U.backborder + gap_max_y; iy++)
        {
            for (uint ix = 0 + gap_min_x + temph_U.leftborder; ix < block_x + temph_U.leftborder + gap_max_x; ix++)
            {
                dUx = C1_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 1] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C2_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 2] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix - 1])
                    + C3_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 3] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix - 2])
                    + C4_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix + 4] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix - 3]);

                dUy = C1_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 1) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C2_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 2) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + (iy - 1) * temph_U.length_x + ix])
                    + C3_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 3) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + (iy - 2) * temph_U.length_x + ix])
                    + C4_4 * (h_U[iz * temph_U.length_y * temph_U.length_x + (iy + 4) * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + (iy - 3) * temph_U.length_x + ix]);

                dUz = C1_4 * (h_U[(iz + 1) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[iz * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C2_4 * (h_U[(iz + 2) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[(iz - 1) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C3_4 * (h_U[(iz + 3) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[(iz - 2) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix])
                    + C4_4 * (h_U[(iz + 4) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix] - h_U[(iz - 3) * temph_U.length_y * temph_U.length_x + iy * temph_U.length_x + ix]);

                uint id0 = (iz - temph_U.topborder + h_VW.topborder) * h_VW.length_y * h_VW.length_x + (iy - temph_U.backborder + h_VW.backborder) * h_VW.length_x + ix - temph_U.leftborder + h_VW.leftborder;
                uint id1 = (iz - temph_U.topborder) * block_y * block_x + (iy - temph_U.backborder) * block_x + ix - temph_U.leftborder;

                h_V[id0] = dUx / Pa.dx + h_PHIx_U_x[id1];
                h_K[id0] = dUy / Pa.dy + h_PHIy_U_y[id1];
                h_W[id0] = dUz / Pa.dz + h_PHIz_U_z[id1];

            }
        }
    }

}

/*------------------------------------------------------------------------
function: StepShotGather

Information:
一步记录检波器位置的波场值

Inputs:
Pa：正演参数
h_U：当前时刻的波场U
h_SG: 炮集
h_re：检波器位置
nstep：循环到第几个时刻
rn：检波器个数
------------------------------------------------------------------------*/
void StepShotGather(AFDP3D Pa,
                    float *h_U,
                    float *h_SG,
                    RL *h_re,
                    uint nstep,
                    const Partition& pt)
{
    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint id = 0;

    H_Border temph_U = pt.geth_U();

    uint rn = pt.getRL_num();
    vector<tuple<uint, uint, uint>> vec = pt.getRL();
    vector<tuple<uint, uint, uint>>::iterator begin = vec.begin();
    tuple<uint, uint, uint> temp;
    for (uint ir = 0; ir < rn; ir++)
    {
        temp = *(begin + ir);
        id = (get<2>(temp) + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (get<1>(temp) + temph_U.backborder) * temph_U.length_x + get<0>(temp) + temph_U.leftborder;
        //id = h_re[ir].Rz * nny * nnx + h_re[ir].Ry * nnx + h_re[ir].Rx;
        h_SG[ir * Pa.Nt + nstep] = h_U[id];
    }
}

/*------------------------------------------------------------------------
function: StepResidual

Information:
求取残差波场

Inputs:
Pa：正演参数
h_TrueWF：观测波场
h_CurrWF: 当前模型对应的正演波场
h_re：检波器位置
nstep：循环到第几个时刻
rn：检波器个数
------------------------------------------------------------------------*/
void StepResidual(AFDP3D Pa,
                  float *h_TrueWF,
                  float *h_CurrWF,
                  float *h_ResWF,
                  uint Rn)
{
    uint id = 0;
    for (uint n = 0; n < Rn; n++)
    {
        for (uint m = 0; m < Pa.Nt; m++)
        {
            id = n * Pa.Nt + m;
            h_ResWF[id] = h_TrueWF[id] - h_CurrWF[id];
        }
    }
}

/*------------------------------------------------------------------------
function: StepResidualConj

Information:
求取残差反传波场震源

Inputs:
Pa：正演参数
h_TrueWF：观测波场
h_CurrWF: 当前模型对应的正演波场
h_re：检波器位置
nstep：循环到第几个时刻
rn：检波器个数
------------------------------------------------------------------------*/
void StepResidualConj(AFDP3D Pa,
                   CPUVs *plan,
                   uint Rn,
                   const Partition& pt)
{
    uint id = 0;

    ofstream fout;
    if(pt.getrank() == 0)
    {
        fout.open("h_fft_in00.txt", ios_base::out);
    }
    else if(pt.getrank() == 1)
    {
        fout.open("h_fft_in01.txt", ios_base::out);
    }

    ofstream fout1;
    if(pt.getrank() == 0)
    {
        fout1.open("h_fft_in10.txt", ios_base::out);
    }
    else if(pt.getrank() == 1)
    {
        fout1.open("h_fft_in11.txt", ios_base::out);
    }



    ofstream fout2;
    if(pt.getrank() == 0)
    {
        fout2.open("h_fft_out00.txt", ios_base::out);
    }
    else if(pt.getrank() == 1)
    {
        fout2.open("h_fft_out01.txt", ios_base::out);
    }

    ofstream fout3;
    if(pt.getrank() == 0)
    {
        fout3.open("h_fft_out10.txt", ios_base::out);
    }
    else if(pt.getrank() == 1)
    {
        fout3.open("h_fft_out11.txt", ios_base::out);
    }



    ofstream fout4;
    if(pt.getrank() == 0)
    {
        fout4.open("h_fft_in20.txt", ios_base::out);
    }
    else if(pt.getrank() == 1)
    {
        fout4.open("h_fft_in21.txt", ios_base::out);
    }

    ofstream fout5;
    if(pt.getrank() == 0)
    {
        fout5.open("h_fft_out20.txt", ios_base::out);
    }
    else if(pt.getrank() == 1)
    {
        fout5.open("h_fft_out21.txt", ios_base::out);
    }


//    ofstream fout1;
//    if(pt.getrank() == 0)
//    {
//        fout1.open("h_fft_in00.txt", ios_base::out | ios_base::app);
//    }
//    else
//    {
//        fout1.open("h_fft_in01.txt", ios_base::out | ios_base::app);
//    }




    // 在这个子函数中，我们将使用FFTW来进行运算
    for (uint n = 0; n < Rn; n++)
    {


        memset((void *)plan->fft_in, 0, sizeof(fftw_complex) * Pa.Nt);
        memset((void *)plan->fft_out, 0, sizeof(fftw_complex) * Pa.Nt);

        for (uint m = 0; m < Pa.Nt; m++)
        {
            id = n * Pa.Nt + m;
            plan->fft_in[m][0] = plan->h_TrueWF[id] - plan->h_CurrWF[id];
            plan->fft_in[m][1] = 0.0f;
        }

        if(pt.getrank() == 0)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout << plan->fft_in[m][0] << endl;
            }

        }

        if(pt.getrank() == 1)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout << plan->fft_in[m][0] << endl;
            }
        }

        fftw_plan p = fftw_plan_dft_1d(Pa.Nt,
            plan->fft_in,
            plan->fft_out,
            FFTW_FORWARD,
            FFTW_ESTIMATE);

        fftw_execute(p);
        fftw_destroy_plan(p);

        if(pt.getrank() == 0)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout1 << plan->fft_in[m][0] << endl;
            }
            //fout1 << plan->fft_in[n][0] << endl;
        }

        if(pt.getrank() == 1)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout1 << plan->fft_in[m][0] << endl;
            }

        }


        if(pt.getrank() == 0)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout2 << plan->fft_out[m][0] << endl;
            }
            //fout2 << plan->fft_out[n][0] << endl;
        }

        if(pt.getrank() == 1)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout2 << plan->fft_out[m][0] << endl;
            }
            //fout2 << plan->fft_out[n][0] << endl;
        }

        // 对fft的结果求共轭
        for (uint m = 0; m < Pa.Nt; m++)
        {
            plan->fft_out[m][0] = plan->fft_out[m][0] * (1.0f / Pa.Nt);
            plan->fft_out[m][1] = plan->fft_out[m][1] * (-1.0f / Pa.Nt);
        }


        if(pt.getrank() == 0)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout3 << plan->fft_out[m][0] << endl;
            }

        }

        if(pt.getrank() == 1)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout3 << plan->fft_out[m][0] << endl;
            }

        }



        // 对上述结果求取反变换
        fftw_plan pp = fftw_plan_dft_1d(Pa.Nt,
            plan->fft_out,
            plan->fft_in,
            FFTW_BACKWARD,
            FFTW_ESTIMATE);

        fftw_execute(pp);


        if(pt.getrank() == 0)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout4 << plan->fft_in[m][0] << endl;
            }

        }

        if(pt.getrank() == 1)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout4 << plan->fft_in[m][0] << endl;
            }

        }


        if(pt.getrank() == 0)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout5 << plan->fft_out[m][0] << endl;
            }
        }

        if(pt.getrank() == 1)
        {
            for(int m = 0; m < Pa.Nt; ++m)
            {
                fout5 << plan->fft_out[m][0] << endl;
            }
        }


        for (uint m = 0; m < Pa.Nt; m++)
        {
            id = n * Pa.Nt + m;
            plan->h_ResWF[id] = plan->fft_in[m][0];
        }

        fftw_destroy_plan(pp);
    }


}


/*------------------------------------------------------------------------
function: PostProcessGrad

Information:
梯度后处理

Inputs:
Pa：正演参数
h_Grad：梯度
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
void PostProcessGrad(AFDP3D Pa,
                     float *h_Grad,
                     float *h_Vp,
                     const Partition& pt)
{
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id1 = 0, id2 = 0;

    H_Border temph_Vp = pt.geth_Vp();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_y = pt.getinteriorLength_y();
    uint interiorLength_z = pt.getinteriorLength_z();

    uint length_x = pt.getblockLength_x();
    uint length_y = pt.getblockLength_y();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_y = pt.getinteriormin_y() - pt.getindexmin_y();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();

    uint interiormin_y = pt.getinteriormin_y();
    uint interiormin_z = pt.getinteriormin_z();

    for (uint iz = 0; iz < interiorLength_z; iz++)
    {
        for (uint iy = 0; iy < interiorLength_y; iy++)
        {
            for (uint ix = 0; ix < interiorLength_x; ix++)
            {
                id1 = iz * interiorLength_y * interiorLength_x +
                    iy * interiorLength_x + ix;
                id2 = (iz + gap_z + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder + gap_y) * temph_Vp.length_x + ix + gap_x + temph_Vp.leftborder;

                h_Grad[id1] = (powf((float)(interiormin_z + iz - Pa.PMLz), PowerGrad) /
                    (h_Vp[id2] * sqrtf(h_Vp[id2]))) * h_Grad[id1];

                if (interiormin_z + iz < Pa.PMLz + 5)
                {
                    h_Grad[id1] = 0.0f;
                }
            }
        }
    }
}

/*------------------------------------------------------------------------
function: UpdateVpInnerGrid

Information:
更新内部网格的速度值

Inputs:
Pa：正演参数
h_Grad：梯度
h_Vp: 当前速度的平方
e: 步长
------------------------------------------------------------------------*/
void UpdateVpInnerGrid(AFDP3D Pa,
                       float *h_Vp,
                       float *h_Grad,
                       float e,
                       const Partition& pt)
{
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id1 = 0, id2 = 0;

    H_Border temph_Vp =pt.geth_Vp();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_y = pt.getinteriorLength_y();
    uint interiorLength_z = pt.getinteriorLength_z();

//    uint length_x = pt.getblockLength_x();
//    uint length_y = pt.getblockLength_y();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_y = pt.getinteriormin_y() - pt.getindexmin_y();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();

    for (uint iz = 0; iz < interiorLength_z; iz++)
    {
        for (uint iy = 0; iy < interiorLength_y; iy++)
        {
            for (uint ix = 0; ix < interiorLength_x; ix++)
            {
                id1 = iz * interiorLength_y * interiorLength_x +
                    iy * interiorLength_x + ix;
                id2 = (iz + gap_z + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + gap_y + temph_Vp.backborder) * temph_Vp.length_x + ix + gap_x + temph_Vp.leftborder;

                h_Vp[id2] = powf((e * h_Grad[id1] + sqrtf(fabs(h_Vp[id2]))), 2.0f);
            }
        }
    }
}

/*------------------------------------------------------------------------
function: UpdateVpPML

Information:
更新PML边界内网格的速度值

Inputs:
Pa：正演参数
h_Grad：梯度
h_Vp: 当前速度的平方
e: 步长
------------------------------------------------------------------------*/
void UpdateVpPML(AFDP3D Pa,
                 float *h_Vp,
                 float *h_Grad,
                 float e,
                 const Partition& pt)
{
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
//    uint nny = Pa.Ny + 2 * Pa.PMLy;
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    uint id = 0;
    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    H_Border temph_Vp = pt.geth_Vp();

    uint min_x = pt.getindexmin_x();
    uint min_y = pt.getindexmin_y();
    uint min_z = pt.getindexmin_z();
    uint max_x = pt.getindexmax_x();
    uint max_y = pt.getindexmax_y();
    uint max_z = pt.getindexmax_z();

    for (uint iz = 0; iz < block_z; iz++)
    {
        for (uint iy = 0; iy < block_y; iy++)
        {
            for (uint ix = 0; ix < block_x; ix++)
            {
                id = (iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder;
                int temp_z = iz + min_z;
                int temp_y = iy + min_y;
                int temp_x = ix + min_x;

                // z方向
                if(temp_z < Pa.PMLz)
                {
                    if(!temph_Vp.bottomborder)
                        h_Vp[id] = h_Vp[(Pa.PMLz - min_z + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                    else
                        h_Vp[id] = h_Vp[(block_z + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                }
                if(temp_z > Pa.PMLz + Pa.Nz - 1)
                {
                    if(!temph_Vp.topborder)
                        h_Vp[id] = h_Vp[(Pa.PMLz + Pa.Nz - 1 - min_z + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                    else
                        h_Vp[id] = h_Vp[(iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                }
                // y方向
                if(temp_y < Pa.PMLy)
                {
                    if(!temph_Vp.frontborder)
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (Pa.PMLy - min_y + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                    else
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (block_y + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                }
                if(temp_y > Pa.PMLy + Pa.Ny - 1)
                {
                    if(!temph_Vp.backborder)
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (Pa.PMLy + Pa.Ny - 1 - min_y + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                    else
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + ix + temph_Vp.leftborder];
                }
                // x方向
                if(temp_x < Pa.PMLx)
                {
                    if(!temph_Vp.rightborder)
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + (Pa.PMLx - min_x + temph_Vp.leftborder)];
                    else
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + block_x + temph_Vp.leftborder];
                }
                if(temp_x > Pa.PMLx + Pa.Nx - 1)
                {
                    if(!temph_Vp.leftborder)
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + (Pa.PMLx + Pa.Nx - 1 - min_x + temph_Vp.leftborder)];
                    else
                        h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x];
                }


//                // z方向
//                if (iz < Pa.PMLz)
//                {
//                    h_Vp[id] = h_Vp[Pa.PMLz * nny * nnx + iy * nnx + ix];
//                }
//                if (iz > Pa.Nz + Pa.PMLz - 1)
//                {
//                    h_Vp[id] = h_Vp[(Pa.Nz + Pa.PMLz - 1) * nny * nnx + iy * nnx + ix];
//                }
//                // y方向
//                if (iy < Pa.PMLy)
//                {
//                    h_Vp[id] = h_Vp[iz * nny * nnx + Pa.PMLy * nnx + ix];
//                }
//                if (iy > Pa.Ny + Pa.PMLy - 1)
//                {
//                    h_Vp[id] = h_Vp[iz * nny * nnx + (Pa.Ny + Pa.PMLy - 1) * nnx + ix];
//                }
//                // x方向
//                if (ix < Pa.PMLx)
//                {
//                    h_Vp[id] = h_Vp[iz * nny * nnx + iy * nnx + Pa.PMLx];
//                }
//                if (ix > Pa.Nx + Pa.PMLx - 1)
//                {
//                    h_Vp[id] = h_Vp[iz * nny * nnx + iy * nnx + Pa.Nx + Pa.PMLx - 1];
//                }
            }
        }
    }
}

/*------------------------------------------------------------------------
function: AddResidual

Information:
残差反传加震源

Inputs:
Pa：正演参数
h_ResWF：残差波场
h_U_r: 反传波场
h_re: 检波器位置
nstep: 时刻
Rn: 检波器个数
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
void AddResidual(AFDP3D Pa,
                 float *h_ResWF,
                 float *h_U_r,
                 RL *h_re,
                 uint nstep,
                 uint Rn,
                 float *h_Vp,
                 const Partition &pt,
                 int it)
{
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id = 0;

    H_Border temph_Vp =pt.geth_Vp();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_y = pt.getindexmin_y();
    uint indexmin_z = pt.getindexmin_z();

    vector<tuple<uint, uint, uint>> vec = pt.getRL();
    tuple<uint, uint, uint> temp;
    vector<tuple<uint, uint, uint>>::iterator begin = vec.begin();


    for(uint ir = 0; ir < vec.size(); ++ir)
    {
        temp = *(begin + ir);

        //if(pt.getrank() == 0)
        //fouth_u << get<0>(temp) << " " << get<1>(temp) << " " << get<2>(temp) << endl;
        h_U_r[get<2>(temp) * block_y * block_x + get<1>(temp) * block_x + get<0>(temp)] += -1.0f * h_ResWF[ir * Pa.Nt + nstep]
                * (Pa.dt * Pa.dt * h_Vp[(get<2>(temp) + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (get<1>(temp) + temph_Vp.backborder) * temph_Vp.length_x + get<0>(temp) + temph_Vp.leftborder]);


    }


begin = vec.begin();
    if(pt.getrank() == 0 && it == 1)
    {
        ofstream fout("h_U_r00.txt");
        for(int ir = 0; ir < Rn; ++ir)
        {
            temp = *(begin + ir);
            fout << get<2>(temp) + indexmin_z << " " << get<1>(temp) + indexmin_y << " " << get<0>(temp) + indexmin_x << endl;
            fout << h_U_r[get<2>(temp) * block_y * block_x + get<1>(temp) * block_x + get<0>(temp)] << " ";
        }
    }

    if(pt.getrank() == 1 && it == 1)
    {
        ofstream fout("h_U_r01.txt");
        for(int ir = 0; ir < Rn; ++ir)
        {
            temp = *(begin + ir);
            fout << get<2>(temp) + indexmin_z << " " << get<1>(temp) + indexmin_y << " " << get<0>(temp) + indexmin_x << endl;
            fout << h_U_r[get<2>(temp) * block_y * block_x + get<1>(temp) * block_x + get<0>(temp)] << " ";
        }
    }

//    for (uint ir = 0; ir < Rn; ir++)
//    {
//        id = h_re[ir].Rz * nny * nnx +
//            h_re[ir].Ry * nnx + h_re[ir].Rx;

//        h_U_r[id] += -1.0f * h_ResWF[ir * Pa.Nt + nstep]
//            * (Pa.dt * Pa.dt * h_Vp[id]);
//    }

}

/*------------------------------------------------------------------------
function: ComplexConjugate

Information:
求取复数向量中每一维的共轭

Inputs:
h_data1: 原始复向量
h_data2: 共轭后的复向量
------------------------------------------------------------------------*/
void ComplexConjugate(Complex *h_data1,
                      Complex *h_data2,
                      uint len)
{
    for (uint n = 0; n < len; n++)
    {
        h_data2[n].x = h_data1[n].x;
        h_data2[n].y = -1.0f * h_data1[n].y;
    }
}

/*------------------------------------------------------------------------
function: ComplexConjugate

Information:
求取复数向量中每一维的模的平方

Inputs:
h_data: 原始复向量
h_abs: 复向量每一维的模的平方
------------------------------------------------------------------------*/
void CalSquareAbs(float *h_abs,
                  Complex *h_data,
                  uint len)
{
    for (uint n = 0; n < len; n++)
    {
        h_abs[n] = h_data[n].x * h_data[n].x
            + h_data[n].y * h_data[n].y;
    }
}

/*------------------------------------------------------------------------
function: ComplexMulConst

Information:
复数向量乘以常数

Inputs:
h_data1: 原始复向量
h_data2: 运算后的复向量
a：常复数
------------------------------------------------------------------------*/
void ComplexMulConst(Complex *h_data1,
                     Complex *h_data2,
                     Complex a,
                     uint len)
{
    float temp1 = 0.0f;
    float temp2 = 0.0f;

    for (uint n = 0; n < len; n++)
    {
        temp1 = h_data1[n].x * a.x - h_data1[n].y * a.y;
        temp2 = h_data1[n].x * a.y + h_data1[n].y * a.x;

        h_data2[n].x = temp1;
        h_data2[n].y = temp2;
    }
}

/*------------------------------------------------------------------------
function: StepCalGradFreq

Information:
一步在频率域计算梯度

Inputs:
Pa：正演参数
h_Grad：梯度
h_WF_f：正传频率域波场
h_WF_b：反传频率域波场
------------------------------------------------------------------------*/
void StepCalGradFreq(AFDP3D Pa,
                     float *h_Grad,
                     Complex *h_WF_f,
                     Complex *h_WF_b,
                     const Partition& pt)
{
    uint id = 0;

    H_Border temph_U = pt.geth_U();

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_y = pt.getinteriorLength_y();
    uint interiorlength_z = pt.getinteriorLength_z();

    uint length_x = pt.getblockLength_x();
    uint length_y = pt.getblockLength_y();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_y = pt.getinteriormin_y() - pt.getindexmin_y();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();

    for(uint iz = 0; iz < interiorlength_z; ++iz)
    {
        for(uint iy = 0; iy < interiorlength_y; ++iy)
        {
            for(uint ix = 0; ix < interiorlength_x; ++ix)
            {
                id = iz * interiorlength_y * interiorlength_x + iy * interiorlength_x + ix;

                h_Grad[id] += -1.0f * (h_WF_f[id].x * h_WF_b[id].x
                                        - h_WF_f[id].y * h_WF_b[id].y);
            }
        }
    }

//    for (uint iz = 0; iz < Pa.Nz; iz++)
//    {
//        for (uint iy = 0; iy < Pa.Ny; iy++)
//        {
//            for (uint ix = 0; ix < Pa.Nx; ix++)
//            {
//                id = iz * Pa.Ny * Pa.Nx + iy * Pa.Nx + ix;

//                h_Grad[id] += -1.0f * (h_WF_f[id].x * h_WF_b[id].x
//                    - h_WF_f[id].y * h_WF_b[id].y);
//            }
//        }
//    }
}

/*------------------------------------------------------------------------
function: StepCalGradFreq

Information:
一步计算频率域的波场值

Inputs:
Pa：正演参数
h_U: 时间域波场
freq：频率
nt: 时刻
h_U_freq: 频率域波场
------------------------------------------------------------------------*/
void StepCalWFFreq(AFDP3D Pa,
                   float *h_U,
                   float freq,
                   uint nt,
                   Complex *h_U_freq,
                   const Partition& pt)
{
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    H_Border temph_U = pt.geth_U();

    uint id1 = 0, id2 = 0;

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_y = pt.getinteriorLength_y();
    uint interiorlength_z = pt.getinteriorLength_z();

    uint length_x = pt.getblockLength_x();
    uint length_y = pt.getblockLength_y();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_y = pt.getinteriormin_y() - pt.getindexmin_y();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();

    for (uint iz = 0; iz < interiorlength_z; iz++)
    {
        for (uint iy = 0; iy < interiorlength_y; iy++)
        {
            for (uint ix = 0; ix < interiorlength_x; ix++)
            {
                id1 = iz * interiorlength_y * interiorlength_x + iy * interiorlength_x + ix;
                id2 = (iz + gap_z + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + gap_y + temph_U.backborder) * temph_U.length_x + ix + gap_x + temph_U.leftborder;

                h_U_freq[id1].x += h_U[id2] *
                    cosf(freq * Pa.dt * (float)(nt));
                h_U_freq[id1].y += -1.0f * h_U[id2] *
                    sinf(freq * Pa.dt * (float)(nt));
            }
        }
    }


//    for (uint iz = 0; iz < Pa.Nz; iz++)
//    {
//        for (uint iy = 0; iy < Pa.Ny; iy++)
//        {
//            for (uint ix = 0; ix < Pa.Nx; ix++)
//            {
//                id1 = iz * Pa.Ny * Pa.Nx + iy * Pa.Nx + ix;
//                id2 = (iz + Pa.PMLz) * nny * nnx +
//                    (iy + Pa.PMLy) * nnx + ix + Pa.PMLx;

//                h_U_freq[id1].x += h_U[id2] *
//                    cosf(freq * Pa.dt * (float)(nt));
//                h_U_freq[id1].y += -1.0f * h_U[id2] *
//                    sinf(freq * Pa.dt * (float)(nt));
//            }
//        }
//    }
}

/*------------------------------------------------------------------------
function: CalSGFreq

Information:
计算频率域炮集

Inputs:
Pa：正演参数
Nt：正演总共的时间步数
Rn：检波器个数
Freq：频率
h_SG：时间域炮集
h_SG_freq：频率域炮集
------------------------------------------------------------------------*/
void CalSGFreq(AFDP3D Pa,
               uint Nt,
               uint Rn,
               float Freq,
               float *h_SG,
               Complex *h_SG_freq)
{
    for (uint ir = 0; ir < Rn; ir++)
    {
        for (uint it = 0; it < Nt; it++)
        {
            h_SG_freq[ir].x += h_SG[ir * Nt + it] *
                cosf(Freq * it * Pa.dt);

            h_SG_freq[ir].x += -1.0f * h_SG[ir * Nt + it] *
                sinf(Freq * it * Pa.dt);

        }
    }
}

/*------------------------------------------------------------------------
function: ComplexDotMulReal

Information:
计算两个复向量点乘的实部

Inputs:
h_data：输出的实部数据
h_data_f1：输入复向量1
h_data_f2：输入复向量2
len：向量维数
------------------------------------------------------------------------*/
void ComplexDotMulReal(float *h_data,
                       Complex *h_data_f1,
                       Complex *h_data_f2,
                       uint len)
{
    for (uint n = 0; n < len; n++)
    {
        h_data[n] = h_data_f1[n].x * h_data_f2[n].x
            - h_data_f1[n].y * h_data_f2[n].y;
    }
}

/*------------------------------------------------------------------------
function: CalTrueWF

Information:
求取观测波场

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
sgs_t: 观测波场
------------------------------------------------------------------------*/
void CalTrueWF(AFDP3D Pa,
               IP *ip,
               CPUVs *plan,
               float *sgs_t,
               const Partition& pt, const MPI_Comm &mycomm)
{
    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnz = Pa.Nz + 2 * Pa.PMLz;

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_y = pt.getindexmin_y();
    uint indexmin_z = pt.getindexmin_z();

    uint indexmax_x = pt.getindexmax_x();
    uint indexmax_y = pt.getindexmax_y();
    uint indexmax_z = pt.getindexmax_z();

    H_Border temph_Vp = pt.geth_Vp();
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();

    uint pos_x = pt.get_in_blockPosition_x();
    uint pos_y = pt.get_in_blockPosition_y();
    uint pos_z = pt.get_in_blockPosition_z();

    uint RL_num = pt.getRL_num();

    float Wavelet = 0.0f;

    int in_rank = pt.get_in_rank();
    int rank = pt.getrank();

    // 给速度赋值
    memset((void *)plan->h_Vp,			0,	sizeof(float) * temph_Vp.length_x * temph_Vp.length_y * temph_Vp.length_z);

    for(int i = 0; i < block_z; ++i)
    {
        for(int j = 0; j < block_y; ++j)
        {
            memcpy(plan->h_Vp + (i + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (j + temph_Vp.backborder) * temph_Vp.length_x + temph_Vp.leftborder,	ip->TrueVp + i * block_y * block_x + j * block_x,	block_x * sizeof(float));
        }
    }

//    for(int i = 0; i < block_z; ++i)
//    {
//        for(int j = 0; j < block_y; ++j)
//        {
//            for(int k = 0; k < block_x; ++k)
//            {
//                cout << *(plan->h_Vp + (i + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (j + temph_Vp.backborder) * temph_Vp.length_x + temph_Vp.leftborder) << " ";
//            }
//        }
//    }

    uint is = rank % ip->ShotN;

//    for (uint is = 0; is < ip->ShotN; is++)
//    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_U_y,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_K_y,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * temph_U.length_z * temph_U.length_y * temph_U.length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_K,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);

        memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);

        if(RL_num)
        {
            memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        }

        // 给检波器的位置信息赋值
        //memcpy(plan->h_re,	ip->St[is].re, ip->St[is].rn * sizeof(RL));

        // 对时间点进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
            for(int i = 0; i < block_z; ++i)
            {
                for(int j = 0; j < block_y; ++j)
                {
                    // 波场时刻转换
                    memcpy(plan->h_U_past + i * block_y * block_x + j * block_x, plan->h_U_now + (i + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (j + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, sizeof(float) * block_x);
                    memcpy(plan->h_U_now + (i + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (j + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next + i * block_y * block_x + j * block_x, sizeof(float) * block_x);
                }
            }

//            if(it == 2)
//            {
//                if(pt.getrank() == 1)
//                for(int i = 0; i < block_z; ++i)
//                {
//                    for(int j = 0; j < block_y; ++j)
//                    {
//                        for(int k = 0; k < block_x; ++k)
//                            if(*(plan->h_U_now + (i + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (j + temph_U.backborder) * temph_U.length_x + temph_U.leftborder + k))
//                                cout << *(plan->h_U_now + (i + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (j + temph_U.backborder) * temph_U.length_x + temph_U.leftborder + k)<< "end" << endl ;
//                    }
//                }
//            }
            dataTransport(plan->h_U_now, pt, STEP_U, it, mycomm);

//            MPI_Barrier(MPI_COMM_WORLD);
//            if(rank == 0)
//            {
//                //MPI_Barrier(MPI_COMM_WORLD);
//                for(int i = 0; i < temph_U.length_z * temph_U.length_x * temph_U.length_y; ++i)
//                    if(*(plan->h_U_now + i))
//                    cout << *(plan->h_U_now + i) << " ";
//                //MPI_Barrier(MPI_COMM_WORLD);
//            }

//cout << "RRR" << in_rank << endl;


            if(RL_num)
            {
                // 一步记录炮集
                StepShotGather(Pa, plan->h_U_now, plan->h_TrueWF,
                    plan->h_re, it, pt);
            }




            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIy_U_y, plan->h_PHIz_U_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

//            MPI_Barrier(MPI_COMM_WORLD);
//            if(rank == ROOT_ID)
//            {
//                cout << "it=" << it << endl;
//            }



            // 一步更新波场V，K和W
            StepVKW(Pa, plan->h_U_now, plan->h_V, plan->h_K, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIy_U_y, plan->h_PHIz_U_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);



            dataTransport(plan->h_V, pt, STEP_V, it, mycomm);//send data
            dataTransport(plan->h_K, pt, STEP_K, it, mycomm);//send data
            dataTransport(plan->h_W, pt, STEP_W, it, mycomm);//send data





            // 一步更新V和W的卷积项
            StepPHIVKW(Pa, plan->h_V, plan->h_K, plan->h_W,
                plan->h_PHIx_V_x, plan->h_PHIy_K_y, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_K, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIy_K_y, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_By, plan->h_Bz,
                plan->h_Vp, pt, it);

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);
            //cout << Wavelet << endl;

            if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sy >= indexmin_y && ip->St[is].s.Sy <= indexmax_y && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
            {
                AddSource(Pa, plan->h_U_next, ip->St[is].s, Wavelet, plan->h_Vp, pt);
            }


        }

        if(RL_num)
        {
            // 输出炮集
            memcpy(sgs_t + is * (Pa.Nt * RL_num),
                plan->h_TrueWF,
                Pa.Nt * RL_num * sizeof(float));
        }

        if(rank == 0)
        {
            ofstream fout("sgs_t00.txt");
            for(int i = 0; i < RL_num; ++i)
            {
                for(int j = 0; j < Pa.Nt; ++j)
                {
                    fout << sgs_t[i * Pa.Nt + j] << " ";
                }
                fout << endl;
            }
        }
        if(rank == 1)
        {
            ofstream fout("sgs_t01.txt");
            for(int i = 0; i < RL_num; ++i)
            {
                for(int j = 0; j < Pa.Nt; ++j)
                {
                    fout << sgs_t[i * Pa.Nt + j] << " ";
                }
                fout << endl;
            }
        }

    //}//shotN

}

/*------------------------------------------------------------------------
function: CalGrad

Information:
求取梯度

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
sgs_t: 观测波场
sgs_c: 当前参数正演波场
sgs_r: 残差波场
It: 迭代到第几步
------------------------------------------------------------------------*/
void CalGrad(AFDP3D Pa,
             IP *ip,
             CPUVs *plan,
             float *sgs_t,
             float *sgs_c,
             float *sgs_r,
             uint It,
             const Partition& pt,
             const MPI_Comm& mycomm)
{
    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    H_Border temph_Vp = pt.geth_Vp();
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_y = pt.getindexmin_y();
    uint indexmin_z = pt.getindexmin_z();

    uint indexmax_x = pt.getindexmax_x();
    uint indexmax_y = pt.getindexmax_y();
    uint indexmax_z = pt.getindexmax_z();

    uint interior_min_z = pt.getinteriormin_z();
    uint interior_min_y = pt.getinteriormin_y();
    uint interior_min_x = pt.getinteriormin_x();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_y = pt.getinteriorLength_y();
    uint interiorLength_z = pt.getinteriorLength_z();

    int rank = pt.getrank();
    int in_rank = pt.get_in_rank();

    int rank_size = pt.getsize();
    int in_rank_size = pt.get_in_size();

    uint RecNum = pt.geth_Coord_length();
    int RL_num = pt.getRL_num();

    float Wavelet = 0.0f;
    Complex a;
    a.x = 0.0f;
    a.y = 0.0f;

    if(interiorLength_x && interiorLength_y && interiorLength_z)
    {
        // 把梯度变量空间初始化
        memset((void *)plan->h_Grad,		0,	sizeof(float) * interiorLength_z * interiorLength_y * interiorLength_x);
    }


    // 给速度赋值
    memset((void *)plan->h_Vp,			0,	sizeof(float) * temph_Vp.length_z * temph_Vp.length_y * temph_Vp.length_x);
    for(int i = 0; i < block_z; ++i)
    {
        for(int j = 0; j < block_y; ++j)
        {
            memcpy(plan->h_Vp + (temph_Vp.topborder + i) * temph_Vp.length_y * temph_Vp.length_x + (j + temph_Vp.backborder) * temph_Vp.length_x + temph_Vp.leftborder, ip->CurrVp + i * block_y * block_x + j * block_x, block_x * sizeof(float));
        }
    }

    MPI_Status status_send, status_recv;
    uint is = rank % ip->ShotN;
    // 对每一炮进行迭代
//    for (uint is = 0; is < ip->ShotN; is++)
//    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_U_y,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_K_y,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * temph_U.length_z * temph_U.length_y * temph_U.length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_K,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_PHIx_U_x_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_U_y_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_U_z_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIx_V_x_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_K_y_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_W_z_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_past_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_now_r,		0,	sizeof(float) * temph_U.length_z * temph_U.length_y * temph_U.length_x);
        memset((void *)plan->h_U_next_r,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_V_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_K_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_W_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);

        memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);

        if(RL_num)
        {
            memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_CurrWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_ResWF,		0,	sizeof(float) * RL_num * Pa.Nt);

            memset((void *)plan->h_Obj,			0,	sizeof(float) * RL_num * Pa.Nt);
        }

        if(interiorLength_z && interiorLength_y && interiorLength_x)
        {
            memset((void *)plan->h_U_freq,		0,	sizeof(Complex) * interiorLength_z * interiorLength_y * interiorLength_x * ip->FreqN);
            memset((void *)plan->h_U_freq_r,	0,	sizeof(Complex) * interiorLength_z * interiorLength_y * interiorLength_x * ip->FreqN);
        }

        // 给检波器位置赋值
        //memcpy(plan->h_re,	ip->St[is].re,	ip->St[is].rn * sizeof(RL));

        // 对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
//            MPI_Barrier(MPI_COMM_WORLD);
//            if(rank == ROOT_ID)
//            {
//                cout << "it=" << it << endl;
//            }

            for(int iz = 0; iz < block_z; ++iz)
            {
                for(int iy = 0; iy < block_y; ++iy)
                {
                    // 波场时刻转换
                    memcpy(plan->h_U_past + iz * block_y * block_x + iy * block_x, plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, sizeof(float) * block_x);
                    memcpy(plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next + iz * block_y * block_x + iy * block_x, sizeof(float) * block_x);
                }
            }

            dataTransport(plan->h_U_now, pt, STEP_U, it, mycomm);//send data

            if(RL_num)
            {
                // 一步记录炮集
                StepShotGather(Pa, plan->h_U_now, plan->h_CurrWF,
                    plan->h_re, it, pt);
            }


            if(interiorLength_x && interiorLength_y && interiorLength_z)
            {
                // 一步计算频率域波场
                for (uint nf = 0; nf < ip->FreqN; nf++)
                {
                    StepCalWFFreq(Pa, plan->h_U_now, StartFreq + (float)nf * InterFreq,
                        it, plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x, pt);
                }
            }


            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIy_U_y, plan->h_PHIz_U_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场V，K和W
            StepVKW(Pa, plan->h_U_now, plan->h_V, plan->h_K, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIy_U_y, plan->h_PHIz_U_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            dataTransport(plan->h_V, pt, STEP_V, it, mycomm);//send data
            dataTransport(plan->h_K, pt, STEP_K, it, mycomm);//send data
            dataTransport(plan->h_W, pt, STEP_W, it, mycomm);//send data

            // 一步更新V和W的卷积项
            StepPHIVKW(Pa, plan->h_V, plan->h_K, plan->h_W,
                plan->h_PHIx_V_x, plan->h_PHIy_K_y, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_K, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIy_K_y, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_By, plan->h_Bz,
                plan->h_Vp, pt, it);

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sy >= indexmin_y && ip->St[is].s.Sy <= indexmax_y && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
            {
                AddSource(Pa, plan->h_U_next, ip->St[is].s, Wavelet, plan->h_Vp, pt);
            }
        }

        if(RL_num)
        {
            // 输出当前模型正演的炮集
            memcpy(sgs_c + is * (Pa.Nt * RL_num),
                plan->h_CurrWF,
                Pa.Nt * RL_num * sizeof(float));

            // 读入观测波场的炮集
            memcpy(plan->h_TrueWF,
                sgs_t + is * (Pa.Nt * RL_num),
                Pa.Nt * RL_num * sizeof(float));

            // 计算残差波场
            StepResidualConj(Pa, plan, RL_num, pt);

            // 输出残差炮集
            memcpy(sgs_r + is * (Pa.Nt * RL_num),
                plan->h_ResWF,
                Pa.Nt * RL_num * sizeof(float));
        }


        if(RL_num)
        {
            if(pt.getrank() == 0)
            {
                ofstream fout_h("h_CurrWF00.txt");
                for(int i = 0; i < RL_num; ++i)
                {
                    for(int j = 0; j < Pa.Nt; ++j)
                    {
                        fout_h << *(plan->h_CurrWF + i * Pa.Nt + j) << " ";
                    }
                }
                fout_h << endl;
            }

            if(pt.getrank() == 1)
            {
                ofstream fout_h("h_CurrWF01.txt");
                for(int i = 0; i < RL_num; ++i)
                {
                    for(int j = 0; j < Pa.Nt; ++j)
                    {
                        fout_h << *(plan->h_CurrWF + i * Pa.Nt + j) << " ";
                    }
                }
                fout_h << endl;
            }
        }

        if(RL_num)
        {
            if(pt.getrank() == 0)
            {
                ofstream fout_h("h_TrueWF00.txt");
                for(int i = 0; i < RL_num; ++i)
                {
                    for(int j = 0; j < Pa.Nt; ++j)
                    {
                        fout_h << *(plan->h_TrueWF + i * Pa.Nt + j) << " ";
                    }
                }
                fout_h << endl;
            }

            if(pt.getrank() == 1)
            {
                ofstream fout_h("h_TrueWF01.txt");
                for(int i = 0; i < RL_num; ++i)
                {
                    for(int j = 0; j < Pa.Nt; ++j)
                    {
                        fout_h << *(plan->h_TrueWF + i * Pa.Nt + j) << " ";
                    }
                }
                fout_h << endl;
            }
        }


        if(RL_num)
        {
            if(pt.getrank() == 0)
            {
                ofstream fout_h("h_ResWF0.txt");
                for(int i = 0; i < RL_num; ++i)
                {
                    for(int j = 0; j < Pa.Nt; ++j)
                    {
                        fout_h << *(plan->h_ResWF + i * Pa.Nt + j) << " ";
                    }
                }
                fout_h << endl;
            }

            if(pt.getrank() == 1)
            {
                ofstream fout_h("h_ResWF1.txt");
                for(int i = 0; i < RL_num; ++i)
                {
                    for(int j = 0; j < Pa.Nt; ++j)
                    {
                        fout_h << *(plan->h_ResWF + i * Pa.Nt + j) << " ";
                    }
                }
                fout_h << endl;
            }
        }



        if(rank == 0)
        {
            ofstream fout_hvp("h_vp0.txt");
            for(int iz = 0; iz < block_z; ++iz)
            {
                for(int iy = 0; iy < block_y; ++iy)
                {
                    for(int ix = 0; ix < block_x; ++ix)
                    {
                        fout_hvp << *(plan->h_Vp + (iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder) << " ";
                    }
                }
            }
            fout_hvp << endl;
        }

        if(rank == 1)
        {
            ofstream fout_hvp("h_vp1.txt");
            for(int iz = 0; iz < block_z; ++iz)
            {
                for(int iy = 0; iy < block_y; ++iy)
                {
                    for(int ix = 0; ix < block_x; ++ix)
                    {
                        fout_hvp << *(plan->h_Vp + (iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + ix + temph_Vp.leftborder) << " ";
                    }
                }
            }
            fout_hvp << endl;
        }

        // 残差反传，对时间的循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
//            MPI_Barrier(MPI_COMM_WORLD);
//            if(rank == ROOT_ID)
//            {
//                cout << "iit=" << it << endl;
//            }
            for(int iz = 0; iz < block_z; ++iz)
            {
                for(int iy = 0; iy < block_y; ++iy)
                {
                    // 波场时刻转换
                    memcpy(plan->h_U_past_r + iz * block_y * block_x + iy * block_x, plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, sizeof(float) * block_x);
                    memcpy(plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next_r + iz * block_y * block_x + iy * block_x, sizeof(float) * block_x);
                }
            }

//            if(pt.getrank() == 1)
//            for(int i = 0; i < block_z; ++i)
//            {
//                for(int j = 0; j < block_y; ++j)
//                {
//                    for(int k = 0; k < block_x; ++k)
//                    {
//                        if(*(plan->h_U_now_r))
//                            cout << *(plan->h_U_now_r + i * block_y * block_x + j * block_x) << " ";
//                    }
//                }
//            }

            dataTransport(plan->h_U_now_r, pt, STEP_U, it, mycomm);//send data

            if(interiorLength_x && interiorLength_y && interiorLength_z)
            {
                // 一步计算频率域波场
                for (uint nf = 0; nf < ip->FreqN; nf++)
                {
                    StepCalWFFreq(Pa, plan->h_U_now_r, StartFreq + (float)nf * InterFreq,
                        it, plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x, pt);
                }
            }


            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now_r, plan->h_PHIx_U_x_r,
                plan->h_PHIy_U_y_r, plan->h_PHIz_U_z_r,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场V，K和W
            StepVKW(Pa, plan->h_U_now_r, plan->h_V_r, plan->h_K_r, plan->h_W_r,
                plan->h_PHIx_U_x_r, plan->h_PHIy_U_y_r, plan->h_PHIz_U_z_r,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            dataTransport(plan->h_V, pt, STEP_V, it, mycomm);//send data
            dataTransport(plan->h_K, pt, STEP_K, it, mycomm);//send data
            dataTransport(plan->h_W, pt, STEP_W, it, mycomm);//send data

            // 一步更新V和W的卷积项
            StepPHIVKW(Pa, plan->h_V_r, plan->h_K_r, plan->h_W_r,
                plan->h_PHIx_V_x_r, plan->h_PHIy_K_y_r, plan->h_PHIz_W_z_r,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next_r, plan->h_U_now_r, plan->h_U_past_r,
                plan->h_V_r, plan->h_K_r, plan->h_W_r, plan->h_PHIx_V_x_r,
                plan->h_PHIy_K_y_r, plan->h_PHIz_W_z_r,
                plan->h_Bx, plan->h_By, plan->h_Bz,
                plan->h_Vp, pt, it);

            if(it == 1)
            {
                if(rank == 0)
                {
                    ofstream fout("h_U_next_r0.txt");
                    ofstream fout1("h_U_now_r0.txt");
                    for(int iz = 0; iz < block_z; ++iz)
                    {
                        for(int iy = 0; iy < block_y; ++iy)
                        {
                            for(int ix = 0; ix < block_x; ++ix)
                            {
                                fout << *(plan->h_U_next_r + iz * block_y * block_x + iy * block_x + ix) << " ";
                                fout1 << *(plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + ix + temph_U.leftborder) << " ";
                            }
                        }
                    }
                }

                if(rank == 1)
                {
                    ofstream fout("h_U_next_r1.txt");
                    ofstream fout1("h_U_now_r1.txt");
                    for(int iz = 0; iz < block_z; ++iz)
                    {
                        for(int iy = 0; iy < block_y; ++iy)
                        {
                            for(int ix = 0; ix < block_x; ++ix)
                            {
                                fout << *(plan->h_U_next_r + iz * block_y * block_x + iy * block_x + ix) << " ";
                                fout1 << *(plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + ix + temph_U.leftborder) << " ";
                            }
                        }
                    }
                }
            }

            if(RL_num)
            {
                // 加残差反传震源
                AddResidual(Pa, plan->h_ResWF, plan->h_U_next_r,
                    plan->h_re, it, RL_num, plan->h_Vp, pt, it);
            }

            if(it == 1)
            {
                if(rank == 0)
                {
                    ofstream fout("h_U_next_r00.txt");
                    ofstream fout1("h_U_now_r00.txt");
                    for(int iz = 0; iz < block_z; ++iz)
                    {
                        for(int iy = 0; iy < block_y; ++iy)
                        {
                            for(int ix = 0; ix < block_x; ++ix)
                            {
                                fout << *(plan->h_U_next_r + iz * block_y * block_x + iy * block_x + ix) << " ";
                                fout1 << *(plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + ix + temph_U.leftborder) << " ";
                            }
                        }
                    }
                }

                if(rank == 1)
                {
                    ofstream fout("h_U_next_r01.txt");
                    ofstream fout1("h_U_now_r01.txt");
                    for(int iz = 0; iz < block_z; ++iz)
                    {
                        for(int iy = 0; iy < block_y; ++iy)
                        {
                            for(int ix = 0; ix < block_x; ++ix)
                            {
                                fout << *(plan->h_U_next_r + iz * block_y * block_x + iy * block_x + ix) << " ";
                                fout1 << *(plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + ix + temph_U.leftborder) << " ";
                            }
                        }
                    }
                }
            }
        }


        if(pt.getrank() == 0)
        {
            ofstream fout("h_U_freq0.txt");
            ofstream fout1("h_U_freq_r0.txt");
            for(int nf = 0; nf < ip->FreqN; ++nf)
            {
                for(int iz = 0; iz < interiorLength_z; ++iz)
                {
                    for(int iy = 0; iy < interiorLength_y; ++iy)
                    {
                        for(int ix = 0; ix < interiorLength_x; ++ix)
                        {
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                        }
                    }
                }
            }
        }

        if(pt.getrank() == 1)
        {
            ofstream fout("h_U_freq1.txt");
            ofstream fout1("h_U_freq_r1.txt");
            for(int nf = 0; nf < ip->FreqN; ++nf)
            {
                for(int iz = 0; iz < interiorLength_z; ++iz)
                {
                    for(int iy = 0; iy < interiorLength_y; ++iy)
                    {
                        for(int ix = 0; ix < interiorLength_x; ++ix)
                        {
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                        }
                    }
                }
            }
        }

        if(pt.getrank() == 2)
        {
            ofstream fout("h_U_freq2.txt");
            ofstream fout1("h_U_freq_r2.txt");
            for(int nf = 0; nf < ip->FreqN; ++nf)
            {
                for(int iz = 0; iz < interiorLength_z; ++iz)
                {
                    for(int iy = 0; iy < interiorLength_y; ++iy)
                    {
                        for(int ix = 0; ix < interiorLength_x; ++ix)
                        {
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                        }
                    }
                }
            }
        }

        if(pt.getrank() == 3)
        {
            ofstream fout("h_U_freq3 .txt");
            ofstream fout1("h_U_freq_r3.txt");
            for(int nf = 0; nf < ip->FreqN; ++nf)
            {
                for(int iz = 0; iz < interiorLength_z; ++iz)
                {
                    for(int iy = 0; iy < interiorLength_y; ++iy)
                    {
                        for(int ix = 0; ix < interiorLength_x; ++ix)
                        {
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout << (plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->x << " ";
                            fout1 << (plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix)->y << " ";
                        }
                    }
                }
            }
        }

        // 计算当前炮下的梯度，对每一个频率进行循环
        for (uint nf = 0; nf < ip->FreqN; nf++)
        {
            a.x = -1.0f * powf((float)(StartFreq + nf * InterFreq), 2.0f);
            a.y = 0.0f;

            if(interiorLength_x && interiorLength_y && interiorLength_z)
            {
                // 第1步，乘以-w^2
                ComplexMulConst(plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x,
                    plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x,
                    a,
                    interiorLength_z * interiorLength_y * interiorLength_x);

                // 第2步，计算负梯度方向
                StepCalGradFreq(Pa, plan->h_Grad,
                    plan->h_U_freq + nf * interiorLength_z * interiorLength_y * interiorLength_x,
                    plan->h_U_freq_r + nf * interiorLength_z * interiorLength_y * interiorLength_x, pt);
            }
        }

        float *temp_h_Grad = new float[interiorLength_x * interiorLength_y * interiorLength_z];
        if(interiorLength_x && interiorLength_y && interiorLength_z)
        {
            if(is == 0)
            {
                for(int i = 0; i < ip->ShotN - 1; ++i)
                {
                    MPI_Recv(temp_h_Grad, interiorLength_z * interiorLength_y * interiorLength_x, MPI_FLOAT, MPI_ANY_SOURCE, STEP_GRAD, MPI_COMM_WORLD, &status_recv);
                    for(int iz = 0; iz < interiorLength_z; ++iz)
                    {
                        for(int iy = 0; iy < interiorLength_y; ++iy)
                        {
                            for(int ix = 0; ix < interiorLength_x; ++ix)
                            {
                                //id1 = iz * interiorLength_x + ix;
                                plan->h_Grad[iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix] += temp_h_Grad[iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix];
                            }
                        }
                    }
                }

                for(int i = 1; i < ip->ShotN; ++i)
                {
                    MPI_Send(plan->h_Grad, interiorLength_z * interiorLength_y * interiorLength_x, MPI_FLOAT, rank + i, STEP_GRAD, MPI_COMM_WORLD);
                }
            }
            else
            {
                MPI_Send(plan->h_Grad, interiorLength_z * interiorLength_y * interiorLength_x, MPI_FLOAT, rank - (rank % ip->ShotN), STEP_GRAD, MPI_COMM_WORLD);
                MPI_Recv(plan->h_Grad, interiorLength_z * interiorLength_y * interiorLength_x, MPI_FLOAT, rank - (rank % ip->ShotN), STEP_GRAD, MPI_COMM_WORLD, &status_recv);
            }
        }
        delete temp_h_Grad;

        int max_rank_RL = in_rank, min_rank_RL = in_rank;

        if(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx)
        {//cout << rank << endl;
            min_rank_RL = in_rank;
            if(min_rank_RL != ROOT_ID)
                MPI_Send(&min_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MIN_RL, mycomm);
            //cout << "AAAA" << min_rank_RL << " " << max_rank_RL << "AAAA" << endl;
        }
        if(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1)
        {//cout << rank << endl;
            max_rank_RL = in_rank;
            if(max_rank_RL != ROOT_ID)
                MPI_Send(&max_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MAX_RL, mycomm);
        }
//cout << "test" << rank << endl;
        if(in_rank == ROOT_ID)
        {
            if(!(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx))
                MPI_Recv(&min_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MIN_RL, mycomm, &status_recv);
            if(!(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1))
                MPI_Recv(&max_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MAX_RL, mycomm, &status_recv);
        }

        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&max_rank_RL, 1, MPI_INT, ROOT_ID, mycomm);
        MPI_Bcast(&min_rank_RL, 1, MPI_INT, ROOT_ID, mycomm);

//cout << "max_RL" << max_rank_RL << endl;
        float *buf_obj = new float[1];
        // 求取目标函数
        if(RL_num)
        {
            if(in_rank != min_rank_RL)
            {
                MPI_Recv(buf_obj, 1, MPI_FLOAT, in_rank - 1, STEP_OBJ, mycomm, &status_recv);
                ip->ObjIter[It] = *buf_obj;
                for (uint m = 0; m < RL_num * Pa.Nt; m++)
                {
                    //cout << *ip->ObjIter << endl;
                    ip->ObjIter[It] += 0.5f * powf(plan->h_ResWF[m], 2.0f);//ci shi It = 0// zhe li suan de zhi shi jubu de han shu
                }
                if(in_rank != max_rank_RL)
                {
                    MPI_Send(&ip->ObjIter[It], 1, MPI_FLOAT, in_rank + 1, STEP_OBJ, mycomm);
                }
            }
            else
            {
                for (uint m = 0; m < RL_num * Pa.Nt; m++)
                {
                    //cout << *ip->ObjIter << endl;
                    ip->ObjIter[It] += 0.5f * powf(plan->h_ResWF[m], 2.0f);//ci shi It = 0// zhe li suan de zhi shi jubu de han shu
                }
                //cout << "rank=" << rank << " " << ip->ObjIter[It] << endl;
                if(in_rank != max_rank_RL)
                    MPI_Send(&ip->ObjIter[It], 1, MPI_FLOAT, in_rank + 1, STEP_OBJ, mycomm);
            }
        }

        MPI_Bcast(&ip->ObjIter[It], 1, MPI_FLOAT, max_rank_RL, mycomm);

//cout << "zzz" << rank << endl;
        if(in_rank == ROOT_ID)
        {
            if(rank != ROOT_ID)
            {
                //cout << "zzz" << rank << endl;
                MPI_Send(&ip->ObjIter[It], 1, MPI_FLOAT, ROOT_ID, STEP_OBJ, MPI_COMM_WORLD);
            }
            else
            {
                for(int i = 1; i < ip->ShotN; ++i)
                {

                    //cout << "ooooooooo" << ip->ObjIter[It] << endl;
                    MPI_Recv(buf_obj, 1, MPI_FLOAT, MPI_ANY_SOURCE, STEP_OBJ, MPI_COMM_WORLD, &status_recv);
                    //cout << "ooooooooo" << *buf_obj << endl;
                    ip->ObjIter[It] += *buf_obj;
                    //cout << "ooooooooo" << ip->ObjIter[It] << endl;
                }
            }
        }
//cout << "test" << rank << endl;
        delete buf_obj;

//        if(rank >= min_rank_RL && rank <= max_rank_RL)
//cout << *ip->ObjIter << " " << rank << endl;

//        if(rank != ROOT_ID)
//        {
//            MPI_Isend(&ip->ObjIter[It], 1, MPI_FLOAT, ROOT_ID, STEP_OBJ, MPI_COMM_WORLD, &request_send);
//        }
//        else
//        {
//            for(uint i = 0; i < rank_size - 1; ++i)
//            {
//                MPI_Irecv(buf_obj, 1, MPI_FLOAT, MPI_ANY_SOURCE, STEP_OBJ, MPI_COMM_WORLD, &request_recv_OBJ);
//                MPI_Wait(&request_recv_OBJ, &status_recv);
//                //cout << *buf_obj << endl;
//                ip->ObjIter[It] += *buf_obj;
//            }
//        }

        //MPI_Barrier(MPI_COMM_WORLD);
//        if(rank == max_rank_RL)
//        {
//            cout << ip->ObjIter[It] << " "  << rank << endl;
//        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&ip->ObjIter[It], 1, MPI_FLOAT, ROOT_ID, MPI_COMM_WORLD);

        //MPI_Barrier(MPI_COMM_WORLD);
//        if(rank == ROOT_ID)
//        {
//            cout << "ObjIter=" << " " << ip->ObjIter[It] << endl;
//        }



//        // 计算目标函数值
//        for (uint m = 0; m < RL_num * Pa.Nt; m++)
//        {
//            ip->ObjIter[It] += 0.5f * powf(plan->h_ResWF[m], 2.0f);
//        }
//    }//shotN

    for(int iz = 0; iz < interiorLength_z; ++iz)
    {
        for(int iy = 0; iy < interiorLength_y; ++iy)
        {
            // 输出梯度
            memcpy(ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x, plan->h_Grad + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x, interiorLength_x * sizeof(float));
        }
    }

    if(pt.getrank() == 0)
    {
        ofstream fout("ip->GradVp0.txt");
        for(int iz = 0; iz < interiorLength_z; ++iz)
        {
            for(int iy = 0; iy < interiorLength_y; ++iy)
            {
                for(int ix = 0; ix < interiorLength_x; ++ix)
                {
                    fout << *(ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix) << " ";
                }
            }
        }
    }

    if(pt.getrank() == 1)
    {
        ofstream fout("ip->GradVp1.txt");
        for(int iz = 0; iz < interiorLength_z; ++iz)
        {
            for(int iy = 0; iy < interiorLength_y; ++iy)
            {
                for(int ix = 0; ix < interiorLength_x; ++ix)
                {
                    fout << *(ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix) << " ";
                }
            }
        }
    }

    if(pt.getrank() == 2)
    {
        ofstream fout("ip->GradVp2.txt");
        for(int iz = 0; iz < interiorLength_z; ++iz)
        {
            for(int iy = 0; iy < interiorLength_y; ++iy)
            {
                for(int ix = 0; ix < interiorLength_x; ++ix)
                {
                    fout << *(ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix) << " ";
                }
            }
        }
    }

    if(pt.getrank() == 3)
    {
        ofstream fout("ip->GradVp3.txt");
        for(int iz = 0; iz < interiorLength_z; ++iz)
        {
            for(int iy = 0; iy < interiorLength_y; ++iy)
            {
                for(int ix = 0; ix < interiorLength_x; ++ix)
                {
                    fout << *(ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x + ix) << " ";
                }
            }
        }
    }

}

/*------------------------------------------------------------------------
function: CalStepLength

Information:
计算迭代步长

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
sgs_t: 观测波场
sgs_c: 当前正演波场
e: 试探步长
------------------------------------------------------------------------*/
void CalStepLength(AFDP3D Pa,
                   IP *ip,
                   CPUVs *plan,
                   float *sgs_t,
                   float *sgs_c,
                   float e,
                   const Partition& pt,
                   const MPI_Comm& mycomm)
{
    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_y = pt.getinteriorLength_y();
    uint interiorLength_z = pt.getinteriorLength_z();

    uint interior_min_x = pt.getinteriormin_x();
    uint interior_min_y = pt.getinteriormin_y();
    uint interior_min_z = pt.getinteriormin_z();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_y = pt.getindexmin_y();
    uint indexmin_z = pt.getindexmin_z();

    uint indexmax_x = pt.getindexmax_x();
    uint indexmax_y = pt.getindexmax_y();
    uint indexmax_z = pt.getindexmax_z();

    H_Border temph_Vp = pt.geth_Vp();
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();

    MPI_Request request_send_U, request_send_V, request_send_W, request_recv_U, request_recv_V, request_recv_W, request_send, request_recv_MAX;
    MPI_Status status_send, status_recv;

    int rank = pt.getrank();
    int in_rank = pt.get_in_rank();

    int rank_size = pt.getsize();
    int in_rank_size = pt.get_in_size();

    uint sumblock_x = pt.getsumBlock_x();
    //uint sumblock_z = pt.getsumBlock_z();

    uint RL_num = pt.getRL_num();
    float *max_buf = new float[rank_size];

    // 给相关参数赋值为0
    memset((void *)plan->h_Vp,			0,	sizeof(float) * temph_Vp.length_z * temph_Vp.length_y * temph_Vp.length_x);
    memset((void *)plan->h_SumFenmu,	0,	sizeof(float));
    memset((void *)plan->h_SumFenzi,	0,	sizeof(float));

    float Wavelet = 0.0f;
    float MV = 0;

    // 将梯度进行归一化
    if(interiorLength_z && interiorLength_y && interiorLength_x)
    {

        MV = MaxValue(ip->GradVp, interiorLength_z * interiorLength_y * interiorLength_x);
    }
    //cout << MV << " " << pt.getrank() << endl;


    if(rank != ROOT_ID)
    {
        MPI_Send(&MV, 1, MPI_FLOAT, ROOT_ID, STEP_MAX, MPI_COMM_WORLD);
    }
    else
    {
        max_buf[rank_size - 1] = MV;
        if(rank_size > 1)
        {
            for(uint i = 0; i < rank_size - 1; ++i)
            {
                MPI_Recv(max_buf + i, 1, MPI_FLOAT, MPI_ANY_SOURCE, STEP_MAX, MPI_COMM_WORLD, &status_recv);
            }
        }

        MV = MaxValue(max_buf, rank_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MV, 1, MPI_FLOAT, ROOT_ID, MPI_COMM_WORLD);

    if(rank == ROOT_ID)
    cout << "max" << MV << endl;


    if (MV > 1.0e-10f)
    {
        for (uint m = 0; m < interiorLength_z * interiorLength_y * interiorLength_x; m++)
        {
            ip->GradVp[m] /= MV;
        }
    }

    // 生成试探速度模型
    if(interiorLength_z && interiorLength_y && interiorLength_x)
    {

        memset((void *)plan->h_Grad, 0, sizeof(float) * interiorLength_z * interiorLength_y * interiorLength_x);
        for(int iz = 0; iz < interiorLength_z; ++iz)
        {
            for(int iy = 0; iy < interiorLength_y; ++iy)
            {
                memcpy(plan->h_Grad + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x, ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x, sizeof(float) * interiorLength_x);
            }
        }
    }

    for(int iz = 0; iz < block_z; ++iz)
    {
        for(int iy = 0; iy < block_y; ++iy)
        {
            memcpy(plan->h_Vp + (iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x, ip->CurrVp + iz * block_y * block_x + iy * block_x, sizeof(float) * block_x);
        }
    }

    UpdateVpInnerGrid(Pa, plan->h_Vp, plan->h_Grad, e, pt);

    vector<uint> trans_h_Vp = pt.gettrans_h_Vp();
    auto begin = trans_h_Vp.begin();

    if(*(begin + TOP))
    {
        dataTransport_Vp(plan->h_Vp, pt, BOTTOM_TO_TOP, Pa, mycomm);
    }
    if(*(begin + LEFT))
    {
        dataTransport_Vp(plan->h_Vp, pt, RIGHT_TO_LEFT, Pa, mycomm);
    }
    if(*(begin + BOTTOM))
    {
        dataTransport_Vp(plan->h_Vp, pt, TOP_TO_BOTTOM, Pa, mycomm);
    }
    if(*(begin + RIGHT))
    {
        dataTransport_Vp(plan->h_Vp, pt, LEFT_TO_RIGHT, Pa, mycomm);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    dataGather(plan->h_Vp, pt, STEP_VP, mycomm);

    UpdateVpPML(Pa, plan->h_Vp, plan->h_Grad, e, pt);




    //uint is = rank % ip->ShotN;



    int *ranks = new int[pt.get_in_sumBlock_x()];
    //int cpu_index_z = (Pa.PMLz + 2) / pt.getblockLength_z();
    int max_rank_RL = in_rank, min_rank_RL = in_rank;
    float buf[2];// = {fenzi, fenmu};
    // 求取最终的步长

    if(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx)
    {
        min_rank_RL = in_rank;
        if(min_rank_RL != ROOT_ID)
            MPI_Send(&min_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MIN_RL, mycomm);
    }
    if(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1)
    {
        max_rank_RL = in_rank;
        if(max_rank_RL != ROOT_ID)
            MPI_Send(&max_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MAX_RL, mycomm);
    }

    if(in_rank == ROOT_ID)
    {
        if(!(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx))
            MPI_Recv(&min_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MIN_RL, mycomm, &status_recv);
        if(!(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1))
            MPI_Recv(&max_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MAX_RL, mycomm, &status_recv);
    }
    MPI_Barrier(mycomm);
    //cout << rank << "in" << endl;
    MPI_Bcast(&max_rank_RL, 1, MPI_INT, ROOT_ID, mycomm);
    MPI_Bcast(&min_rank_RL, 1, MPI_INT, ROOT_ID, mycomm);

    if(rank == ROOT_ID)
    {
        cout << "max_rank" << max_rank_RL << endl;
        cout << "min_rank" << min_rank_RL << endl;
    }

    uint is = rank % ip->ShotN;

    // 对炮进行循环
//    for (uint is = 0; is < ip->ShotN; is++)
//    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_U_y,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIy_K_y,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * temph_U.length_z * temph_U.length_y * temph_U.length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_y * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_K,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);
        memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_y * h_VW.length_x);

        memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);

        if(RL_num)
        {
            memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_CurrWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_ResWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_TrailWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        }

        // 给检波器位置赋值
        //memcpy(plan->h_re,	ip->St[is].re, ip->St[is].rn * sizeof(RL));


        // 对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if(rank == ROOT_ID)
            {
                cout << "stepit=" << it << endl;
            }

            for(uint iz = 0; iz < block_z; ++iz)
            {
                for(uint iy = 0; iy < block_y; ++iy)
                {
                    // 波场时刻转换
                    memcpy(plan->h_U_past + iz * block_y * block_x + iy * block_x, plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, sizeof(float) * block_x);
                    memcpy(plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_y * temph_U.length_x + (iy + temph_U.backborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next + iz * block_y * block_x + iy * block_x, sizeof(float) * block_x);
                }
            }

            dataTransport(plan->h_U_now, pt, STEP_U, it, mycomm);//send data
            MPI_Barrier(mycomm);

            // 一步记录炮集
            StepShotGather(Pa, plan->h_U_now, plan->h_TrailWF,
                plan->h_re, it, pt);

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIy_U_y, plan->h_PHIz_U_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场V，K和W
            StepVKW(Pa, plan->h_U_now, plan->h_V, plan->h_K, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIy_U_y, plan->h_PHIz_U_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            dataTransport(plan->h_V, pt, STEP_V, it, mycomm);//send data
            dataTransport(plan->h_K, pt, STEP_K, it, mycomm);//send data
            dataTransport(plan->h_W, pt, STEP_W, it, mycomm);//send data

            // 一步更新V和W的卷积项
            StepPHIVKW(Pa, plan->h_V, plan->h_K, plan->h_W,
                plan->h_PHIx_V_x, plan->h_PHIy_K_y, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_By, plan->h_Bz, pt, it);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_K, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIy_K_y, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_By, plan->h_Bz,
                plan->h_Vp, pt, it);

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sy >= indexmin_y && ip->St[is].s.Sy <= indexmax_y && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
            {
                AddSource(Pa, plan->h_U_next, ip->St[is].s, Wavelet, plan->h_Vp, pt);
            }
        }

        if(RL_num)
        {
            // 重新读取当前模型对应的波场炮集
            memcpy(plan->h_CurrWF,
                sgs_c + is * (Pa.Nt * RL_num),
                Pa.Nt * RL_num * sizeof(float));

            // 计算上述两个炮集的残差
            StepResidual(Pa, plan->h_TrailWF, plan->h_CurrWF,
                plan->h_TrailWF, RL_num);

            // 重新读取观测炮集
            memcpy(plan->h_TrueWF,
                sgs_t + is * (Pa.Nt * RL_num),
                Pa.Nt * RL_num * sizeof(float));

            // 计算当前模型对应的波场炮集域观测炮集的差
            StepResidual(Pa, plan->h_TrueWF, plan->h_CurrWF,
                plan->h_ResWF, RL_num);
        }

        // 对频率进行循环
        for (uint nf = 0; nf < ip->FreqN; nf++)
        {
            // 初始化变量
            memset((void *)plan->h_ResTrail_freq,	0,	sizeof(Complex) * RL_num);
            memset((void *)plan->h_ResCurr_freq,	0,	sizeof(Complex) * RL_num);

            if(RL_num)
            {
                // 计算频率域炮集
                CalSGFreq(Pa, Pa.Nt, RL_num,
                    (float)(StartFreq + nf * InterFreq),
                    plan->h_TrailWF, plan->h_ResTrail_freq);

                CalSGFreq(Pa, Pa.Nt, RL_num,
                    (float)(StartFreq + nf * InterFreq),
                    plan->h_ResWF, plan->h_ResCurr_freq);

                // 计算步长公式中的分子
                memset((void *)plan->h_Temp,	0,	sizeof(float) * RL_num);

                // 求取h_ResCurr_freq的共轭
                ComplexConjugate(plan->h_ResCurr_freq,
                    plan->h_ResCurr_freq,
                    RL_num);

                // 求取乘积的实部
                ComplexDotMulReal(plan->h_Temp,
                    plan->h_ResCurr_freq,
                    plan->h_ResTrail_freq,
                    RL_num);

                if(in_rank != min_rank_RL)
                {
                    MPI_Recv(buf, 2, MPI_FLOAT, in_rank - 1, STEP_FEN, mycomm, &status_recv);
                    *plan->h_SumFenzi = buf[0];
                    *plan->h_SumFenmu = buf[1];


                    for (uint ir = 0; ir < RL_num; ir++)
                    {
                        *plan->h_SumFenzi += plan->h_Temp[ir];
                    }

                    // 计算步长公式中的分母
                    memset((void *)plan->h_Temp,	0,	sizeof(float) * RL_num);

                    // 求取乘积的实部
                    CalSquareAbs(plan->h_Temp,
                        plan->h_ResTrail_freq,
                        RL_num);

                    for (uint ir = 0; ir < RL_num; ir++)
                    {
                        *plan->h_SumFenmu += plan->h_Temp[ir];
                    }

                    buf[0] = *plan->h_SumFenzi;
                    buf[1] = *plan->h_SumFenmu;

                    if(in_rank != max_rank_RL)
                    {
                        MPI_Send(buf, 2, MPI_FLOAT, in_rank + 1, STEP_FEN, mycomm);
                    }
                    else
                    {
                        ip->Alpha = (*plan->h_SumFenzi / *plan->h_SumFenmu) * e;
                    }
                }
                else
                {
                    for (uint ir = 0; ir < RL_num; ir++)
                    {
                        *plan->h_SumFenzi += plan->h_Temp[ir];
                    }

                    // 计算步长公式中的分母
                    memset((void *)plan->h_Temp,	0,	sizeof(float) * RL_num);

                    // 求取乘积的实部
                    CalSquareAbs(plan->h_Temp,
                        plan->h_ResTrail_freq,
                        RL_num);

                    for (uint ir = 0; ir < RL_num; ir++)
                    {
                        *plan->h_SumFenmu += plan->h_Temp[ir];
                    }

                    buf[0] = *plan->h_SumFenzi;
                    buf[1] = *plan->h_SumFenmu;
                    if(in_rank != max_rank_RL)
                        MPI_Send(buf, 2, MPI_FLOAT, in_rank + 1, STEP_FEN, mycomm);
                    else
                        ip->Alpha = (*plan->h_SumFenzi / *plan->h_SumFenmu) * e;
                }
            }

            MPI_Barrier(mycomm);
            MPI_Bcast(&ip->Alpha, 1, MPI_FLOAT, max_rank_RL, mycomm);
            if(rank == ROOT_ID)
                cout << "alpha=" << ip->Alpha << endl;
        }
    //}//shotN

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ip->Alpha, 1, MPI_FLOAT, ROOT_ID, MPI_COMM_WORLD);


//    if(RL_num)
//    {
//        if(in_rank != min_rank_RL)
//        {
//            MPI_Recv(buf, 2, MPI_FLOAT, in_rank - 1, STEP_FEN, mycomm, &status_recv);
//            *plan->h_SumFenzi = buf[0];
//            *plan->h_SumFenmu = buf[1];
//            for (uint m = 0; m < Pa.Nt * RL_num; m++)
//            {
//                *plan->h_SumFenzi += plan->h_SumResCurr[m] * plan->h_SumResTrial[m];
//                *plan->h_SumFenmu += plan->h_SumResTrial[m] * plan->h_SumResTrial[m];
//            }
//            buf[0] = *plan->h_SumFenzi;
//            buf[1] = *plan->h_SumFenmu;

//            if(in_rank != max_rank_RL)
//            {
//                MPI_Send(buf, 2, MPI_FLOAT, in_rank + 1, STEP_FEN, mycomm);
//            }
//            else
//            {
//                ip->Alpha = (*plan->h_SumFenzi / *plan->h_SumFenmu) * e;
//            }
//        }
//        else
//        {
//            for (uint m = 0; m < Pa.Nt * RL_num; m++)
//            {
//                *plan->h_SumFenzi += plan->h_SumResCurr[m] * plan->h_SumResTrial[m];
//                *plan->h_SumFenmu += plan->h_SumResTrial[m] * plan->h_SumResTrial[m];
//            }
//            buf[0] = *plan->h_SumFenzi;
//            buf[1] = *plan->h_SumFenmu;
//            if(in_rank != max_rank_RL)
//                MPI_Send(buf, 2, MPI_FLOAT, in_rank + 1, STEP_FEN, mycomm);
//            else
//                ip->Alpha = (*plan->h_SumFenzi / *plan->h_SumFenmu) * e;
//        }
//    }




    //ip->Alpha = *plan->h_SumFenzi / *plan->h_SumFenmu;
}


/*------------------------------------------------------------------------
function: PreProcess

Information:
下次迭代前的预处理

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
------------------------------------------------------------------------*/
void PreProcess(AFDP3D Pa,
                IP *ip,
                CPUVs *plan,
                const Partition& pt,
                const MPI_Comm& mycomm)
{
    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint nny = Pa.Ny + 2 * Pa.PMLy;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    int interiorLength_x = pt.getinteriorLength_x();
    int interiorLength_y = pt.getinteriorLength_y();
    int interiorLength_z = pt.getinteriorLength_z();

    if(interiorLength_x == 0 || interiorLength_y == 0 || interiorLength_z == 0)
        return;

    int interior_min_x = pt.getinteriormin_x();
    int interior_min_y = pt.getinteriormin_y();
    int interior_min_z = pt.getinteriormin_z();

    H_Border temph_Vp = pt.geth_Vp();

    for(int iz = 0; iz < interiorLength_z; ++iz)
    {
        for(int iy = 0; iy < interiorLength_y; ++iy)
        {
            memcpy(plan->h_Grad + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x, ip->GradVp + iz * interiorLength_y * interiorLength_x + iy * interiorLength_x, sizeof(float) * interiorLength_x);
        }
    }

    for(int iz = 0; iz < block_z; ++iz)
    {
        for(int iy = 0; iy < block_y; ++iy)
        {
            memcpy(plan->h_Vp + (iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + temph_Vp.leftborder, ip->CurrVp + iz * block_y * block_x + iy * block_x, sizeof(float) * block_x);
        }
    }


    // 更新速度
    UpdateVpInnerGrid(Pa, plan->h_Vp, plan->h_Grad, ip->Alpha, pt);

    vector<uint> trans_h_Vp = pt.gettrans_h_Vp();
    auto begin = trans_h_Vp.begin();
    if(*(begin + TOP))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, BOTTOM_TO_TOP, Pa, mycomm);
    }
    if(*(begin + LEFT))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, RIGHT_TO_LEFT, Pa, mycomm);
    }
    if(*(begin + BOTTOM))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, TOP_TO_BOTTOM, Pa, mycomm);
    }
    if(*(begin + RIGHT))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, LEFT_TO_RIGHT, Pa, mycomm);
    }

    dataGather(plan->h_Vp, pt, STEP_VP, mycomm);

    UpdateVpPML(Pa, plan->h_Vp, plan->h_Grad, ip->Alpha, pt);

    for(int iz = 0; iz < block_z; ++iz)
    {
        for(int iy = 0; iy < block_y; ++iy)
        {
            memcpy(ip->CurrVp + iz * block_y * block_x + iy * block_x, plan->h_Vp + (iz + temph_Vp.topborder) * temph_Vp.length_y * temph_Vp.length_x + (iy + temph_Vp.backborder) * temph_Vp.length_x + temph_Vp.leftborder, sizeof(float) * block_x);
        }
    }


    ip->Alpha = 0.0f;
}
