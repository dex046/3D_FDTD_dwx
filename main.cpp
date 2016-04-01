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

int main(int argc, char ** argv)
{
    ofstream fout("result.txt");
    struct tm *ptr;
    time_t lt;
    lt = time(NULL);
    fout << (ctime(&lt)) << endl;
    fout << "******************************************" << endl;

    int rank, p_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);

    uint cpu_x = 1, cpu_y = 8, cpu_z = 1;
    uint shot_x = 2, shot_y = 1, shot_z = 2;



    // 有限差分参数
    AFDP3D *Pa;
    Pa = new AFDP3D[1];
    memset((void *)Pa, 0, sizeof(AFDP3D));

    Pa->dt = 0.001f;
    Pa->dx = 10.0f;
    Pa->dy = 10.0f;
    Pa->dz = 10.0f;
    Pa->f0 = 10.0f;
    Pa->Nt = 1000;////////1000
    Pa->Nx = 200;
    Pa->Ny = 8;
    Pa->Nz = 100;
    Pa->PMLx = 50;
    Pa->PMLy = 50;
    Pa->PMLz = 50;
    uint nnx = Pa->Nx + 2 * Pa->PMLx;
    uint nny = Pa->Ny + 2 * Pa->PMLy;
    uint nnz = Pa->Nz + 2 * Pa->PMLz;


    H_Border temph_U(3, 3, 4, 4, 4, 3);//top left bottom right front back
    H_Border temph_VW(4, 4, 3, 3, 3, 4);

    // 反演参数
    IP *ip;
    ip = new IP[1];
    memset((void *)ip, 0, sizeof(IP));
    ip->ShotN = 1;
    ip->IterN = 1;
    ip->FreqN = 6;
    ip->Alpha = 0.0f;

    ip->St = new Shot[ip->ShotN];
    memset((void *)ip->St, 0, ip->ShotN * sizeof(Shot));

    // 炮信息
    for (uint is = 0; is < ip->ShotN; is++)
    {
        ip->St[is].rn = 200;
        ip->St[is].s.Sz = Pa->PMLz + 2;
        ip->St[is].s.Sy = Pa->PMLy + (uint)(Pa->Ny / 2);
        ip->St[is].s.Sx = Pa->PMLx + (uint)(Pa->Nx / 2);

//        ip->St[is].re = new RL[ip->St[is].rn];
//        memset((void *)ip->St[is].re, 0, sizeof(RL) * ip->St[is].rn);

//        // 检波器位置
//        for (uint ir = 0; ir < ip->St[is].rn; ir++)
//        {
//            ip->St[is].re[ir].Rz = Pa->PMLz + 2;
//            ip->St[is].re[ir].Ry = Pa->PMLy + (uint)(Pa->Ny / 2);
//            ip->St[is].re[ir].Rx = Pa->PMLx + ir;
//        }
    }



    int group_in_size = p_size / ip->ShotN;
    int *ranks = new int[group_in_size];
    int temp = rank % ip->ShotN;
    for(int i = 0; i < group_in_size; ++i)
    {
        ranks[i] = temp + i * ip->ShotN;
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(rank == ROOT_ID)
//    cout << "now1" <<endl;
    MPI_Group comm_group, newgroup;
    MPI_Comm mycomm;

    MPI_Comm_group(MPI_COMM_WORLD, &comm_group);

    MPI_Group_incl(comm_group, group_in_size, ranks, &newgroup);
    MPI_Comm_create(MPI_COMM_WORLD, newgroup, &mycomm);

    Partition pt(Pa, ip, nnx, nny, nnz, cpu_x, cpu_y, cpu_z, shot_x, shot_y, shot_z, temph_U, temph_VW, 8, rank, p_size, rank / ip->ShotN, group_in_size);
//H_Border temph_Vp = pt.geth_Vp();
//cout << temph_Vp.topborder << " " << temph_Vp.leftborder << " " << temph_Vp.bottomborder << " " << temph_Vp.rightborder << " " << temph_Vp.frontborder << " " << temph_Vp.backborder << endl;
//    if(1)
//    {
//        cout << pt.get_in_rank() << " " << pt.in_isfirstblock_z() << " " << pt.in_islastblock_z() << endl;
//    }
    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_y = pt.getinteriorLength_y();
    uint interiorlength_z = pt.getinteriorLength_z();

//    if(rank == 3)
//    {
//        cout << pt.geth_Vp().backborder << endl;
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    //if(rank == 0)
//    cout << interiorlength_y << endl;

    uint RL_num = pt.getRL_num();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmax_x = pt.getindexmax_x();

    uint indexmin_y = pt.getindexmin_y();
    uint indexmax_y = pt.getindexmax_y();

    uint indexmin_z = pt.getindexmin_z();
    uint indexmax_z = pt.getindexmax_z();

    try
    {
        ip->TrueVp = new float[block_z * block_y * block_x];
        memset((void *)ip->TrueVp, 0, sizeof(float) * block_z * block_y * block_x);
        ip->CurrVp = new float[block_z * block_y * block_x];
        memset((void *)ip->CurrVp, 0, sizeof(float) * block_z * block_y * block_x);

        ip->GradVp = new float[interiorlength_z * interiorlength_y * interiorlength_x];
        memset((void *)ip->GradVp, 0, sizeof(float) * interiorlength_z * interiorlength_y * interiorlength_x);
        ip->ObjIter = new float[ip->IterN];
        memset((void *)ip->ObjIter, 0, sizeof(float) * ip->IterN);
    }
    catch(const std::bad_alloc& e)
    {
        cout << "Allocation failed: " << e.what() << endl;
    }


    // 读入真实速度和初始速度
    char TrueVp[] = "TrueVp-3D.sgy";
    char InitVp[] = "InitVp-3D.sgy";

    ReadData(TrueVp, ip->TrueVp, pt, 0);
    ReadData(InitVp, ip->CurrVp, pt, 0);

//    if(rank == 0)
//    {

//        ofstream fout("True0.txt");
//        ofstream fout1("CurrVp0.txt");
//        for(int iz = 0; iz < block_z; ++iz)
//        {
//            for(int iy = 0; iy < block_y; ++iy)
//            {
//                for(int ix = 0; ix < block_x; ++ix)
//                {
//                    fout << *(ip->TrueVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                    fout1 << *(ip->CurrVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                }
//            }
//        }
//        fout.flush();
//        fout.close();

//        fout1.flush();
//        fout1.close();
//    }

//    if(rank == 1)
//    {
//        ofstream fout("True1.txt");
//        ofstream fout1("CurrVp1.txt");
//        for(int iz = 0; iz < block_z; ++iz)
//        {
//            for(int iy = 0; iy < block_y; ++iy)
//            {
//                for(int ix = 0; ix < block_x; ++ix)
//                {
//                    fout << *(ip->TrueVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                    fout1 << *(ip->CurrVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                }
//            }
//        }
//        fout.flush();
//        fout.close();

//        fout1.flush();
//        fout1.close();
//    }

//    if(rank == 2)
//    {
//        ofstream fout("True2.txt");
//        ofstream fout1("CurrVp2.txt");
//        for(int iz = 0; iz < block_z; ++iz)
//        {
//            for(int iy = 0; iy < block_y; ++iy)
//            {
//                for(int ix = 0; ix < block_x; ++ix)
//                {
//                    fout << *(ip->TrueVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                    fout1 << *(ip->CurrVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                }
//            }
//        }
//        fout.flush();
//        fout.close();

//        fout1.flush();
//        fout1.close();
//    }

//    if(rank == 3)
//    {/*cout << block_x << " " << block_y << " " << block_z << endl;*/
//        ofstream fout("True3.txt");
//        ofstream fout1("CurrVp3.txt");
//        for(int iz = 0; iz < block_z; ++iz)
//        {
//            for(int iy = 0; iy < block_y; ++iy)
//            {
//                for(int ix = 0; ix < block_x; ++ix)
//                {
//                    fout << *(ip->TrueVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                    fout1 << *(ip->CurrVp + iz * block_y * block_x + iy * block_x + ix) << " ";
//                }
//            }
//        }
//        fout.flush();
//        fout.close();

//        fout1.flush();
//        fout1.close();
//    }


//    int *temp_iz = new int[1];
//    for (uint iz = 0; iz < block_z; iz++)
//    {
//        for (uint iy = 0; iy < block_y; iy++)
//        {
//            for (uint ix = 0; ix < block_x; ix++)
//            {
//                *temp_iz = iz + indexmin_z;
//                if (*temp_iz < Pa->PMLz + (uint)(Pa->Nz / 2))
//                {
//                    ip->TrueVp[iz * block_y * block_x + iy * block_x + ix] = powf(2000.0f, 2.0f);
//                }
//                else
//                {
//                    ip->TrueVp[iz * block_y * block_x + iy * block_x + ix] = powf(4000.0f, 2.0f);
//                }

//                if (*temp_iz < Pa->PMLz + 5)
//                {
//                    ip->CurrVp[iz * block_y * block_x + iy * block_x + ix] = powf(2000.0f, 2.0f);
//                }
//                else if (*temp_iz <= Pa->PMLz + Pa->Nz)
//                {
//                    ip->CurrVp[iz * block_y * block_x + iy * block_x + ix] = powf((2000.0f + 2000.0f * (*temp_iz - Pa->PMLz - 5) / (Pa->Nz - 5)), 2.0f);
//                }
//                else
//                {
//                    ip->CurrVp[iz * block_y * block_x + iy * block_x + ix] = powf(4000.0f, 2.0f);
//                }
//            }
//        }
//    }
//    delete temp_iz;

    // 全局变量
    CPUVs *plan;
    plan = new CPUVs[1];
    memset((void *)plan, 0, sizeof(CPUVs));

    float *sgs_t = NULL, *sgs_c = NULL, *sgs_r = NULL;

    if(RL_num)
    {
        try
        {
            sgs_t = new float[ip->ShotN * Pa->Nt * RL_num];
            memset((void *)sgs_t, 0, sizeof(float) * ip->ShotN * Pa->Nt * RL_num);

            sgs_c = new float[ip->ShotN * Pa->Nt * RL_num];
            memset((void *)sgs_c, 0, sizeof(float) * ip->ShotN * Pa->Nt * RL_num);

            sgs_r = new float[ip->ShotN * Pa->Nt * RL_num];
            memset((void *)sgs_r, 0, sizeof(float) * ip->ShotN * Pa->Nt * RL_num);
        }
        catch(const std::bad_alloc& e)
        {
            cout << "Allocation failed: " << e.what() << endl;
        }
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    if(rank == ROOT_ID)
//    cout << "now3" <<endl;


    // 给全局变量开辟空间
    MallocVariables(*Pa, ip, plan, pt);

    // 求取NPML的参数
    GenerateNPML(*Pa, plan, pt);


//    MPI_Barrier(MPI_COMM_WORLD);
//    if(rank == ROOT_ID)
//    cout << "now5" <<endl;


    if(rank == ROOT_ID)
    {
        cout << "********************************************************************************************" << endl;
        cout << "(" << shot_x << " * " << shot_y << " * " << shot_z << ") * " << ip->ShotN << endl;
        cout << "********************************************************************************************" << endl;
        cout << "Doing 3D Hybrid Full Waveform Inversion ..." << endl;
        cout << "Time domain Forward modeling and frequency domain inversion" << endl;
        cout << "Parameters of Inversion are as follows:" << endl;
        cout << "\tNx = " << Pa->Nx << endl;
        cout << "\tNy = " << Pa->Ny << endl;
        cout << "\tNz = " << Pa->Nz << endl;
        cout << "\tdx = " << Pa->dx << "m" << endl;
        cout << "\tdy = " << Pa->dy << "m" << endl;
        cout << "\tdz = " << Pa->dz << "m" << endl;
        cout << "\tNt = " << Pa->Nt << endl;
        cout << "\tdt = " << Pa->dt << "s" << endl;
        cout << "\tNpml = " << Pa->PMLx << endl;
    #ifndef _FROM_TXT_
        cout << "\tf0 = " << Pa->f0 << endl;
    #endif
        cout << "\tNshot = " << ip->ShotN << endl;
        cout << "\tIteration number = " << ip->IterN << endl;
        cout << "\tUsing " << ip->FreqN << " frequencies" << endl;
        cout << "\tf_min = " << (StartFreq / 2.0f / PI) << "Hz" << endl;
        cout << "\tf_int = " << (InterFreq / 2.0f / PI) << "Hz" << endl;
        cout << "\tf_max = " << ((StartFreq + (ip->FreqN - 1) * InterFreq) / (2.0f * PI)) << "Hz" << endl;
        //cout << "---------------------------------------------------------" << endl;
    #ifdef _FROM_SGY_
        cout << "\tReading the wavelet from sgy" << endl;
    #endif
        cout << "********************************************************************************************" << endl;

        fout << "********************************************************************************************" << endl;
        fout << "(" << shot_x << " * " << shot_y << " * " << shot_z << ") * " << ip->ShotN << endl;
        fout << "********************************************************************************************" << endl;
        fout << "Doing 3D Hybrid Full Waveform Inversion ..." << endl;
        fout << "Time domain Forward modeling and frequency domain inversion" << endl;
        fout << "Parameters of Inversion are as follows:" << endl;
        fout << "\tNx = " << Pa->Nx << endl;
        fout << "\tNy = " << Pa->Ny << endl;
        fout << "\tNz = " << Pa->Nz << endl;
        fout << "\tdx = " << Pa->dx << "m" << endl;
        fout << "\tdy = " << Pa->dy << "m" << endl;
        fout << "\tdz = " << Pa->dz << "m" << endl;
        fout << "\tNt = " << Pa->Nt << endl;
        fout << "\tdt = " << Pa->dt << "s" << endl;
        fout << "\tNpml = " << Pa->PMLx << endl;
    #ifndef _FROM_TXT_
        fout << "\tf0 = " << Pa->f0 << endl;
    #endif
        fout << "\tNshot = " << ip->ShotN << endl;
        fout << "\tIteration number = " << ip->IterN << endl;
        fout << "\tUsing " << ip->FreqN << " frequencies" << endl;
        fout << "\tf_min = " << (StartFreq / 2.0f / PI) << "Hz" << endl;
        fout << "\tf_int = " << (InterFreq / 2.0f / PI) << "Hz" << endl;
        fout << "\tf_max = " << ((StartFreq + (ip->FreqN - 1) * InterFreq) / (2.0f * PI)) << "Hz" << endl;
        //cout << "---------------------------------------------------------" << endl;
    #ifdef _FROM_SGY_
        fout << "\tReading the wavelet from sgy" << endl;
    #endif
        fout << "********************************************************************************************" << endl;
    }

    clock_t begin, duration;

    if(rank == ROOT_ID)
    {
        // 求取观测波场
//        cout << "\tCalculating the observed data..." << endl;

        fout << "\tCalculating the observed data..." << endl;
    }

    begin = clock();
    CalTrueWF(*Pa, ip, plan, sgs_t, pt, mycomm);
    duration = clock() - begin;

    if(rank == ROOT_ID)
    {
//        cout << "\tCalculating the observed data used:\t" << duration / CLOCKS_PER_SEC << "s" << endl;

        fout << "\tCalculating the observed data used:\t" << duration / CLOCKS_PER_SEC << "s" << endl;
    }


    // 开始迭代
    float e = 12.0f;

    for (uint It = 0; It < ip->IterN; It++)
    {
        if(rank == ROOT_ID)
        {
            // 求取梯度
            cout << "\n\tDoing the " << It << "th iteration" << endl;
            cout << "\tCalculating the Gradient..." << endl;

            fout << "\n\tDoing the " << It << "th iteration" << endl;
            fout << "\tCalculating the Gradient..." << endl;
        }

        begin = clock();
        CalGrad(*Pa, ip, plan, sgs_t, sgs_c, sgs_r, It, pt, mycomm);

        // 梯度后处理
        PostProcessGrad(*Pa, ip->GradVp, plan->h_Vp, pt);

        // 求取步长
        CalStepLength(*Pa, ip, plan, sgs_t, sgs_c, e, pt, mycomm);
        duration = clock() - begin;

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == ROOT_ID)
        {
            cout << "\tObjective function value:\t" << ip->ObjIter[It] << endl;
            cout << "\tStep length:\t" << ip->Alpha << endl;
            cout << "\tThe " << It << "th iteration used " << duration / CLOCKS_PER_SEC << "s" << endl;

            fout << "\tObjective function value:\t" << ip->ObjIter[It] << endl;
            fout << "\tStep length:\t" << ip->Alpha << endl;
            fout << "\tThe " << It << "th iteration used " << duration / CLOCKS_PER_SEC << "s" << endl;
        }

        //cout << rank << endl;
        // 下一步迭代预处理
        PreProcess(*Pa, ip, plan, pt, mycomm);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == ROOT_ID)
    {
//        cout << "\tWriting data to .sgy" << endl;

        fout << "\tWriting data to .sgy" << endl;
    }

    const char * const TrueSg = "TrueSG.sgy";
    const char * const GradVp = "GradientVp.sgy";
    const char * const InvertedVp = "InvertedVp.sgy";

    fopen(TrueSg, "wb");
    fopen(GradVp, "wb");
    fopen(InvertedVp, "wb");

//    write_sgs_t_Data(TrueSg, (usht)Pa->Nt, (usht)ip->St[0].rn, (usht)(Pa->dt * 1000000), sgs_t, pt, *Pa, 1);
//    WriteData(GradVp, Pa->Nz, Pa->Nx * Pa->Ny, Pa->dz * 1000, ip->GradVp, 0, pt, *Pa, WRITE_INTER);
//    WriteData(InvertedVp, nnz, nnx * nny, Pa->dz * 1000, ip->CurrVp, 0, pt, *Pa, WRITE_ALL);


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;

}
