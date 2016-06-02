/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/

#include "DataTran.h"

void dataTran_temp()
{

}

void copydatatobuf(float *data, float *buf, const Partition& pt, uint transportlen_side, int tag, int flag)
{
    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint topborder = 0;
    uint leftborder = 0;
    //uint bottomborder = 0;
    uint rightborder = 0;

    uint backborder = 0;
    uint frontborder = 0;

    uint length_x = block_x;
    uint length_y = block_y;
    //uint length_z = block_z;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        //bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;
        backborder = temph_U.backborder;
        frontborder = temph_U.frontborder;

        length_x = leftborder + block_x + rightborder;
        length_y = backborder + block_y + frontborder;
        //length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_V || tag == STEP_W || tag == STEP_K)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        //bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;
        backborder = temph_VW.backborder;;
        frontborder = temph_VW.frontborder;

        length_x = leftborder + block_x + rightborder;
        length_y = backborder + block_y + frontborder;
        //length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        //bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;
        backborder = temph_Vp.backborder;
        frontborder = temph_Vp.frontborder;

        length_x = leftborder + block_x + rightborder;
        length_y = backborder + block_y + frontborder;
        //length_z = topborder + block_z + bottomborder;
    }
    else
    {

    }

    if(flag == RIGHT_TO_LEFT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * transportlen_side + j * transportlen_side, data + (i + topborder) * length_y * length_x + (j + backborder) * length_x + leftborder, sizeof(float) * transportlen_side);
            }
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * transportlen_side + j * transportlen_side, data + (i + topborder) * length_y * length_x + (j + backborder) * length_x + leftborder + block_x - transportlen_side, sizeof(float) * transportlen_side);
            }
        }
    }
    if(flag == TOP_TO_BOTTOM)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * block_x + j * block_x, data + (block_z + topborder + i - transportlen_side) * length_y * length_x + (j + backborder) * length_x + leftborder, sizeof(float) * block_x);
            }//cout << 1324 <<endl;
        }
//        if(pt.getrank() == 1)
//        {
//            for(int i = 0; i < transportlen_side * block_y * block_x; ++i)
//            {
//                if(*(buf + i))
//                    cout << *(buf + i) << " ";
//            }
//            cout << endl;
//        }
    }
    if(flag == BOTTOM_TO_TOP)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * block_x + j * block_x, data + (topborder + i) * length_y * length_x + (j + backborder) * length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
    if(flag == FRONT_TO_BACK)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < transportlen_side; ++j)
            {
                memcpy(buf + i * transportlen_side * block_x + j * block_x, data + (topborder + i) * length_y * length_x + (j + backborder) * length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
    if(flag == BACK_TO_FRONT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < transportlen_side; ++j)
            {
                memcpy(buf + i * transportlen_side * block_x + j * block_x, data + (topborder + i) * length_y * length_x + (backborder + block_y + j - transportlen_side) * length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
}
void copydatatobuf_Vp(float *data, float *buf, const Partition& pt, uint transportlen_side, int flag, const AFDP3D &Pa)
{
    int gap_x = 0;
    int gap_y = 0;
    int gap_z = 0;

    H_Border temph_Vp = pt.geth_Vp();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_y = pt.getindexmin_y();
    uint indexmin_z = pt.getindexmin_z();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint topborder = temph_Vp.topborder;
    uint leftborder = temph_Vp.leftborder;
//    uint bottomborder = temph_Vp.bottomborder;
//    uint rightborder = temph_Vp.rightborder;
    uint backborder = temph_Vp.backborder;
//    uint frontborder = temph_Vp.frontborder;

    if(flag == RIGHT_TO_LEFT)
    {
        gap_x = Pa.PMLx - indexmin_x;
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * transportlen_side + j * transportlen_side, data + (i + topborder) * temph_Vp.length_y * temph_Vp.length_x + (j + backborder) * temph_Vp.length_x + leftborder + gap_x, sizeof(float) * transportlen_side);
            }
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {
        gap_x = Pa.PMLx + Pa.Nx - 1 - indexmin_x;
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * transportlen_side + j * transportlen_side, data + (i + topborder) * temph_Vp.length_y * temph_Vp.length_x + (j + backborder) * temph_Vp.length_x + leftborder + gap_x, sizeof(float) * transportlen_side);
            }
        }
    }
    if(flag == TOP_TO_BOTTOM)
    {
        gap_z = Pa.PMLz + Pa.Nz - 1 - indexmin_z;
        for(uint i = 0; i < transportlen_side; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * block_x + j * block_x, data + (topborder + i + gap_z) * temph_Vp.length_y * temph_Vp.length_x + (backborder + j) * temph_Vp.length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
    if(flag == BOTTOM_TO_TOP)
    {
        gap_z = Pa.PMLz - indexmin_z;
        for(uint i = 0; i < transportlen_side; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(buf + i * block_y * block_x + j * block_x, data + (topborder + i + gap_z) * temph_Vp.length_y * temph_Vp.length_x + (backborder + j) * temph_Vp.length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
    if(flag == FRONT_TO_BACK)
    {
        gap_y = Pa.PMLy - indexmin_y;
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < transportlen_side; ++j)
            {
                memcpy(buf + i * transportlen_side * block_x + j * block_x, data + (topborder + i) * temph_Vp.length_y * temph_Vp.length_x + (backborder + j + gap_y) * temph_Vp.length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
    if(flag == BACK_TO_FRONT)
    {
        gap_y = Pa.PMLy + Pa.Ny - 1 - indexmin_y;
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < transportlen_side; ++j)
            {
                memcpy(buf + i * transportlen_side * block_x + j * block_x, data + (topborder + i) * temph_Vp.length_y * temph_Vp.length_x + (backborder + j + gap_y) * temph_Vp.length_x + leftborder, sizeof(float) * block_x);
            }
        }
    }
}



void copybuftodata(float *buf, float *data, const Partition& pt, uint transportlen_side, int tag, int flag)
{
    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint topborder = 0;
    uint leftborder = 0;
    //uint bottomborder = 0;
    uint rightborder = 0;
    uint frontborder= 0;
    uint backborder = 0;

    uint length_x = block_x;
    uint length_y = block_y;
    //uint length_z = block_z;
    //transportlen_side = 0;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        //bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;
        backborder = temph_U.backborder;
        frontborder = temph_U.frontborder;

        length_x = leftborder + block_x + rightborder;
        length_y = backborder + block_y + frontborder;
        //length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_V || tag == STEP_W || tag == STEP_K)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        //bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;
        backborder = temph_VW.backborder;;
        frontborder = temph_VW.frontborder;

        length_x = leftborder + block_x + rightborder;
        length_y = backborder + block_y + frontborder;
        //length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        //bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;
        backborder = temph_Vp.backborder;
        frontborder = temph_Vp.frontborder;

        length_x = leftborder + block_x + rightborder;
        length_y = backborder + block_y + frontborder;
        //length_z = topborder + block_z + bottomborder;
    }
    else
    {

    }

    if(flag == TOP_TO_BOTTOM)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(data + i * length_y * length_x + (j + backborder) * length_x + leftborder, buf + i * block_y * block_x + j * block_x, sizeof(float) * block_x);
            }
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {/*cout << block_z * transportlen_side << endl;
        for(int i = 0; i < block_z * transportlen_side; ++i)
            cout << *(buf + i);*/
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(data + (i + topborder) * length_y * length_x + (j + backborder) * length_x, buf + i * block_y * transportlen_side + j * transportlen_side, sizeof(float) * transportlen_side);
            }
        }

    }
    if(flag == BOTTOM_TO_TOP)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(data + (i + topborder + block_z) * length_y * length_x + (j + backborder) * length_x + leftborder, buf + i * block_y * block_x + j * block_x, sizeof(float) * block_x);
            }
        }
    }
    if(flag == RIGHT_TO_LEFT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < block_y; ++j)
            {
                memcpy(data + (i + topborder) * length_y * length_x + (j + backborder) * length_x + leftborder + block_x, buf + i * transportlen_side * block_y + j * transportlen_side, sizeof(float) * transportlen_side);
            }
        }
    }
    if(flag == FRONT_TO_BACK)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < transportlen_side; ++j)
            {
                memcpy(data + (i + topborder) * length_y * length_x + (j + backborder + block_y) * length_x + leftborder, buf + i * block_x * transportlen_side + j * block_x, sizeof(float) * block_x);
            }
        }
    }
    if(flag == BACK_TO_FRONT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            for(uint j = 0; j < transportlen_side; ++j)
            {
                memcpy(data + (i + topborder) * length_y * length_x + j * length_x + leftborder, buf + i * block_x * transportlen_side + j * block_x, sizeof(float) * block_x);
            }
        }
    }
}
void dataGather(float *data, const Partition& pt, int tag, const MPI_Comm& mycomm)
{
    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint topborder = 0;
    uint leftborder = 0;
    uint bottomborder = 0;
    uint rightborder = 0;
    uint frontborder = 0;
    uint backborder = 0;

    uint transportlen_side = 0;
    uint length = 0;

//    uint length_x = block_x;
//    uint length_y = block_y;
//    uint length_z = block_z;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;
        backborder = temph_U.backborder;
        frontborder = temph_U.frontborder;

//        length_x = leftborder + block_x + rightborder;
//        length_y = backborder + block_y + frontborder;
//        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_V || tag == STEP_W || tag == STEP_K)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;
        backborder = temph_VW.backborder;;
        frontborder = temph_VW.frontborder;

//        length_x = leftborder + block_x + rightborder;
//        length_y = backborder + block_y + frontborder;
//        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;
        backborder = temph_Vp.backborder;
        frontborder = temph_Vp.frontborder;

//        length_x = leftborder + block_x + rightborder;
//        length_y = backborder + block_y + frontborder;
//        length_z = topborder + block_z + bottomborder;
    }
    else
    {

    }

    MPI_Status status;
    //MPI_Request request;


    if(topborder && !pt.in_isfirstblock_z())
    {
        length = topborder * block_y * block_x;
        float *buf = new float[length];
        MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + TOP_TO_BOTTOM, mycomm, &status);

//        for(int i = 0; i < length; ++i)
//            cout << *(buf + i);

        transportlen_side = topborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, TOP_TO_BOTTOM);

        delete[] buf;
    }

    if(leftborder && !pt.in_isfirstblock_x())
    {

        length = block_z * block_y * leftborder;
        float *buf = new float[length];
        MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + LEFT_TO_RIGHT, mycomm, &status);

        transportlen_side = leftborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, LEFT_TO_RIGHT);

        delete[] buf;
    }

    if(bottomborder && !pt.in_islastblock_z())
    {
        length = bottomborder * block_y * block_x;
        float *buf = new float[length];
        MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BOTTOM_TO_TOP, mycomm, &status);
        transportlen_side = bottomborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, BOTTOM_TO_TOP);

        delete[] buf;
    }

    if(rightborder && !pt.in_islastblock_x())
    {
        length = block_z * block_y * rightborder;
        float *buf = new float[length];
        MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + RIGHT_TO_LEFT, mycomm, &status);
        transportlen_side = rightborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, RIGHT_TO_LEFT);

        delete[] buf;
    }

    if(frontborder && !pt.in_islastblock_y())
    {
        length = block_z * frontborder * block_x;
        float *buf = new float[length];
        MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + FRONT_TO_BACK, mycomm, &status);
        transportlen_side = frontborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, FRONT_TO_BACK);

        delete[] buf;
    }

//    if(pt.getrank() == 3)
//    cout << "EEEE" << backborder << endl;
    if(backborder && !pt.in_isfirstblock_y())
    {
        length = block_z * backborder * block_x;
        float *buf = new float[length];
        MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BACK_TO_FRONT, mycomm, &status);
        transportlen_side = backborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, BACK_TO_FRONT);

        delete[] buf;
    }
}


void dataTransport(float *data, const Partition& pt, int tag, int it, const MPI_Comm &mycomm)
{

    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    MPI_Status status;

    uint pos_x = pt.get_in_blockPosition_x();
    uint pos_y = pt.get_in_blockPosition_y();
    uint pos_z = pt.get_in_blockPosition_z();

//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << pt.getrank() << " " << pos_y << endl;

    //int rank = pt.getrank();
    int in_rank = pt.get_in_rank();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();

    uint sumBlock_x = pt.get_in_sumBlock_x();
    uint sumBlock_y = pt.get_in_sumBlock_y();
    //uint sumBlock_z = pt.get_in_sumBlock_z();

    uint topborder = 0;
    uint leftborder = 0;
    uint bottomborder = 0;
    uint rightborder = 0;
    uint frontborder = 0;
    uint backborder = 0;

    uint transportlen_side = 0;
    uint length = 0;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;
        frontborder = temph_U.frontborder;
        backborder = temph_U.backborder;
    }
    else if(tag == STEP_V || tag == STEP_W || tag == STEP_K)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;
        frontborder = temph_VW.frontborder;
        backborder = temph_VW.backborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;
        frontborder = temph_Vp.frontborder;
        backborder = temph_Vp.backborder;
    }
    else
    {

    }




    if(pos_z % 2)
    {
        if(!pt.in_islastblock_z())
        {
            //cout << pt.getrank() + pt.getsumBlock_x() << "www" << pt.getrank() << endl;
            uint transportlength_side = topborder;
            length = transportlength_side * block_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_side, tag, TOP_TO_BOTTOM);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + sumBlock_x * sumBlock_y, tag + TOP_TO_BOTTOM, mycomm);

            delete[] buf;
        }
        if(!pt.in_isfirstblock_z())
        {
            uint transportlength_side = bottomborder;
            length = transportlength_side * block_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_side, tag, BOTTOM_TO_TOP);

//            if(it == 50 && pt.getrank() == 2 && tag == STEP_U)
//                //cout << tag << endl;
//            for(int i = 0; i < length; ++i)
//                cout << *(buf + i) << " ";
            MPI_Send(buf, length, MPI_FLOAT, in_rank - sumBlock_x * sumBlock_y, tag + BOTTOM_TO_TOP, mycomm);
//cout << "SSS" << in_rank - sumBlock_x * sumBlock_y << endl;
            delete[] buf;
        }
        ///////////////////
        if(!pt.in_islastblock_z())
        {
            length = bottomborder * block_x * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BOTTOM_TO_TOP, mycomm, &status);
            transportlen_side = bottomborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, BOTTOM_TO_TOP);

            delete[] buf;
        }
        if(!pt.in_isfirstblock_z())
        {
            length = topborder * block_x * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + TOP_TO_BOTTOM, mycomm, &status);

    //        for(int i = 0; i < length; ++i)
    //            cout << *(buf + i);

            transportlen_side = topborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, TOP_TO_BOTTOM);

            delete[] buf;
        }
    }
    else
    {
        if(!pt.in_islastblock_z())
        {
            length = bottomborder * block_x * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BOTTOM_TO_TOP, mycomm, &status);
            transportlen_side = bottomborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, BOTTOM_TO_TOP);

            delete[] buf;

//            if(it == 50 && pt.getrank() == 0 && tag == STEP_U)
//                //cout << tag << endl;
//            for(int i = 0; i < length; ++i)
//                if(*(buf + i))
//                cout << i % block_x << " ";

        }
        if(!pt.in_isfirstblock_z())
        {
            length = topborder * block_x * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + TOP_TO_BOTTOM, mycomm, &status);

    //        for(int i = 0; i < length; ++i)
    //            cout << *(buf + i);

            transportlen_side = topborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, TOP_TO_BOTTOM);

            delete[] buf;
        }
        //////////////////
        if(!pt.in_islastblock_z())
        {
            //cout << pt.getrank() + pt.getsumBlock_x() << "www" << pt.getrank() << endl;
            uint transportlen_side = topborder;
            length = transportlen_side * block_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlen_side, tag, TOP_TO_BOTTOM);

            MPI_Send(buf, length, MPI_FLOAT, in_rank + sumBlock_x * sumBlock_y, tag + TOP_TO_BOTTOM, mycomm);


            delete[] buf;
        }
        if(!pt.in_isfirstblock_z())
        {
            uint transportlen_side = bottomborder;
            length = transportlen_side * block_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlen_side, tag, BOTTOM_TO_TOP);
            MPI_Send(buf, length, MPI_FLOAT, in_rank - sumBlock_x * sumBlock_y, tag + BOTTOM_TO_TOP, mycomm);

            delete[] buf;
        }
    }

    if(pos_y % 2)
    {
        if(!pt.in_islastblock_y())
        {
            //cout << pt.getrank() + pt.getsumBlock_x() << "www" << pt.getrank() << endl;
            uint transportlength_y = backborder;
            length = block_z * transportlength_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_y, tag, BACK_TO_FRONT);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + sumBlock_x, tag + BACK_TO_FRONT, mycomm);

            delete[] buf;
        }
        if(!pt.in_isfirstblock_y())
        {
            uint transportlength_y = frontborder;
            length = block_z * transportlength_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_y, tag, FRONT_TO_BACK);

//            if(it == 50 && pt.getrank() == 2 && tag == STEP_U)
//                //cout << tag << endl;
//            for(int i = 0; i < length; ++i)
//                cout << *(buf + i) << " ";
            MPI_Send(buf, length, MPI_FLOAT, in_rank - sumBlock_x, tag + FRONT_TO_BACK, mycomm);

            delete[] buf;
        }
        //////////////////
        if(!pt.in_islastblock_y())
        {
            length = block_z * frontborder * block_x;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + FRONT_TO_BACK, mycomm, &status);
            transportlen_side = frontborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, FRONT_TO_BACK);

            delete[] buf;
        }
        ///////////////////
        if(!pt.in_isfirstblock_y())
        {
            length = block_z * backborder * block_x;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BACK_TO_FRONT, mycomm, &status);

    //        for(int i = 0; i < length; ++i)
    //            cout << *(buf + i);

            transportlen_side = backborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, BACK_TO_FRONT);

            delete[] buf;
        }
    }
    else
    {
        if(!pt.in_islastblock_y())
        {
            length = block_z * frontborder * block_x;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + FRONT_TO_BACK, mycomm, &status);
            transportlen_side = frontborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, FRONT_TO_BACK);

            delete[] buf;

//            if(it == 50 && pt.getrank() == 0 && tag == STEP_U)
//                //cout << tag << endl;
//            for(int i = 0; i < length; ++i)
//                if(*(buf + i))
//                cout << i % block_x << " ";

        }
        if(!pt.in_isfirstblock_y())
        {
            length = block_z * backborder * block_x;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BACK_TO_FRONT, mycomm, &status);

    //        for(int i = 0; i < length; ++i)
    //            cout << *(buf + i);

            transportlen_side = backborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, BACK_TO_FRONT);

            delete[] buf;
        }
        //////////////////
        if(!pt.in_islastblock_y())
        {
            //cout << pt.getrank() + pt.getsumBlock_x() << "www" << pt.getrank() << endl;
            uint transportlength_y = backborder;
            length = block_z * transportlength_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_y, tag, BACK_TO_FRONT);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + sumBlock_x, tag + BACK_TO_FRONT, mycomm);

            delete[] buf;
        }
        if(!pt.in_isfirstblock_y())
        {
            uint transportlength_y = frontborder;
            length = block_z * transportlength_y * block_x;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_y, tag, FRONT_TO_BACK);
            MPI_Send(buf, length, MPI_FLOAT, in_rank - sumBlock_x, tag + FRONT_TO_BACK, mycomm);

            delete[] buf;
        }
    }

    if(pos_x % 2)
    {
        if(!pt.in_islastblock_x())
        {

            uint transportlength_x = leftborder;
            length = transportlength_x * block_z * block_y;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_x, tag, LEFT_TO_RIGHT);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + 1, tag + LEFT_TO_RIGHT, mycomm);

            delete[] buf;
        }
        if(!pt.in_isfirstblock_x())
        {

            uint transportlength_x = rightborder;
            length = transportlength_x * block_z * block_y;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_x, tag, RIGHT_TO_LEFT);
            MPI_Send(buf, length, MPI_FLOAT, in_rank - 1, tag + RIGHT_TO_LEFT, mycomm);
//cout << "RRR" << pt.get_in_rank() - 1 << endl;


            delete[] buf;
        }


        if(!pt.in_islastblock_x())
        {
            length = rightborder * block_z * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + RIGHT_TO_LEFT, mycomm, &status);
            transportlen_side = rightborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, RIGHT_TO_LEFT);

            delete[] buf;
        }
        if(!pt.in_isfirstblock_x())
        {//cout << pt.getrank() << endl;
            length = leftborder * block_z * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + LEFT_TO_RIGHT, mycomm, &status);
//cout << "www" << pt.get_in_rank() << endl;
//            if(pt.getrank() == 1)
//                for(int i = 0; i < length; ++i)
//                    cout << *(buf + i);

            transportlen_side = leftborder;
            //cout << transportlen_side << endl;
            copybuftodata(buf, data, pt, transportlen_side, tag, LEFT_TO_RIGHT);
//            for(int i = 0; i < length; ++i)
//                cout << *(data + i) << " ";

            delete[] buf;
        }
    }
    else
    {
        if(!pt.in_islastblock_x())
        {
            length = rightborder * block_z * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + RIGHT_TO_LEFT, mycomm, &status);



            transportlen_side = rightborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, RIGHT_TO_LEFT);

//            if(pt.getrank() == 0 && tag == STEP_U)
//            {
//                //cout << tag << endl << endl;
//                cout << "it=" << it << endl;
//                for(int i = 0; i < length; ++i)
//                {
//                    if(*(buf + i))
//                        cout << *(buf + i) << " ";
//                }
//                cout << endl;
//            }
//cout << "RRR" << pt.get_in_rank() + 1 << endl;
            delete[] buf;
        }
        if(!pt.in_isfirstblock_x())
        {
            length = leftborder * block_z * block_y;
            float *buf = new float[length];
            MPI_Recv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + LEFT_TO_RIGHT, mycomm, &status);

            transportlen_side = leftborder;
            copybuftodata(buf, data, pt, transportlen_side, tag, LEFT_TO_RIGHT);

            delete[] buf;
        }


        if(!pt.in_islastblock_x())
        {

            uint transportlength_x = leftborder;
            length = transportlength_x * block_y * block_z;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_x, tag, LEFT_TO_RIGHT);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + 1, tag + LEFT_TO_RIGHT, mycomm);
//cout << "SSS" << pt.get_in_rank() + 1 << endl;
            delete[] buf;
        }
        if(!pt.in_isfirstblock_x())
        {

            uint transportlength_x = rightborder;
            length = transportlength_x * block_y * block_z;
            float *buf = new float[length];
            copydatatobuf(data, buf, pt, transportlength_x, tag, RIGHT_TO_LEFT);//cout << "RRR" << pt.get_in_rank() << endl;
            MPI_Send(buf, length, MPI_FLOAT, in_rank - 1, tag + RIGHT_TO_LEFT, mycomm);

            delete[] buf;
        }
    }


//    if(!pt.islastblock_z())
//    {


//    }
//    if(!pt.islastblock_x())
//    {
//        //cout << "toright" << pt.getrank() << endl;
//        uint transportlength_x = 0;
//        switch(tag)
//        {
//        case STEP_U: transportlength_x = temph_U.leftborder;break;
//        case STEP_V: transportlength_x = temph_VW.leftborder;break;
//        case STEP_W: transportlength_x = temph_VW.leftborder;break;
//        };

//        uint length = transportlength_x * pt.getblockLength_z();
//        float *buf = new float[length];
//        copydatatobuf(data, buf, pt, transportlength_x, tag, LEFT_TO_RIGHT);
//        MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + 1, tag + LEFT_TO_RIGHT, MPI_COMM_WORLD);
//    }
//    if(!pt.isfirstblock_z())
//    {
//        //cout << "totop" << pt.getrank() << endl;

//    }
//    if(!pt.isfirstblock_x())
//    {
//        //cout << "toleft" << pt.getrank() << endl;
//        uint transportlength_x = 0;
//        switch(tag)
//        {
//        case STEP_U: transportlength_x = temph_U.rightborder;break;
//        case STEP_V: transportlength_x = temph_VW.rightborder;break;
//        case STEP_W: transportlength_x = temph_VW.rightborder;break;
//        };

//        uint length = transportlength_x * pt.getblockLength_z();
//        float *buf = new float[length];
//        copydatatobuf(data, buf, pt, transportlength_x, tag, RIGHT_TO_LEFT);
//        MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - 1, tag + RIGHT_TO_LEFT, MPI_COMM_WORLD);
//    }
}


void dataTransport_Vp(float *data, const Partition& pt, int tag, const AFDP3D &Pa, const MPI_Comm& mycomm)
{
    //uint rank = pt.getrank();
    uint in_rank = pt.get_in_rank();

    uint blockPosition_x = pt.get_in_blockPosition_x();
    uint blockPosition_y = pt.get_in_blockPosition_y();
    uint blockPosition_z = pt.get_in_blockPosition_z();

    uint sumBlock_x = pt.get_in_sumBlock_x();
    uint sumBlock_y = pt.get_in_sumBlock_y();
    uint sumBlock_z = pt.get_in_sumBlock_z();

    uint block_x = pt.getblockLength_x();
    uint block_y = pt.getblockLength_y();
    uint block_z = pt.getblockLength_z();
    //cout << rank << endl;
//    if(rank == 15)
//    {
//        cout << "ooooooooooo" << endl;
//    }

    int transportlength_side = 1;

    if(tag == RIGHT_TO_LEFT)
    {
        for(uint i = 1; i <= blockPosition_x; ++i)
        {
            uint length = transportlength_side * block_z * block_y;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_side, RIGHT_TO_LEFT, Pa);
            MPI_Send(buf, length, MPI_FLOAT, in_rank - i, STEP_VP + RIGHT_TO_LEFT, mycomm);

            delete[] buf;
        }
    }
    if(tag == LEFT_TO_RIGHT)
    {
        for(uint i = 1; i < sumBlock_x - blockPosition_x; ++i)
        {
            uint length = transportlength_side * block_z * block_y;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_side, LEFT_TO_RIGHT, Pa);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + i, STEP_VP + LEFT_TO_RIGHT, mycomm);

            delete[] buf;
        }
    }
    if(tag == BOTTOM_TO_TOP)
    {
        for(uint i = 1; i <= blockPosition_z ; ++i)
        {
            uint length = transportlength_side * block_x * block_y;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_side, BOTTOM_TO_TOP, Pa);
            MPI_Send(buf, length, MPI_FLOAT, in_rank - (i * sumBlock_x * sumBlock_y), STEP_VP + BOTTOM_TO_TOP, mycomm);

            delete[] buf;
        }
    }
    if(tag == TOP_TO_BOTTOM)
    {//cout << "wwwwww" << rank << endl;
        for(uint i = 1; i < sumBlock_z - blockPosition_z; ++i)
        {
            uint length = transportlength_side * block_x * block_y;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_side, TOP_TO_BOTTOM, Pa);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + (i * sumBlock_x * sumBlock_y), STEP_VP + TOP_TO_BOTTOM, mycomm);

            delete[] buf;
        }
    }
    if(tag == BACK_TO_FRONT)
    {//cout << "EEE" << in_rank + (1 * sumBlock_x) << endl;
        for(uint i = 1; i < sumBlock_y - blockPosition_y; ++i)
        {
            uint length = transportlength_side * block_x * block_z;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_side, BACK_TO_FRONT, Pa);
            MPI_Send(buf, length, MPI_FLOAT, in_rank + (i * sumBlock_x), STEP_VP + BACK_TO_FRONT, mycomm);

            delete[] buf;
        }
    }
    if(tag == FRONT_TO_BACK)
    {
        for(uint i = 1; i <= blockPosition_y; ++i)
        {
            uint length = transportlength_side * block_x * block_z;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_side, FRONT_TO_BACK, Pa);
            MPI_Send(buf, length, MPI_FLOAT, in_rank - (i * sumBlock_x), STEP_VP + FRONT_TO_BACK, mycomm);

            delete[] buf;
        }
    }
}

