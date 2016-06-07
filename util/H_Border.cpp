/******************************************
 * author:dwx
 ******************************************/
#include "Partition.h"
H_Border::H_Border()
{

}
H_Border::H_Border(uint topborder, uint leftborder, uint bottomborder, uint rightborder, uint frontborder, uint backborder)
{
    this->topborder = topborder;
    this->leftborder = leftborder;
    this->bottomborder = bottomborder;
    this->rightborder = rightborder;

    this->frontborder = frontborder;
    this->backborder = backborder;
}
H_Border::H_Border(uint length_x, uint legnth_z, uint topborder, uint leftborder, uint bottomborder, uint rightborder, uint frontborder, uint backborder)
{
    this->length_x = length_x;
    this->length_y = length_y;
    this->length_z = length_z;

    this->topborder = topborder;
    this->leftborder = leftborder;
    this->bottomborder = bottomborder;
    this->rightborder = rightborder;

    this->frontborder = frontborder;
    this->backborder = backborder;
}

