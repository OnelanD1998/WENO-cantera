#include "myField.H"
#include <fstream>
#include <string>
#include <sstream>
#include "dictionary.H"
#include "timeTerm.H"
#include "mpi.h"


#ifndef MESH_H
#define MESH_H

class mesh
:
public dictionary
{
public:
    int wenoR;
    int xNode;
    int yNode;
    int zNode;

    double xLeft;
    double xRight;
    double yLeft;
    double yRight;
    double zLeft;
    double zRight;

    double delX;
    double delY;
    double delZ;

    //-- for parallel
    int rank;
    int size;
    int len;

    double* xField;
    double* yField;
    double* zField;

    timeTerm  meshTime_;

public:
    mesh();

    ~mesh();
};

mesh::mesh( )
:
dictionary("meshDic"),
wenoR(this->lookUp<int >("wenoR")),
xNode(this->lookUp<int>("xNode")),
yNode(this->lookUp<int>("yNode")),
zNode(this->lookUp<int>("zNode")),
xLeft(this->lookUp<double>("xLeft")),
xRight(this->lookUp<double>("xRight")),
yLeft(this->lookUp<double>("yLeft")),
yRight(this->lookUp<double>("yRight")),
zLeft(this->lookUp<double>("zLeft")),
zRight(this->lookUp<double>("zRight"))
{
    char version[MPI_MAX_LIBRARY_VERSION_STRING];

    //MPI::Init();
    this->rank = MPI::COMM_WORLD.Get_rank();
    this->size = MPI::COMM_WORLD.Get_size();
    MPI_Get_library_version(version, & (this->len));

    xNode = xNode/this->size; //-- split mesh

    xField = new double[xNode + 2*wenoR];
    yField = new double[yNode + 2*wenoR];
    zField = new double[zNode + 2*wenoR];

    double subLength = (xRight-xLeft)/size;
    xLeft = this->rank * subLength + xLeft;
    xRight = this->rank * subLength + subLength;

    delX = (xRight - xLeft)/xNode;
    delY = (yRight - yLeft)/yNode;
    delZ = (zRight - zLeft)/zNode;
   
    for(int i=0; i< xNode + 2*wenoR; i++)
        xField[i] = xLeft + (i - wenoR + 0.5)*delX;
    for(int j=0; j< yNode + 2*wenoR; j++)
        yField[j] = yLeft + (j - wenoR + 0.5)*delY;
    for(int k=0; k< zNode + 2*wenoR; k++)
        zField[k] = zLeft + (k - wenoR + 0.5)*delZ;

    std::cout << "processor is " << rank << " and the total processor = " << this->size
              << ", we use the order of " << wenoR <<std::endl;
    MPI::COMM_WORLD.Barrier();

}

mesh::~mesh()
{
    if(xField != NULL)
        delete [] xField;
    xField = NULL;
    
    if(yField != NULL)
        delete [] yField;
    yField = NULL;

    if(zField != NULL)
        delete [] zField;
    zField = NULL;
}


#endif