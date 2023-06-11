#include "mpi.h"

#ifndef MYFIELD_H
#define MYFIELD_H


template<class T>
class myPtrList
{
    T*** x;
public:
    unsigned i_;
    unsigned j_;
    unsigned k_;
public:
    myPtrList()
    :
    i_(0),
    j_(0),
    k_(0)
    {this->setSize(0,0,0);};
    myPtrList(unsigned columns, unsigned rows, unsigned floors)
    :
    i_(columns),
    j_(rows),
    k_(floors)
    {
        this->setSize(columns, rows, floors);
    }
    
    myPtrList(const myPtrList<T> &copyer)
    {
        i_ = copyer.i_;
        j_ = copyer.j_;
        k_ = copyer.k_;
        this->setSize(this->i_, this->j_, this->k_);
        for(int i = 0; i<this->i_; i++)
        {
            for(int j = 0; j<this->j_; j++)
            {
                for(int k = 0; k<this->k_; k++)
                {
                    x[i][j][k] = copyer.dataRead()[i][j][k];
                }
            }
        }

    }

    ~myPtrList()
    {
        if(x != nullptr)
        {
            for(int i=0 ; i<this->i_; i++)
            {
                for(int j=0 ; j<this->j_; j++)
                {
                    delete [] x[i][j];
                }
                delete [] x[i];
            }
            delete [] x;
        }
        x = nullptr;
    }

    T*** data(){return this->x;}

    T*** dataRead() const
    {return this->x;}

    void setSize(unsigned columns, unsigned rows, unsigned floors)
    {
        x = new T**[columns];
        for(unsigned i=0; i<columns; i++)
        {
            x[i] = new T*[rows];
            for(unsigned j=0; j<rows; j++)
                x[i][j] = new T[floors];
        }
    }

    inline T** operator[](const int n)
    {
        return x[n];
    }

    void operator ==(myPtrList<T> ptr1);

    friend myPtrList<double> operator+(myPtrList<double> ptr1,  myPtrList<double> ptr2);
    friend myPtrList<double> operator-(myPtrList<double> ptr1,  myPtrList<double> ptr2);
    friend myPtrList<double> operator*(myPtrList<double> ptr1,  myPtrList<double> ptr2);
    friend myPtrList<double> operator/(myPtrList<double> ptr1,  myPtrList<double> ptr2);
    friend myPtrList<double> operator*(double ptr1, myPtrList<double> ptr2);

    friend double maxValue(myPtrList<double> ptr);
    friend double minValue(myPtrList<double> ptr);

    friend myPtrList<double> sqrt(myPtrList<double> ptr);
    friend myPtrList<double> pos0(myPtrList<double> ptr);
    friend myPtrList<double> neg(myPtrList<double> ptr);
};

typedef myPtrList<double> doubleField; 


template<class T>
void myPtrList<T>::operator== (myPtrList<T> ptr1)
{
    for(int i=0; i<this->i_; i++)
    {
        for(int j=0; j<this->j_; j++)
        {
            for(int k=0; k<this->k_; k++)
            {
                this->x[i][j][k] = ptr1[i][j][k];
            }
        }
    }
}


//-------------------- friend functions ------------------------//
myPtrList<double> operator+(myPtrList<double> ptr1,  myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr1.i_, ptr1.j_, ptr1.k_);
    for(int i=0; i<ptr1.i_; i++)
    {
        for(int j=0; j<ptr1.j_; j++)
        {
            for(int k=0; k<ptr1.k_; k++)
            {
                tempPtr[i][j][k] = ptr1[i][j][k] + ptr2[i][j][k];
            }
        }
    }
    return tempPtr;
}

myPtrList<double> operator-(myPtrList<double> ptr1,  myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr1.i_, ptr1.j_, ptr1.k_);
    for(int i=0; i<ptr1.i_; i++)
    {
        for(int j=0; j<ptr1.j_; j++)
        {
            for(int k=0; k<ptr1.k_; k++)
            {
                tempPtr[i][j][k] = ptr1[i][j][k] - ptr2[i][j][k];
            }
        }
    }
    return tempPtr;
}

myPtrList<double> operator*(myPtrList<double> ptr1,  myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr1.i_, ptr1.j_, ptr1.k_);
    for(int i=0; i<ptr1.i_; i++)
    {
        for(int j=0; j<ptr1.j_; j++)
        {
            for(int k=0; k<ptr1.k_; k++)
            {
                tempPtr[i][j][k] = ptr1[i][j][k] * ptr2[i][j][k];
            }
        }
    }
    return tempPtr;
}

myPtrList<double> operator/(myPtrList<double> ptr1,  myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr1.i_, ptr1.j_, ptr1.k_);
    for(int i=0; i<ptr1.i_; i++)
    {
        for(int j=0; j<ptr1.j_; j++)
        {
            for(int k=0; k<ptr1.k_; k++)
            {
                tempPtr[i][j][k] = ptr1[i][j][k] / (ptr2[i][j][k] + 1e-9);
            }
        }
    }
    return tempPtr;
}

myPtrList<double> operator*(double ptr1, myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr2.i_, ptr2.j_, ptr2.k_);
    for(int i=0; i<ptr2.i_; i++)
    {
        for(int j=0; j<ptr2.j_; j++)
        {
            for(int k=0; k<ptr2.k_; k++)
            {
                tempPtr[i][j][k] = ptr1 + ptr2[i][j][k];
            }
        }
    }
    return tempPtr;
}


//------------------------------- friend functions ---------------------------//
double maxValue(myPtrList<double> ptr)
{

    double tempDou = -10000;
    double tempTMax = 0.;
    for(int i=0; i<ptr.i_; i++)
    {
        for(int j=0; j<ptr.j_; j++)
        {
            for(int k=0; k<ptr.k_; k++)
            {
                if(tempDou < ptr[i][j][k]) tempDou = ptr[i][j][k];
            }
        }
    }
    // if(MPI::Is_initialized())
    // {
    //     MPI::COMM_WORLD.Reduce(&tempDou, &tempTMax, 1, MPI::DOUBLE, MPI::MAX, 0);
    //    // MPI::COMM_WORLD.Barrier();

    //     MPI::COMM_WORLD.Bcast(&tempTMax, 1, MPI::DOUBLE, 0);      
    //    // MPI::COMM_WORLD.Barrier(); 
    // }

    return tempTMax;
}

double minValue(myPtrList<double> ptr)
{
    double tempDou = 10000000;
    double tempTMin = 0.;
    for(int i=0; i<ptr.i_; i++)
    {
        for(int j=0; j<ptr.j_; j++)
        {
            for(int k=0; k<ptr.k_; k++)
            {
                if(tempDou > ptr[i][j][k]) tempDou = ptr[i][j][k];
            }
        }
    }

   // MPI::COMM_WORLD.Reduce(&tempDou, &tempTMin, 1, MPI::DOUBLE, MPI::MIN, 0);

    //MPI::COMM_WORLD.Bcast(&tempTMin, 1, MPI::DOUBLE, 0);      
   // MPI::COMM_WORLD.Barrier(); 

    return tempTMin;
}

myPtrList<double> sqrt(myPtrList<double> ptr)
{
    myPtrList<double> tempPtr(ptr.i_, ptr.j_, ptr.k_);
    for(int i=0; i<ptr.i_; i++)
    {
        for(int j=0; j<ptr.j_; j++)
        {
            for(int k=0; k<ptr.k_; k++)
            {
                tempPtr[i][j][k] = std::sqrt(ptr[i][j][k]);
            }
        }
    }
    return tempPtr;
}

myPtrList<double> pos0(myPtrList<double> ptr)
{
    myPtrList<double> tempPtr(ptr.i_, ptr.j_, ptr.k_);
    for(int i=0; i<ptr.i_; i++)
    {
        for(int j=0; j<ptr.j_; j++)
        {
            for(int k=0; k<ptr.k_; k++)
            {
                if(ptr[i][j][k] >= 0)
                    tempPtr[i][j][k] = 1.;
                else
                    tempPtr[i][j][k] = 0.;
            }
        }
    }
    return tempPtr;
}

myPtrList<double> neg(myPtrList<double> ptr)
{
    myPtrList<double> tempPtr(ptr.i_, ptr.j_, ptr.k_);
    for(int i=0; i<ptr.i_; i++)
    {
        for(int j=0; j<ptr.j_; j++)
        {
            for(int k=0; k<ptr.k_; k++)
            {
                if(ptr[i][j][k] < 0)
                    tempPtr[i][j][k] = 1.;
                else
                    tempPtr[i][j][k] = 0.;
            }
        }
    }
    return tempPtr;
}



//------------------------------- static functions ---------------------------//
static myPtrList<double> min(myPtrList<double> ptr1, myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr2.i_, ptr2.j_, ptr2.k_);
    for(int i=0; i<ptr2.i_; i++)
    {
        for(int j=0; j<ptr2.j_; j++)
        {
            for(int k=0; k<ptr2.k_; k++)
            {
                if(ptr1[i][j][k] > ptr2[i][j][k]) tempPtr[i][j][k] = ptr2[i][j][k];
                else tempPtr[i][j][k] = ptr1[i][j][k];
            }
        }
    }
    return tempPtr;
}

static myPtrList<double> mag(myPtrList<double> ptr1)
{
    myPtrList<double> tempPtr(ptr1.i_, ptr1.j_, ptr1.k_);
    for(int i=0; i<ptr1.i_; i++)
    {
        for(int j=0; j<ptr1.j_; j++)
        {
            for(int k=0; k<ptr1.k_; k++)
            {
                tempPtr[i][j][k] = std::fabs(ptr1[i][j][k]);
            }
        }
    }
    return tempPtr;
}


static myPtrList<double> max(myPtrList<double> ptr1, myPtrList<double> ptr2)
{
    myPtrList<double> tempPtr(ptr2.i_, ptr2.j_, ptr2.k_);
    for(int i=0; i<ptr2.i_; i++)
    {
        for(int j=0; j<ptr2.j_; j++)
        {
            for(int k=0; k<ptr2.k_; k++)
            {
                if(ptr1[i][j][k] > ptr2[i][j][k]) tempPtr[i][j][k] = ptr1[i][j][k];
                else tempPtr[i][j][k] = ptr2[i][j][k];
            }
        }
    }
    return tempPtr;
}

#endif