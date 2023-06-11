#include "mesh.H"
#include "dictionary.H"
#include <sstream>
#include <string>

#ifndef VARIABELFIELD_H
#define VARIABELFIELD_H

class variableField
{
public:
    std::string fieldName_;

    mesh& mesh_;
    int xTotNode_;
    int yTotNode_;
    int zTotNode_;
    doubleField field_;

    std::vector<std::string> boundaryType;
    std::vector<double> boundaryValue;
    //dictionary fieldFile_;

public:
    variableField(mesh& meshDom);

    variableField(std::string fieldName, mesh& meshDom);

    ~variableField();


    //-- reconstruct face flux in x direction
    doubleField wenoRecX_neg();
    doubleField wenoRecX_pos();
    //variableField projectX(doubleField faceField, mesh& meshDom);

    void boundaryCorrect();

    //---------------- operator functions ----------------//
    variableField operator+(variableField var1);
    variableField operator-(variableField var1);
    variableField operator*(variableField var1);
    variableField operator/(variableField var1);
    variableField operator==(variableField var1);

    friend variableField operator*(double dou1, variableField var1);
    friend variableField sqrt(variableField var1);

    inline double** operator[](const int n) 
    {
        return field_[n];
    }
};

//-------------------------- constructor function ---------------------------------------//
variableField::variableField(mesh& meshDom)
:
mesh_(meshDom),
xTotNode_(mesh_.xNode + 2*mesh_.wenoR),
yTotNode_(mesh_.yNode + 2*mesh_.wenoR),
zTotNode_(mesh_.zNode + 2*mesh_.wenoR),
field_(xTotNode_, yTotNode_, zTotNode_),
fieldName_("")
{
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                this->field_[i][j][k] = 0.0;
}

variableField::variableField(std::string fieldName, mesh& meshDom)
:
variableField(meshDom)
{
    
    fieldName = meshDom.meshTime_.timeName() + "/" + fieldName;
    dictionary fieldFile(fieldName);
    //fieldFile_(fieldFile);
    fieldName_ = fieldName;
    
    std::cout << "\n|+++++++++++++++ construct from file ++++++++++++++++++|\n";
    double initialValue = fieldFile.lookUp<double>("initialValue");
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                this->field_[i][j][k] = initialValue;
    std::cout << fieldName << ", INIT V == " << initialValue << std::endl;


    std::cout << "|+++++++++++++++ boundary type ++++++++++++++++++|\n";
    for(int bNum=0; bNum < 6; bNum++)
    {
        std::stringstream fieldNameStream;
        std::string nBoundary;
        fieldNameStream << bNum;
        fieldNameStream >> nBoundary;
        //nBoundary = "boundary" + nBoundary; 
        boundaryType.push_back(fieldFile.lookUp<std::string>("boundary" + nBoundary));

        if(boundaryType[bNum] == "fixedValue")
        {
            boundaryValue.push_back(fieldFile.lookUp<double>("fixedValue" + nBoundary));
            //std::cout << fieldName << ", boundary type == " << boundaryType[bNum] << " and value is " << boundaryValue[bNum] << std::endl;
        }
        else
        {
                boundaryValue.push_back(0.);
                //std::cout << fieldName << ", boundary type == " << boundaryType[bNum] << std::endl;
        }
        
    }
}

variableField::~variableField()
{
}


//----------------------- member functions-----------------------------//
doubleField variableField::wenoRecX_neg()
{
    doubleField faceField(mesh_.xNode + 1, mesh_.yNode, mesh_.zNode);
    for(int i=0; i<mesh_.xNode + 1; i++)
    {
        for(int j=0;j< mesh_.yNode; j++)
        {
            for(int k=0;k< mesh_.zNode; k++)
            {
                std::vector<double> stencil(mesh_.wenoR*2 -1);
                for(int st=0; st<stencil.size(); st++) stencil[st] = this->field_[i+st][j+mesh_.wenoR][k+mesh_.wenoR];

                //-- smooth omdocators
                double b0n = 13./12.*pow(stencil[0] - 2.*stencil[1] + stencil[2], 2.)
                           + 1./4.*pow(stencil[0] - 4.*stencil[1] + 3.*stencil[2], 2.);

                double b1n = 13./12.*pow(stencil[1] - 2.*stencil[2] + stencil[3], 2.)
                           + 1./4.*pow(stencil[1] - stencil[3], 2.);

                double b2n = 13./12.*pow(stencil[2] - 2.*stencil[3] + stencil[4], 2.)
                           + 1./4.*pow(3.*stencil[2] - 4.*stencil[3] + stencil[4], 2.);
                
                double af0n = pow(b0n + 1e-9, -2)/10.;
                double af1n = pow(b1n + 1e-9, -2)/10.*6.;
                double af2n = pow(b2n + 1e-9, -2)/10.*3.;
                double afSum = af0n + af1n + af2n;

                //-- stencils weigths
                double w0n = af0n/afSum,
                       w1n = af1n/afSum,
                       w2n = af2n/afSum;
                
                faceField[i][j][k] = w0n*(2.*stencil[0] - 7.*stencil[1] + 11.*stencil[2])/6.
                                   + w1n*(-stencil[1] + 5.*stencil[2] + 2.*stencil[3])/6.
                                   + w2n*(2.*stencil[2] + 5.*stencil[3] - stencil[4])/6.;

            }
        }
    }

    return faceField;
}

doubleField variableField::wenoRecX_pos()
{
    doubleField faceField(mesh_.xNode + 1, mesh_.yNode, mesh_.zNode);
    for(int i=0; i<mesh_.xNode + 1; i++)
    {
        for(int j=0;j< mesh_.yNode; j++)
        {
            for(int k=0;k< mesh_.zNode; k++)
            {
                std::vector<double> stencil(mesh_.wenoR*2 -1);
                for(int st=0; st<stencil.size(); st++) stencil[st] = this->field_[i+st+1][j+mesh_.wenoR][k+mesh_.wenoR];

                //-- smooth omdocators
                double b0n = 13./12.*pow(stencil[0] - 2.*stencil[1] + stencil[2], 2.)
                           + 1./4.*pow(stencil[0] - 4.*stencil[1] + 3.*stencil[2], 2.);

                double b1n = 13./12.*pow(stencil[1] - 2.*stencil[2] + stencil[3], 2.)
                           + 1./4.*pow(stencil[1] - stencil[3], 2.);

                double b2n = 13./12.*pow(stencil[2] - 2.*stencil[3] + stencil[4], 2.)
                           + 1./4.*pow(3.*stencil[2] - 4.*stencil[3] + stencil[4], 2.);
                
                double af0n = pow(b0n + 1e-9, -2)/10.*3.;
                double af1n = pow(b1n + 1e-9, -2)/10.*6.;
                double af2n = pow(b2n + 1e-9, -2)/10.;
                double afSum = af0n + af1n + af2n;

                //-- stencils weigths
                double w0n = af0n/afSum,
                       w1n = af1n/afSum,
                       w2n = af2n/afSum;
                
                faceField[i][j][k] = w0n*(-stencil[0] + 5.*stencil[1] + 2.*stencil[2])/6.
                                   + w1n*(2.*stencil[1] + 5.*stencil[2] - stencil[3])/6.
                                   + w2n*(11.*stencil[2] - 7.*stencil[3] + 2.*stencil[4])/6.;

            }
        }
    }

    return faceField;
}



void variableField::boundaryCorrect()
{
    //-- boundary in y direction
    for(int k=mesh_.wenoR; k<mesh_.wenoR + mesh_.zNode; k++){
        for(int j=mesh_.wenoR; j<mesh_.wenoR + mesh_.xNode; j++){
            for(int i=0; i<mesh_.wenoR; i++){
                int typeIndex = 0;
                if(boundaryType[typeIndex] == "slipWall")
                    this->field_[j][mesh_.wenoR-i-1][k] = this->field_[j][mesh_.wenoR+i][k];
                else if (boundaryType[typeIndex] == "reflectWall")
                    this->field_[j][mesh_.wenoR-i-1][k] = -this->field_[j][mesh_.wenoR+i][k];
                else if (boundaryType[typeIndex] == "zeroGradient")
                    this->field_[j][mesh_.wenoR-i-1][k] = this->field_[j][mesh_.wenoR][k];
                else if(boundaryType[typeIndex] == "fixedValue")
                    this->field_[j][mesh_.wenoR-i-1][k] = boundaryValue[typeIndex];
                else;
                
                typeIndex = 1;
                if(boundaryType[typeIndex] == "slipWall")
                    this->field_[j][mesh_.wenoR+mesh_.yNode+i][k] = this->field_[j][mesh_.wenoR+mesh_.yNode-1-i][k];
                else if (boundaryType[typeIndex] == "reflectWall")
                    this->field_[j][mesh_.wenoR+mesh_.yNode+i][k] = -this->field_[j][mesh_.wenoR+mesh_.yNode-1-i][k];
                else if (boundaryType[typeIndex] == "zeroGradient")
                    this->field_[j][mesh_.wenoR+mesh_.yNode+i][k] = this->field_[j][mesh_.wenoR+mesh_.yNode-1][k];
                else if(boundaryType[typeIndex] == "fixedValue")
                    this->field_[j][mesh_.wenoR+mesh_.yNode+i][k] = boundaryValue[typeIndex];
                else;
            }
        }
    }   


    //-- boundary in x direction
    double* tempSR1 = new double[mesh_.wenoR], 
          * tempSR2 = new double[mesh_.wenoR];
    for(int pcrIndex = -1; pcrIndex< mesh_.size; pcrIndex++)
    {
        for(int k=mesh_.wenoR; k<mesh_.wenoR + mesh_.zNode; k++){
            for(int j=mesh_.wenoR; j<mesh_.wenoR + mesh_.yNode; j++){
               
                if(mesh_.rank == pcrIndex)
                {
                    for(int i=0; i<mesh_.wenoR; i++){
                        
                        int typeIndex = 3;
                        if(pcrIndex != mesh_.size - 1)  //-- no trans if latest core 
                        {
                            tempSR1[i] = this->field_[mesh_.wenoR+mesh_.xNode-1-i][j][k];    
                        }
                        else if(boundaryType[typeIndex] == "slipWall")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = this->field_[mesh_.wenoR+mesh_.xNode-1-i][j][k];

                        else if (boundaryType[typeIndex] == "reflectWall")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = -this->field_[mesh_.wenoR+mesh_.xNode-1-i][j][k];

                        else if (boundaryType[typeIndex] == "zeroGradient")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = this->field_[mesh_.wenoR+mesh_.xNode-1][j][k];

                        else if(boundaryType[typeIndex] == "fixedValue")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = boundaryValue[typeIndex];
                        else;
                    }
                    if(mesh_.size > 1 && mesh_.rank != mesh_.size - 1) MPI::COMM_WORLD.Send(tempSR1, mesh_.wenoR, MPI::DOUBLE, pcrIndex + 1, 80);
                    if(mesh_.size > 1 && mesh_.rank != mesh_.size - 1) MPI::COMM_WORLD.Recv(tempSR2, mesh_.wenoR, MPI::DOUBLE, pcrIndex + 1, 90);

                    for(int i=0; i<mesh_.wenoR; i++){
                        
                        int typeIndex = 3;
                        if(pcrIndex != mesh_.size - 1)
                        {
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = tempSR2[i];
                        }
                        else if(boundaryType[typeIndex] == "slipWall")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = this->field_[mesh_.wenoR+mesh_.xNode-1-i][j][k];
                        else if (boundaryType[typeIndex] == "reflectWall")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = -this->field_[mesh_.wenoR+mesh_.xNode-1-i][j][k];
                        else if (boundaryType[typeIndex] == "zeroGradient")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = this->field_[mesh_.wenoR+mesh_.xNode-1][j][k];
                        else if(boundaryType[typeIndex] == "fixedValue")
                            this->field_[mesh_.wenoR+mesh_.xNode+i][j][k] = boundaryValue[typeIndex];
                        else;
                    }

                }
                else if(mesh_.rank == pcrIndex + 1)
                {
                    if(mesh_.size > 1 && mesh_.rank != 0) MPI::COMM_WORLD.Recv(tempSR1, mesh_.wenoR, MPI::DOUBLE, pcrIndex, 80);
                    
                    for(int i=0; i<mesh_.wenoR; i++){
                        
                        int typeIndex = 2;
                        if((pcrIndex+1) != 0) 
                        {
                            this->field_[mesh_.wenoR-i-1][j][k] = tempSR1[i];

                            tempSR2[i] = this->field_[mesh_.wenoR+i][j][k];
                        }
                        else if(boundaryType[typeIndex] == "slipWall")
                            this->field_[mesh_.wenoR-i-1][j][k] = this->field_[mesh_.wenoR+i][j][k];
                        else if (boundaryType[typeIndex] == "reflectWall")
                            this->field_[mesh_.wenoR-i-1][j][k] = -this->field_[mesh_.wenoR+i][j][k];
                        else if (boundaryType[typeIndex] == "zeroGradient")
                            this->field_[mesh_.wenoR-i-1][j][k] = this->field_[mesh_.wenoR][j][k];
                        else if(boundaryType[typeIndex] == "fixedValue")
                            this->field_[mesh_.wenoR-i-1][j][k] = boundaryValue[typeIndex];
                        else;
                    }
                    if(mesh_.size > 1 && mesh_.rank != 0) MPI::COMM_WORLD.Send(tempSR2, mesh_.wenoR, MPI::DOUBLE, pcrIndex, 90);
                }
            
            }
        }   
    }
    delete [] tempSR1, tempSR2;


    //-- boundary in z direction
    for(int k=mesh_.wenoR; k<mesh_.wenoR + mesh_.yNode; k++){
        for(int j=mesh_.wenoR; j<mesh_.wenoR + mesh_.xNode; j++){
            for(int i=0; i<mesh_.wenoR; i++){
                int typeIndex = 4;
                if(boundaryType[typeIndex] == "slipWall")
                    this->field_[j][k][mesh_.wenoR-i-1] = this->field_[j][k][mesh_.wenoR+i];
                else if (boundaryType[typeIndex] == "reflectWall")
                    this->field_[j][k][mesh_.wenoR-i-1] = -this->field_[j][k][mesh_.wenoR+i];
                else if (boundaryType[typeIndex] == "zeroGradient")
                    this->field_[j][k][mesh_.wenoR-i-1] = this->field_[j][k][mesh_.wenoR];
                else if(boundaryType[typeIndex] == "fixedValue")
                    this->field_[j][k][mesh_.wenoR-i-1] = boundaryValue[typeIndex];
                else;
                
                typeIndex = 5;
                if(boundaryType[typeIndex] == "slipWall")
                    this->field_[j][k][mesh_.wenoR+mesh_.zNode+i] = this->field_[j][k][mesh_.wenoR+mesh_.zNode-1-i];
                else if (boundaryType[typeIndex] == "reflectWall")
                    this->field_[j][k][mesh_.wenoR+mesh_.zNode+i] = -this->field_[j][k][mesh_.wenoR+mesh_.zNode-1-i];
                else if (boundaryType[typeIndex] == "zeroGradient")
                    this->field_[j][k][mesh_.wenoR+mesh_.zNode+i] = this->field_[j][k][mesh_.wenoR+mesh_.zNode-1];
                else if(boundaryType[typeIndex] == "fixedValue")
                    this->field_[j][k][mesh_.wenoR+mesh_.zNode+i] = boundaryValue[typeIndex];
                else;
            }
        }
    }  

    MPI::COMM_WORLD.Barrier(); 

}



//---------------------------------------- operator function------------------------------------------------//
variableField variableField::operator+(variableField var1)
{
    variableField tempVar(mesh_);
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                tempVar.field_[i][j][k] = this->field_[i][j][k]
                                        + var1.field_[i][j][k];
    return tempVar;
}

variableField variableField::operator-(variableField var1)
{
    variableField tempVar(mesh_);
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                tempVar.field_[i][j][k] = this->field_[i][j][k]
                                        - var1.field_[i][j][k];
    return tempVar;
}

variableField variableField::operator*(variableField var1)
{
    variableField tempVar(mesh_);
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                tempVar.field_[i][j][k] = this->field_[i][j][k]
                                        * var1.field_[i][j][k];
    return tempVar;
}

variableField variableField::operator/(variableField var1)
{
    variableField tempVar(mesh_);
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                tempVar.field_[i][j][k] = this->field_[i][j][k]
                                        / (var1.field_[i][j][k] + 1e-16);
    return tempVar;
}

variableField variableField::operator==(variableField var1)
{
    for(int i=0; i< xTotNode_; i++)
        for(int j=0; j< yTotNode_; j++)
            for(int k=0; k< zTotNode_; k++)
                this->field_[i][j][k] = var1.field_[i][j][k];
    return *this;
}

//---------------- friend function------------------//

variableField operator*(double dou1, variableField var1)
{
    variableField tempVar(var1.mesh_);
    for(int i=0; i< tempVar.xTotNode_; i++)
        for(int j=0; j< tempVar.yTotNode_; j++)
            for(int k=0; k< tempVar.zTotNode_; k++)
                tempVar.field_[i][j][k] = dou1*var1.field_[i][j][k];
    return tempVar;
}

variableField sqrt(variableField var1)
{
    variableField tempVar(var1.mesh_);
    for(int i=0; i< tempVar.xTotNode_; i++)
        for(int j=0; j< tempVar.yTotNode_; j++)
            for(int k=0; k< tempVar.zTotNode_; k++)
                tempVar.field_[i][j][k] = std::sqrt(var1.field_[i][j][k]);
    return tempVar;
}

static variableField projectX(doubleField& faceField, mesh& meshDom)
{
    variableField tempPro(meshDom);
    tempPro == (0.*tempPro);
    //std::cout << "check point 1\n";
    for(int i=meshDom.wenoR; i<meshDom.wenoR+meshDom.xNode; i++)
    {
        for(int j=meshDom.wenoR; j<meshDom.wenoR+meshDom.yNode; j++)
        {
            for(int k=meshDom.wenoR; k<meshDom.wenoR+meshDom.zNode; k++)
            {
                tempPro[i][j][k] = faceField[i-meshDom.wenoR+1][j-meshDom.wenoR][k-meshDom.wenoR]
                                 - faceField[i-meshDom.wenoR][j-meshDom.wenoR][k-meshDom.wenoR];
            }
        }
    }
    return tempPro;
}

#endif