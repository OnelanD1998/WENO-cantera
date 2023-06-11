
#ifndef CANTERATHERMOMIX_H
#define CANTERATHERMOMIX_H

#include "cantera/zerodim.h"
#include "cantera/transport.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "dictionary.H"
#include "mesh.H"
#include "variableField.H"
#include <vector>
#include <string>

class canteraThermoMix
{
public:
    dictionary canteraDic_;
    const std::string canteraMechanism_;
    mesh& mesh_;
    std::shared_ptr<Cantera::Solution> canteraSolution_;
    std::shared_ptr<Cantera::ThermoPhase> canteraGas_;
    std::vector<variableField> speciesField_;
    std::vector<std::string> speciesNames_;

    variableField rho_;
    variableField psi_;
    variableField Hs_;
    variableField Cv_;
    variableField Cp_;

    variableField& p_;
    variableField& T_;

public:
    canteraThermoMix(mesh& meshDo, variableField& p, variableField& T);
    ~canteraThermoMix();

    void initialThermo();
    void correctThermo();

    void speciesBoundary();

    int nSpecie(){ return canteraGas_->nSpecies();}

    variableField& rho(){return this->rho_;}
    variableField& psi(){return this->psi_;}
    variableField& Hs(){return this->Hs_;}
};




//------------------------- construct function ---------------------------------------//
canteraThermoMix::canteraThermoMix(mesh& meshDo, variableField& p, variableField& T)
:
canteraDic_("canteraDic"),
canteraMechanism_(canteraDic_.lookUp<std::string>("canteraMechanism")),
mesh_(meshDo),
canteraSolution_(Cantera::newSolution(canteraMechanism_, "")),
canteraGas_(canteraSolution_->thermo()),
rho_(mesh_),
psi_(mesh_),
Hs_("T", mesh_),
Cv_(mesh_),
Cp_(mesh_),
p_(p),
T_(T)
{
    for(int s=0; s<nSpecie(); s++)
    {
        speciesNames_.push_back(canteraGas_->speciesName(s));

        variableField tempSp(speciesNames_[s], mesh_);
        speciesField_.push_back(tempSp);

        std::cout << "specie " << s << " is " << speciesNames_[s] << std::endl;
    }

    std::cout << "cantera mechanism is " << canteraMechanism_ << std::endl;
}

//------------------------- destruct function ---------------------------------------//
canteraThermoMix::~canteraThermoMix()
{
}


//------------------------------------ member functions ------------------------//
void canteraThermoMix::initialThermo()
{
    for(int i=0; i< p_.xTotNode_; i++){
        for(int j=0; j< p_.yTotNode_; j++){
            for(int k=0; k< p_.zTotNode_; k++){
                //std::vector<double> yTemp;
                double* yTemp = new double[nSpecie()];
                for(int s=0; s<nSpecie(); s++)
                {
                    yTemp[s] = (speciesField_[s][i][j][k]);
                }
                canteraGas_->setState_TPY(T_[i][j][k], p_[i][j][k], yTemp);
                rho_[i][j][k] = canteraGas_->density();
                psi_[i][j][k] = canteraGas_->meanMolecularWeight()/canteraGas_->RT();
                Hs_[i][j][k] = canteraGas_->enthalpy_mass();
                Cv_[i][j][k] = canteraGas_->cv_mass();
                Cp_[i][j][k] = canteraGas_->cp_mass();

                delete [] yTemp;
            }
        }
    }
}

void canteraThermoMix::correctThermo()
{
    //
    if(mesh_.rank == 0)
        std::cout << "To correct the thermoproperties\n ";
        
    for(int i=0; i< p_.xTotNode_; i++){
        for(int j=0; j< p_.yTotNode_; j++){
            for(int k=0; k< p_.zTotNode_; k++){
                double* yTemp = new double[nSpecie()];

                for(int s=0; s<nSpecie(); s++)
                {
                    yTemp[s] = (speciesField_[s][i][j][k]);
                }
                canteraGas_->setState_PY(p_[i][j][k], yTemp);
                canteraGas_->setState_HP(Hs_[i][j][k], p_[i][j][k]); // setState_HP needs (J/kg)

                if(i < p_.mesh_.wenoR || i > p_.xTotNode_-p_.mesh_.wenoR-1
                || j < p_.mesh_.wenoR || j > p_.yTotNode_-p_.mesh_.wenoR-1
                || k < p_.mesh_.wenoR || k > p_.zTotNode_-p_.mesh_.wenoR-1) 
                    rho_[i][j][k] = canteraGas_->density();     //-- boundary

                T_[i][j][k] = canteraGas_->temperature();
                psi_[i][j][k] = canteraGas_->meanMolecularWeight()/canteraGas_->RT();
                //Hs_[i][j][k] = canteraGas_->enthalpy_mass();
                Cv_[i][j][k] = canteraGas_->cv_mass();
                Cp_[i][j][k] = canteraGas_->cp_mass();

                delete [] yTemp;
            }
        }
    }
}


void canteraThermoMix::speciesBoundary()
{
    for(int s=0; s<nSpecie(); s++)
    {
        speciesField_[s].boundaryCorrect();
    }
}
#endif