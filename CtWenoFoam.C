/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ComCellR2DyFoam

Description
    Density-based compressible reaction flow solver based on central-upwind schemes of
    Kurganov and Tadmor with reactions and record cell structure and wave system

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "mpi.h"

#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "myField.H"
#include "mesh.H"
#include "dictionary.H"
#include "variableField.H"
#include "canteraThermoMix.H"

#include <iomanip>
#include <vector>

#include <fstream>


using namespace Cantera;



int main(int argc, char* argv[])
{
    std::cout.precision(9);
    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();

    std::vector<double> rkA, rkB;
    rkA.push_back(1.);
    rkA.push_back(0.25);
    rkA.push_back(2./3.);

    rkB.push_back(0.);
    rkB.push_back(0.75);
    rkB.push_back(1./3.);
    
    //-- initialize
    //#include "initializeField.H"
    mesh meshDomain;
    timeTerm& timeMain(meshDomain.meshTime_);

    variableField U("U", meshDomain);
    variableField V("V", meshDomain);
    variableField W("W", meshDomain);
    variableField p("p", meshDomain);
    variableField T("T", meshDomain);

    canteraThermoMix ctMix(meshDomain, p, T);

    #include "setFields.H"

    //-- boundary set
    if(rank == 0)
        std::cout << "The first correct boundary start\n";
    p.boundaryCorrect();
    U.boundaryCorrect();
    V.boundaryCorrect();
    W.boundaryCorrect();
    T.boundaryCorrect();
   

    ctMix.speciesBoundary();
    ctMix.initialThermo();
    ctMix.correctThermo();

    variableField& rho = ctMix.rho();
    variableField& psi = ctMix.psi();
    variableField& Hs = ctMix.Hs();
    variableField e(Hs - p/rho);
    variableField rhoE(
                  rho*e 
                + 0.5*(U*U+V*V+W*W)*rho
                );
    variableField rhoU(rho*U);
    variableField rhoV(rho*V);
    variableField rhoW(rho*W);

    std::cout << "save old field\n";

    variableField rhoOld(meshDomain);
    variableField rhoUOld(meshDomain);
    variableField rhoVOld(meshDomain);
    variableField rhoWOld(meshDomain);
    variableField rhoEOld(meshDomain);

    doubleField phix = 0.*rho.wenoRecX_neg();
    doubleField phiUx = 0.*rho.wenoRecX_neg();
    doubleField phiEx = 0.*rho.wenoRecX_neg();


    for(int timeStep=0; timeStep < timeMain.overStep(); timeStep++)
    {
        if(rank == 0)
            std::cout << "\n\n||++++++++++++++++++Current step is " << timeStep << " +++++++++++++++++++++++||\n";

        rhoOld == rho;
        rhoUOld == rhoU;
        rhoVOld == rhoV;
        rhoWOld == rhoW;
        rhoEOld == rhoE;

        variableField a_sonic(sqrt(ctMix.Cp_/ctMix.Cv_/psi));
        variableField USum(sqrt(U*U + V*V + W*W) + a_sonic);

        #include "setTimeDel.H"
        MPI::COMM_WORLD.Barrier();
        
        for(int rk=0; rk< rkA.size(); rk++)
        {   
            a_sonic == (sqrt(ctMix.Cp_/ctMix.Cv_/psi));

            std::string fluxScheme(timeMain.timeDic().lookUp<std::string>("fluxScheme"));
            if(fluxScheme == "KNPScheme")
            {
                #include "KNPScheme.H"
            }
            else if(fluxScheme == "LXFScheme")
            {
                #include "LXFScheme.H"
            }
            else if(fluxScheme == "HLLCScheme")
            {
                #include "HLLCScheme.H"
            }
            else return 0;

            MPI::COMM_WORLD.Barrier();

            if(rank == 0)
                std::cout << "++++++++++++++++++solve the consavation value +++++++++++++++++++++++\n";

            rho == rkB[rk]*rhoOld
                + rkA[rk]*(rho - timeMain.timeDelta()*pow(meshDomain.delX ,-1)*projectX(phix, meshDomain));

            rhoU == rkB[rk]*rhoUOld
                  + rkA[rk]*(rhoU - timeMain.timeDelta()
                    *pow(meshDomain.delX ,-1)*projectX(phiUx, meshDomain));

            rhoV == rhoV;

            rhoE == rkB[rk]*rhoEOld
                 + rkA[rk]*(
                rhoE - timeMain.timeDelta()
                *pow(meshDomain.delX ,-1)*projectX(phiEx, meshDomain)
                );

            U == (rhoU/rho);
            U.boundaryCorrect();

            e == (rhoE/rho - 0.5*(U*U + V*V + W*W));
            p == (rho/psi);
            p.boundaryCorrect();

            Hs == (e + p/rho);
            Hs.boundaryCorrect();

            ctMix.correctThermo();
            e == Hs - p/rho;

            MPI::COMM_WORLD.Barrier();
        }
        
        //-- save fields
        if(timeStep % timeMain.saveStep() == 0)
             #include "saveFields.H"

        if(rank == 0)
        {
            std::cout << "Y-H2 =====" << ctMix.speciesField_[1][1][3][3] << std::endl;
            std::cout << "T =====" << rho[1][3][3] << std::endl;
            std::cout << "maximum phix =====" << maxValue(phix) << std::endl;
            std::cout << "maximum phiUx =====" << minValue(phiUx) << std::endl;
            std::cout << "minimum T =====" << minValue(T.field_) << std::endl;
            std::cout << "calculate complete" << std::endl;
        }
        MPI::COMM_WORLD.Barrier();
        

    }
    
    MPI::Finalize();
    return 0;
}

// ************************************************************************* //
