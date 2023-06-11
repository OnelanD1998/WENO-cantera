
    {	
    //   double temp1[101][101][10]; //�Ѿ��ĳɶ�̬������	
    //-- name
        using namespace std;
        stringstream saveFile;
        string fileName;
        saveFile << "output/" << rank << "_" << timeStep << "_a.dat" ;
        saveFile >> fileName;
        //sprintf(str1, "output3D/output%06d_%d.dat",step,myid);
        //-- output
    
        ofstream dataFile(fileName.c_str());
        // dataFile << "TITLE = \"Result (time = 0" << ")\" \n";
        // dataFile << "VARIABLES = \"X\", \"Y\", \"Z\", \"D\", \"U\", \"V\", \"W\", \"P\", \"S1\", \"S2\", \"T\" \n";
        // dataFile << "ZONE i = " << meshDomain.xNode
        //         << ", j = " << meshDomain.yNode
        //         << ", k = " << meshDomain.zNode << ", F = POINT\n";
    
        for(int k=meshDomain.wenoR; k< meshDomain.wenoR + meshDomain.zNode; k++)
            for(int j=meshDomain.wenoR; j< meshDomain.wenoR + meshDomain.yNode; j++)
                for(int i=meshDomain.wenoR; i< meshDomain.wenoR + meshDomain.xNode; i++)
                {

                    //T[i][j][k]=P[i][j][k]/D[i][j][k]/(CV*(GAMMA - 1.0));
                    dataFile << meshDomain.xField[i] << setw(24)
                            <<  meshDomain.yField[j] << setw(24)
                            <<  meshDomain.zField[k] << setw(24)
                            << U[i][j][k] << setw(24)
                            << V[i][j][k] << setw(24)
                            << W[i][j][k] << setw(24)
                            << p[i][j][k] << setw(24)
                            << T[i][j][k] << "\n";
                }
           

        dataFile.close();
    }//end void save
    
    MPI::COMM_WORLD.Barrier();