        double minDt = 1000.0;
        double minGlobalDt = 10000.0;
        for(int k=meshDomain.wenoR; k< meshDomain.wenoR + meshDomain.zNode; k++)
                for(int j=meshDomain.wenoR; j< meshDomain.wenoR + meshDomain.yNode; j++)
                    for(int i=meshDomain.wenoR; i< meshDomain.wenoR + meshDomain.xNode; i++)
                    {
                        double tempLenth = std::min(std::min(meshDomain.delX, meshDomain.delY), meshDomain.delZ);
                        if(minDt > timeMain.CFL()*tempLenth/USum[i][j][k]) minDt = timeMain.CFL()*tempLenth/USum[i][j][k];
                    }

        MPI::COMM_WORLD.Reduce(&minDt, &minGlobalDt, 1, MPI::DOUBLE, MPI::MIN, 0);
        MPI::COMM_WORLD.Barrier();

        MPI::COMM_WORLD.Bcast(&minGlobalDt, 1, MPI::DOUBLE, 0);      
        MPI::COMM_WORLD.Barrier();  

        meshDomain.meshTime_.setTimeDelta(minGlobalDt);

        meshDomain.meshTime_.setTime(meshDomain.meshTime_.timeValue() + meshDomain.meshTime_.timeDelta());
        if(rank == 0)
        {
            std::cout << " minimum time delta is " << meshDomain.meshTime_.timeDelta() << std::endl;
            std::cout << "processor is " << rank << ", current time is " << meshDomain.meshTime_.timeValue() << std::endl;
            
        }
        