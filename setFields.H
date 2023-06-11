for(int k=meshDomain.wenoR; k< meshDomain.wenoR + meshDomain.zNode; k++)
            for(int j=meshDomain.wenoR; j< meshDomain.wenoR + meshDomain.yNode; j++)
                for(int i=meshDomain.wenoR; i< meshDomain.wenoR + meshDomain.xNode; i++)
                {
                    if(meshDomain.xField[i] < 0.5)
                    {
                        p[i][j][k] *= 10.;
                        T[i][j][k] = 348.432;
                        U[i][j][k] = 00.;
                    }
                }