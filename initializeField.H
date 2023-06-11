for(int i=meshDomain.wenoR; i< meshDomain.wenoR + meshDomain.xNode; i++){
        for(int j=meshDomain.wenoR; j< meshDomain.wenoR + meshDomain.yNode; j++){
            for(int k=meshDomain.wenoR; k< meshDomain.wenoR + meshDomain.zNode; k++)
            {
                U[i][j][k] = 0.0;
                V[i][j][k] = 0.0;
                W[i][j][k] = 0.0;
                p[i][j][k] = 101325.0;
                T[i][j][k] = 300.0;

                if(meshDomain.xField[i] > 0.5)
                {
                    p[i][j][k] = 501325.0;
                    T[i][j][k] = 1500.0;
                }
            }
        }
    }

