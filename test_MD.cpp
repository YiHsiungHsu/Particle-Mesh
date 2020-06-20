#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

float ***buildGrid(int numRows, int numCols, int numLevels);
void mass_deposition( int N, int Nthread,double *M, double *x, double *y, double *z, double gs, int GN, int mode, float ****M_grid);

int main(void){
    //mpi initialize
    int NRank, MyRank;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );
    MPI_Comm_size( MPI_COMM_WORLD, &NRank );
    
    double gs = 1.0;
    int GN = 5;
    int N = 2;
    double M[N], x[N], y[N], z[N];
    int Nthread = 1;
    /*
    srand( time(NULL) );//set random seed for creating random number
    for(int n = 0; n<N; n++){
        M[n] = 10.0*(double)rand()/RAND_MAX;//10 is maxium mass
        x[n] = (double)rand()/RAND_MAX*(GN-1);
        y[n] = (double)rand()/RAND_MAX*(GN-1);
        z[n] = (double)rand()/RAND_MAX*(GN-1);
    }*/
    M[0] = 10.0;
    M[1] = 5.0;
    x[0] = 1.5;
    x[1] = 2.7;
    y[0] = 3.8;
    y[1] = 3.2;
    z[0] = 0.4;
    z[1] = 2.1;
    int mode = 2;
    float ***M_grid;
     

    int sendcount = N/NRank;
    int RecvCount = SendCount;
    const int RootRank = 0;
    float *SendBuf = new float [N]; // only relevant for the root rank
    float *RecvBuf = new float [RecvCount];
    
    if(MyRank = RootRank){
        for(int n=0; n<N; n++){
            SendBuf [n] = M[n];
        }
    }
    
    MPI_Scatter( SendBuf, SendCount, MPI_INT,RecvBuf, RecvCount, MPI_FLOAT, RootRank, MPI_COMM_WORLD );
    
    mass_deposition(RecvCount, Nthread, RecvBuf, x, y, z, gs, GN, mode, &M_grid);
    delete [] SendBuf;
    delete [] RecvBuf;
    
    printf( "\nM_grid:\n" );
    for(int k = 0; k<GN; k++){
        printf( "k = %2d \n", k );
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                printf( "  %5.3f", M_grid[i][j][k] );
            }
            printf( "\n" );
        }
    }
    return 0;
}



float ***buildGrid(int numRows, int numCols, int numLevels)
{
    float ***levels;
    levels = (float * * *)malloc(numLevels *sizeof(float *)); //Contains all levels

    int rowIndex, levelIndex;

    for (levelIndex = 0; levelIndex < numLevels; levelIndex++)
    {
        float **level = (float * *)malloc(numRows * sizeof(float *)); //Contains all rows

        for(rowIndex = 0; rowIndex < numRows; rowIndex++)
        {
            level[rowIndex] = (float *)malloc(numCols * sizeof(float)); //Contains all columns
        }

        levels[levelIndex] = level;
    }

    return levels;
}

void mass_deposition(int N, int Nthread,double *M, double *x, double *y, double *z, double gs, int GN, int mode, float ****M_grid)
{
/*
 mode 1: NGP
 mode 2: CIC
 mode 3: TSC
 N is total number of particles.
 gs is grid size.
 GN is total grid number. Note that I assume we have a equilateral box.
 M_grid is the output matrix (a pointer) which gives allocated mass on every grid.
 M_grid has input of zero matrix.
*/
    double m[N][GN][GN][GN]; //allocated mass for every single particle with m[particle][gridx][gridy][gridz]
    double dx, dy, dz;
    double wx, wy, wz; //weighted function
    
    *M_grid = buildGrid(GN,GN,GN);
    // initialize M_grid
# pragma omp parallel num_threads( Nthread )
{
    # pragma omp parallel for collapse( 3 )
    for(int i = 0; i<GN; i++){
        for(int j = 0; j<GN; j++){
            for(int k = 0; k<GN; k++){
                (*M_grid)[i][j][k] = 0.0;
            }
        }
    }
    if(mode == 1)
    {
        for(int n = 0; n<N; n++)
        {
            for(int i = 0; i<GN; i++)
            {
                dx = fabs(x[n]-i*gs);
                
                if(dx<0.5*gs) wx = 1.0;
                else if(dx==0.5*gs) wx = 0.5;
                else wx = 0.0;
                
                for(int j = 0; j<GN; j++)
                {
                    dy = fabs(y[n]-j*gs);
                    
                    if(dy<0.5*gs) wy = 1.0;
                    else if(dy==0.5*gs) wy = 0.5;
                    else wy = 0.0;
                    
                    for(int k = 0; k<GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        
                        if(dz<0.5*gs) wz = 1.0;
                        else if(dz==0.5*gs) wz = 0.5;
                        else wz = 0.0;
                        
                        m[n][i][j][k] = wx*wy*wz*M[n];
                        
                    }
                }
            }
        }
    }
    if(mode == 2)
    {
        for(int n = 0; n<N; n++)
        {
            for(int i = 0; i<GN; i++)
            {
                dx = fabs(x[n]-i*gs);
                
                if(dx<gs) wx = (1.0-dx/gs);
                else wx = 0.0;
                
                for(int j = 0; j<GN; j++)
                {
                    dy = fabs(y[n]-j*gs);
                    
                    if(dy<gs) wy = (1.0-dy/gs);
                    else wy = 0.0;
                    
                    for(int k = 0; k<GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        
                        if(dz<gs) wz = (1.0-dz/gs);
                        else wz = 0.0;
                        
                        m[n][i][j][k] = wx*wy*wz*M[n];
                        
                    }
                }
            }
        }
    }
    if(mode == 3)
    {
        for(int n = 0; n<N; n++)
        {
            for(int i = 0; i<GN; i++)
            {
                dx = fabs(x[n]-i*gs);
                
                if(dx<0.5*gs)
                {
                    wx = 0.75-dx*dx/gs/gs;
                }
                else if(dx>=0.5*gs && dx <= 1.5*gs)
                {
                    wx = 0.5*(1.5-dx/gs)*(1.5-dx/gs);
                }
                else wx = 0.0;
                
                for(int j = 0; j<GN; j++)
                {
                    dy = fabs(y[n]-j*gs);
                    
                    if(dy<0.5*gs)
                    {
                        wy = 0.75-dy*dy/gs/gs;
                    }
                    else if(dy>=0.5*gs && dy <= 1.5*gs)
                    {
                        wy = 0.5*(1.5-dy/gs)*(1.5-dy/gs);
                    }
                    else wy = 0.0;
                    
                    for(int k = 0; k<GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        
                        if(dz<0.5*gs)
                        {
                            wz = 0.75-dz*dz/gs/gs;
                        }
                        else if(dz>=0.5*gs && dz <= 1.5*gs)
                        {
                            wz = 0.5*(1.5-dz/gs)*(1.5-dz/gs);
                        }
                        else wz = 0.0;
                        
                        m[n][i][j][k] = wx*wy*wz*M[n];

                    }
                }
            }
        }
    }
//# pragma omp barrier
# pragma omp parallel for collapse( 4 )
        for(int n = 0; n<N; n++){
            for(int i = 0; i<GN; i++){
                for(int j = 0; j<GN; j++){
                    for(int k = 0; k<GN; k++){
                        (*M_grid)[i][j][k] += (float)m[n][i][j][k];
                    }
                }
            }
        }
}
}
