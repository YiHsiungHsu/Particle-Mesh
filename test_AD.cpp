#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float ***buildGrid(int numRows, int numCols, int numLevels);
void mass_deposition( const int N, double *M, double *x, double *y, double *z, const double gs, const int GN, const int mode, float ****M_grid);
void acceleration_deposition( int N, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode, double *a);

int main(void){
    double gs = 1.0;
    int GN = 5;
    int N = 2;
    double M[N], x[N], y[N], z[N];
    M[0] = 10.0;
    x[0] = 1.3;
    y[0] = 2.2;
    z[0] = 0.8;
    M[1] = 5.0;
    x[1] = 2.3;
    y[1] = 3.2;
    z[1] = 1.8;
    int mode = 1;
    float ***M_grid;
    float ***a_grid = buildGrid(GN,GN,GN);
    mass_deposition(N, M, x, y, z, gs, GN, mode, &M_grid);
    // Note we need to allocate a place to pass the 3-d function.
    double a[N];
    
     for(int i = 0; i<GN; i++){
         for(int j = 0; j<GN; j++){
             for(int k = 0; k<GN; k++){
                 a_grid[i][j][k] = 1.0 + 1.0*(double)i+1.0*(double)j+1.0*(double)k;
                 // provides the values a where we should give the value of real a when implementation.
             }
         }
     }
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
     
    acceleration_deposition( N, a_grid, M_grid, M, x, y, z, gs, GN, mode, a);
    
    for(int n = 0; n<N; n++){
        printf( "a[%2d] = %5.3f \n",n, a[n] );
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

void acceleration_deposition( int N, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode, double *a)
{
/*
 mode 1: NGP
 mode 2: CIC
 mode 3: TSC
 N is total number of particles.
 gs is grid size.
 GN is total grid number. Note that I assume we have a equilateral box.
 a_grid is the input acceleration matrix.
 a is output acceleration of one component for every particle. (a zero matrix pointer)
*/
    for(int n = 0; n<N; n++){
        a[n] = 0.0;
    }
    float m[N][GN][GN][GN]; //allocated mass for every single particle with m[particle][gridx][gridy][gridz]
    double dx, dy, dz;
    double wx, wy, wz; //weighted function
    for(int n = 0; n<N; n++){
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    m[n][i][j][k] = 0.0;
                }
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
                            
                            m[n][i][j][k] = (float)wx*wy*wz*M[n];
                            if(m[n][i][j][k] != 0.0 && M_grid[i][j][k] != 0.0)
                           {
                                a[n] += (double)m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k];
                            }
                            
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
                            
                            m[n][i][j][k] = (float)wx*wy*wz*M[n];
                            if(m[n][i][j][k] != 0.0 && M_grid[i][j][k] != 0.0)
                            {
                                a[n] += (double)m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k];
                            }
                            
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
                            
                            m[n][i][j][k] = (float)wx*wy*wz*M[n];
                            if(m[n][i][j][k] != 0.0 && M_grid[i][j][k] != 0.0)
                            {
                                a[n] += (double)m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k];
                            }
                 
                        }
                    }
                }
            }
        }
}

void mass_deposition(int N, double *M, double *x, double *y, double *z, double gs, int GN, int mode, float ****M_grid)
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
                        (*M_grid)[i][j][k] += (float)m[n][i][j][k];
                        
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
                        (*M_grid)[i][j][k] += (float)m[n][i][j][k];
                        
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

                        (*M_grid)[i][j][k] += (float)m[n][i][j][k];
                    }
                }
            }
        }
    }
}
