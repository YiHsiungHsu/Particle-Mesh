#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void acceleration_deposition( int N, float ***a_grid, double *x, double *y, double *z, double gs, int GN, int mode, double *a)
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
    double dx, dy, dz;
    double wx, wy, wz; //weighted function
    
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
                        
                        a[n] += wx*wy*wz*a_grid[i][j][k];
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
                    for(int k = 0; k<GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        if(dx<=gs && dy<=gs && dz<=gs)
                        {
                            a[n] += (1.0-dx/gs)*(1.0-dy/gs)*(1.0-dz/gs)*a_grid[i][j][k];
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
                if(dx<=0.5*gs)
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
                    if(dy<=0.5*gs)
                    {
                        wy = 0.75-dy*dy/gs/gs;
                    }
                    else if(dy>=0.5*gs && dy <= 1.5*gs)
                    {
                        wy = 0.5*(1.5-dy/gs)*(1.5-dy/gs);
                    }
                    else wy = 0.0;
                    
                    for(int k = 0; k<=GN; k++)
                    {
                        dz = fabs(z[n]-k*gs);
                        if(dz<=0.5*gs)
                        {
                            wz = 0.75-dz*dz/gs/gs;
                        }
                        else if(dz>=0.5*gs && dz <= 1.5*gs)
                        {
                            wz = 0.5*(1.5-dz/gs)*(1.5-dz/gs);
                        }
                        else wz = 0.0;
                        
                        a[n] += wx*wy*wz*a_grid[i][j][k];
                        
                    }
                }
            }
        }
    }
}
