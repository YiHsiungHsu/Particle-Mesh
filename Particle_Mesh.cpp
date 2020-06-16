#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//include any package you need and omp or mpi

float ***buildGrid(const int numRows, const int numCols, const int numLevels); //creat grid points
void mass_deposition( const int N, double *M, double *x, double *y, double *z, const double gs, const int GN, const int mode, float ****M_grid);
void acceleration_deposition( int N, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode, double *a);
void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double gs, const int GN  );

int main(void){
    // constants
    // add any necessary const. here
    const double gs = 1.0; // grid size (distance between every grid point)
    const int GN = 5; //box size. I set equilateral box, tell me if you need to change it.
    const int N = 1; // number of particles
    double M[N], x[N], y[N], z[N];
    const int mode_d = 3; // choose the mode for deposition
    const double t_end = 100.0; // end time
    const double dt = 0.05; //time step size of each step
    // end constants
    
    // initial conditions
    // add any necessary IC for your function
    // I set random mass and random position for now. You can change them if you want but remember removing this note.
    double ti = 0.0; //initial time
    double t = ti;
    srand( time(NULL) );// set random seed for creating random number
    for(int n = 0; n<N; n++){
        M[n] = 10.0*(double)rand()/RAND_MAX;// 10 is maxium mass
        x[n] = (double)rand()/RAND_MAX*(GN-1);
        y[n] = (double)rand()/RAND_MAX*(GN-1);
        z[n] = (double)rand()/RAND_MAX*(GN-1);
    }
    // end IC
    
    while(t <= t_end)
    {
        // mass deposition
        // Note that the output of this is 3 by 3 matrix which from M_grid[0][0][0] to M[GN-1][GN-1][GN-1]
        float ***M_grid;
        mass_deposition(N, M, x, y, z, gs, GN, mode_d, &M_grid);
        // end mass deposition
        
        // calculate potential here
        // Read the output of mass deposition carefully, and please inform me if I need to change the form of mass_grid.
        // end potential
        
        // Gradient of potential
        // tell me the form of the output of potential calculation
        // End Gradient potential
    
        // acceleration deposition here
        // I expect my output to be ax[N], ay[N], az[N]
        double ax[N], ay[N], az[N];
        float ***a_grid = buildGrid(GN,GN,GN);
        //assign a_grid for x here
        acceleration_deposition( int N, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode_d, double *ax);
        //assign a_grid for y here
        acceleration_deposition( int N, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode_d, double *ay);
        //assign a_grid for z here
        acceleration_deposition( int N, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode_d, double *az);
        // end acceleration deopsotion
        
        // calculate jerk
        double axcopy[N], aycopy[N], azcopy[N]; // copy acceleration for next step;
        double aax[N], aay[N], aaz[N]; // jerk for x, y, z
        if(t == ti) // modified with initial time
        {
            for(n = 0; n < N ; n++)
            {
                axcopy[n] = ax[n];
                aycopy[n] = ay[n];
                azcopy[n] = az[n];
            }
        }
        for(n = 0; n < N ; n++)
        {
            aax[n] = (ax[n] - axcopy[n])/dt;
            aay[n] = (ay[n] - aycopy[n])/dt;
            aaz[n] = (az[n] - azcopy[n])/dt;
            axcopy[n] = ax[n];
            aycopy[n] = ay[n];
            azcopy[n] = az[n];
        }
        //end jerk
        
        // Hermite Integral, DKD, KDK
        // Read the output of acceleration deposition and see if there should be any change.
        // end HI, DKD, KDK
        
        // Dump data
        
        // end dump data
        
        t += dt;
    }
    return 0;
}



float ***buildGrid(const int numRows, const int numCols, const int numLevels)
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
    double m[N][N][N][N]; //allocated mass for every single particle with m[particle][gridx][gridy][gridz]
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
    double m[N][N][N][N]; //allocated mass for every single particle with m[particle][gridx][gridy][gridz]
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
                            
                            m[n][i][j][k] = wx*wy*wz*M[n];
                            if(m[n][i][j][k] != 0.0 && M_grid[i][j][k] != 0.0)
                            {
                                a[n] += m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k];
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
                            
                            m[n][i][j][k] = wx*wy*wz*M[n];
                            if(m[n][i][j][k] != 0.0 && M_grid[i][j][k] != 0.0)
                            {
                                a[n] += m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k];
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
                            
                            m[n][i][j][k] = wx*wy*wz*M[n];
                            if(m[n][i][j][k] != 0.0 && M_grid[i][j][k] != 0.0)
                            {
                                a[n] += m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k];
                            }
                 
                        }
                    }
                }
            }
        }
    
    void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double gs, const int GN  )
{
float dt[N];         //timestep based on all particles' properties respectively.
double rx;
double ry;
double rz;
double rvx;
double rvy;
double rvz;
double afx[N];//f stand for final, the predicted acceleration at time t+ts
double afy[N];
double afz[N];
double jfx[N];//j stand for jerk
double jfy[N];
double jfz[N];
double a2x[N]; //second order derivative of acceleration
double a2y[N];
double a2z[N];
double a3x[N]; //third order derivative of acceleration
double a3y[N];
double a3z[N];
double af2x[N];   //second derivative of acceleration at t + dt
double af2y[N];
double af2z[N];
double aaf[N];   //first a of all below stands for absolute value
double aaf2[N];
double ajf[N];
double aa3[N];
double rex = 0;
double rey = 0;
double rez = 0;

  for (int i=0; i<N; i++){ //initialize value to zero
      aaf[i] = 0;
      aaf2[i] = 0;
      ajf[i] = 0;
      aa3[i] = 0;
      afx[i] = 0; 
      afy[i] = 0;
      afz[i] = 0;
      jfx[i] = 0;
      jfy[i] = 0;
      jfz[i] = 0;
      a2x[i] = 0;
      a2y[i] = 0;
      a2z[i] = 0;
      a3x[i] = 0;
      a3y[i] = 0;
      a3z[i] = 0;
      af2x[i] = 0;
      af2y[i] = 0;
      af2z[i] = 0;
                          }

   for (int n=0; n<N; n++){//first position, velocity update (drift and kick)
      x[n] += ( ts*vx[n] + pow(ts, 2)*ax[n]/2 + pow(ts, 3)*jx[n]/6 );
      y[n] += ( ts*vy[n] + pow(ts, 2)*ay[n]/2 + pow(ts, 3)*jy[n]/6 );
      z[n] += ( ts*vz[n] + pow(ts, 2)*az[n]/2 + pow(ts, 3)*jz[n]/6 );
      vx[n] += ( ts*ax[n] + pow(ts, 2)*jx[n]/2 );
      vy[n] += ( ts*ay[n] + pow(ts, 2)*jy[n]/2 );
      vz[n] += ( ts*az[n] + pow(ts, 2)*jz[n]/2 );
			  }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = fabs(x[n]-x[j]);
      ry = fabs(y[n]-y[j]);
      rz = fabs(z[n]-z[j]);
      rex = pow(rx, 2)+pow(sp, 2);
      rey = pow(ry, 2)+pow(sp, 2);
      rez = pow(rz, 2)+pow(sp, 2);
      afx[n] += G*M[j]*rx/sqrt(pow(rex, 3));
      afy[n] += G*M[j]*ry/sqrt(pow(rey, 3));
      afz[n] += G*M[j]*rz/sqrt(pow(rez, 3));
      rvx = fabs(vx[n] - vx[j]);
      rvy = fabs(vy[n] - vy[j]);
      rvz = fabs(vz[n] - vz[j]);
      dafx[n] += G*M[j]*(rvx/sqrt(pow(rex, 3)) + 3*(rvx*rx)*rx/sqrt(pow(rex, 5)));
      dafy[n] += G*M[j]*(rvy/sqrt(pow(rey, 3)) + 3*(rvy*ry)*ry/sqrt(pow(rey, 5)));
      dafz[n] += G*M[j]*(rvz/sqrt(pow(rez, 3)) + 3*(rvz*rz)*rz/sqrt(pow(rez, 5)));
			   }

//high order correction
   for (int n=0; n<N; n++){
      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*jx[n] + 2*jfx[n]) )/(pow(ts, 2)+1e-7);
      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*jy[n] + 2*jfy[n]) )/(pow(ts, 2)+1e-7);
      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*jz[n] + 2*jfz[n]) )/(pow(ts, 2)+1e-7);
      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(jx[n] + jfx[n]) )/(pow(ts, 3)+1e-7);
      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(jy[n] + jfy[n]) )/(pow(ts, 3)+1e-7);
      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(jz[n] + jfz[n]) )/(pow(ts, 3)+1e-7);
                          }

//final acceleration, position and velocity
   for (int n=0; n<N; n++){
      ax[n] += ( ts*dax[n] + pow(ts, 2)*a2x[n]/2 + pow(ts, 3)*a3x[n]/6 );
      ay[n] += ( ts*day[n] + pow(ts, 2)*a2y[n]/2 + pow(ts, 3)*a3y[n]/6 );
      az[n] += ( ts*daz[n] + pow(ts, 2)*a2z[n]/2 + pow(ts, 3)*a3z[n]/6 );
      jx[n] = jfx[n];
      jy[n] = jfy[n];
      jz[n] = jfz[n];
      x[n] += ( pow(ts, 4)*a2x[n]/24 + pow(ts, 5)*a3x[n]/120 );
      y[n] += ( pow(ts, 4)*a2y[n]/24 + pow(ts, 5)*a3y[n]/120 );
      z[n] += ( pow(ts, 4)*a2z[n]/24 + pow(ts, 5)*a3z[n]/120 );
      vx[n] += ( pow(ts, 3)*a2x[n]/6 + pow(ts, 4)*a3x[n]/24 );
      vy[n] += ( pow(ts, 3)*a2y[n]/6 + pow(ts, 4)*a3y[n]/24 );
      vz[n] += ( pow(ts, 3)*a2z[n]/6 + pow(ts, 4)*a3z[n]/24 );
      af2x[n] = a2x[n] + ts*a3x[n];
      af2y[n] = a2y[n] + ts*a3y[n];
      af2z[n] = a2z[n] + ts*a3z[n];
                          }
       
//new timestep to next step
   for (int n=0; n<N; n++){
      aaf[n]  = sqrt((pow(ax[n], 2) + pow(ay[n], 2) + pow(az[n], 2)));
      ajf[n] = sqrt((pow(jx[n], 2) + pow(jy[n], 2) + pow(jz[n], 2)));
      aa3[n]  = sqrt((pow(a3x[n], 2) + pow(a3y[n], 2) + pow(a3z[n], 2)));
      aaf2[n] = sqrt((pow(af2x[n], 2) + pow(af2y[n], 2) + pow(af2z[n], 2)));
      dt[n]   = sqrt(0.01*((aaf[n]*aaf2[n] + pow(ajf[n], 2))/(ajf[n]*aa3[n] + pow(aaf2[n], 2) + 1e-7)));
			  }

   ts = dt[0];

   for (int i=0; i<N; i++){
   if  ( dt[i] < ts )   {
       ts = dt[i];
       if ( ts < 1e-6 ) ts = 1e-7;
			              }
                          }
}
