#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

//include any package you need and omp or mpi

float ***buildGrid(const int numRows, const int numCols, const int numLevels); //creat grid points
void mass_deposition( const int N, int Nthread,double *M, double *x, double *y, double *z, const double gs, const int GN, const int mode, float ****M_grid);
void acceleration_deposition( int N, int Nthread, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode, double *a);
void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G );
void hermiteDKD( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G );
void hermiteKDK( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G );
void Potential( double *rho, double *phi, double G, int BC, int GN, double gs , int Nthread );
void hermiteMPI( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G );
void herMPIDKD( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G );
void herMPIKDK( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G );

int main( int argc, char *argv[] ){

double start, end;

start = MPI_Wtime();

   int NRank, MyRank;

   MPI_Init( &argc, &argv );

   MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );

    // constants
    // add any necessary const. here
    int Nthread = 8;
    const double gs = 1.0; // grid size (distance between every grid point)
    const int GN = 32; //box size. I set equilateral box, tell me if you need to change it.
    const int N = 2; // number of particles
    double M[N], x[N], y[N], z[N];
    const int mode_d = 1; // choose the mode for deposition
    const int mode_h = 4; // choose the mode for hermite
    const double t_end = 10.0; // end time
    const double ts = 0.01; //time step size of each step
    const double G = 0.25/M_PI; //(m3 kg-1 s-2)
    const int BC = 0;// choose boundary condition (0=isolated 1=period)
    int c = 0;
    // end constants
	    
    
    // initial conditions
    // add any necessary IC for your function
    // I set random mass and random position for now. You can change them if you want but remember removing this note.
    double ti = 0.0; //initial time
    double t = ti;
    double vx[N], vy[N], vz[N], jx[N], jy[N], jz[N]; // jerk for x, y, z. j stand for jerk	
    srand( time(NULL) );// set random seed for creating random number
    for(int n = 0; n<N; n++){
        M[n] = 100.0;//*(double)rand()/RAND_MAX;// 10 is maxium mass
        x[n] = 10 + (double) n;//(double)rand()/RAND_MAX*(GN-1);
        y[n] = 10 + (double) n;//(double)rand()/RAND_MAX*(GN-1);
        z[n] = 10 + (double) n;//(double)rand()/RAND_MAX*(GN-1);
        vx[n] = (double) n;// (double)rand()/RAND_MAX*(GN-1);
        vy[n] = (double) n;// (double)rand()/RAND_MAX*(GN-1);
        vz[n] = (double) n;// (double)rand()/RAND_MAX*(GN-1);
        jx[n] = 0;
        jy[n] = 0;
        jz[n] = 0;
    }
    // end IC
    
    //set space
    float ***M_grid;
    int rhoN = GN*GN*GN;
    double rho[rhoN], phi_grid[GN+2][GN+2][GN+2];
    int index;
    int phiN = GN*GN*GN;
    double phi[phiN];
    float phi_dx[GN][GN][GN], phi_dy[GN][GN][GN], phi_dz[GN][GN][GN];
    double ax[N], ay[N], az[N];
    float ***a_grid = buildGrid(GN,GN,GN);
    double Mx = 0.0;
    double My = 0.0;
    double Mz = 0.0;
    FILE *file;
    //end set space
    
    while(t <= t_end)
    {
        // mass deposition
        // Note that the output of this is 3 by 3 matrix which from M_grid[0][0][0] to M[GN-1][GN-1][GN-1]
        mass_deposition(N, Nthread, M, x, y, z, gs, GN, mode_d, &M_grid);
        // end mass deposition
        
        // calculate potential here
	//get rho in row-major form
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    index = k+GN*(j+GN*i);
                    rho[index] = M_grid[i][j][k];
                }
            }
        }
	for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    index = k+GN*(j+GN*i);
                    phi[index] = 0;
                }
            }
        }
        Potential(  rho,  phi , G, BC, GN, gs,Nthread);
        // I need the row-major rho matrix and a row-major phi matrix with all 0
        // it will get the row-major phi matrix return
        for(int i = 1; i<GN+1; i++){
            for(int j = 1; j<GN+1; j++){
                for(int k = 1; k<GN+1; k++){
                    index = (k-1)+GN*((j-1)+GN*(i-1));
                    phi_grid[i][j][k] = phi[index];
                }
            }
        }
        //padding boundary for phi
        for(int l = 0; l<GN+2; l++){
            for(int m = 0; m<GN+2; m++){
                phi_grid[0][l][m] = 2.0*phi_grid[1][l][m] - phi_grid[2][l][m];
                phi_grid[GN+1][l][m] = 2.0*phi_grid[GN][l][m] - phi_grid[GN-1][l][m];
                phi_grid[l][0][m] = 2.0*phi_grid[l][1][m] - phi_grid[l][2][m];
                phi_grid[l][GN+1][m] = 2.0*phi_grid[l][GN][m] - phi_grid[l][GN-1][m];
                phi_grid[l][m][0] = 2.0*phi_grid[l][m][1] - phi_grid[l][m][2];
                phi_grid[l][m][GN+1] = 2.0*phi_grid[l][m][GN] - phi_grid[l][m][GN-1];
            }
        }
        // end potential
        
        // Gradient of potential
        //gradient inside
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    phi_dx[i][j][k] = (phi_grid[i+2][j+1][k+1]-phi_grid[i][j+1][k+1])/(2*gs);
                    phi_dy[i][j][k] = (phi_grid[i+1][j+2][k+1]-phi_grid[i+1][j][k+1])/(2*gs);
                    phi_dz[i][j][k] = (phi_grid[i+1][j+1][k+2]-phi_grid[i+1][j+1][k])/(2*gs);
                }
            }
        }

        // End Gradient potential
    
        // acceleration deposition here
        // I expect my output to be ax[N], ay[N], az[N]
        //assign a_grid for x here
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    a_grid[i][j][k] = phi_dx[i][j][k];
                }
            }
        }

        acceleration_deposition( N, Nthread, a_grid, M_grid, M, x, y, z, gs, GN, mode_d, ax);
        //assign a_grid for y here
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    a_grid[i][j][k] = phi_dy[i][j][k];
                }
            }
        }
        acceleration_deposition( N, Nthread, a_grid, M_grid, M, x, y, z, gs, GN, mode_d, ay);
        //assign a_grid for z here
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    a_grid[i][j][k] = phi_dz[i][j][k];
                }
            }
        }
        acceleration_deposition( N, Nthread, a_grid, M_grid, M, x, y, z, gs, GN, mode_d, az);
        
        // end acceleration deopsotion
        // Hermite Integral, DKD, KDK
        // Read the output of acceleration deposition and see if there should be any change.
        //HI
        if(mode_h == 1) hermite( N, M, x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz, ts, G );
        //DKD Hermite
        if(mode_h == 2) hermiteDKD( N, M, x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz, ts, G );
        //KDK Hermite
        if(mode_h == 3) hermiteKDK( N, M, x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz, ts, G );
	//MPI HI
	if(mode_h == 4) hermiteMPI( N, M, x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz, ts, G );
	//MPI DKD HI
        if(mode_h == 5) herMPIDKD( N, M, x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz, ts, G );
	//MPI KDK HI
	if(mode_h == 6) herMPIKDK( N, M, x, y, z, vx, vy, vz, ax, ay, az, jx, jy, jz, ts, G );
        // end HI, DKD, KDK
        
        //Momentum
        for(int n = 0; n<N; n++)
        {
            Mx += M[n] * vx[n];
            My += M[n] * vy[n];
            Mz += M[n] * vz[n];
        }
        //end Momentum
        
        // Dump data
        if (MyRank == 0){
        if (c%50 == 0){
	file = fopen("Particle_position.txt","ab");
        for(int n = 0; n < N; n++)
        {
            fprintf(file, "%5.5f \t %5.5f \t %5.5f \n", x[n], y[n], z[n]);
        }
        fclose(file);
        file = fopen("Momentum.txt","ab");
        fprintf(file, "%5.5f \t %5.5f \t %5.5f \n", Mx, My, Mz);
        fclose(file);
	}}
            Mx = 0;
            My = 0;
            Mz = 0;

        // end dump data
        c += 1;
        t += ts;
    }
    end = MPI_Wtime();
    if (MyRank == 0){
    file = fopen("elapse_time.txt","ab");
    fprintf(file, "%lf\n" , end -start );
    fclose(file);
    }
    MPI_Finalize();
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
    # pragma omp parallel num_threads( Nthread )
       {
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

void acceleration_deposition( int N, int Nthread, float ***a_grid, float ***M_grid, double *M, double *x, double *y, double *z, double gs, int GN, int mode, double *a)
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
    float m[N][GN][GN][GN]; //allocated mass for every single particle with m[particle][gridx][gridy][gridz]
    double dx, dy, dz;
    double wx, wy, wz; //weighted function
    # pragma omp parallel num_threads( Nthread )
    {
         # pragma omp for
    for(int n = 0; n<N; n++){
        a[n] = 0.0;
    }
        # pragma omp parallel for collapse( 4 )
    for(int n = 0; n<N; n++){
        for(int i = 0; i<GN; i++){
            for(int j = 0; j<GN; j++){
                for(int k = 0; k<GN; k++){
                    m[n][i][j][k] = 0.0;
                }
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
                                a[n] += (double)m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k]/M[n];
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
                                a[n] += (double)m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k]/M[n];
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
                                a[n] += (double)m[n][i][j][k]/M_grid[i][j][k]*a_grid[i][j][k]/M[n];
                            }
                 
                        }
                    }
                }
            }
        }
}
//---------------------------------------------------------------------------------------------
////Function    :  potential  
////Description :  use the DFT for self-gravity to solve the poisson solver to get the potential
////
////Note        :  still need to declare N box size
////                                     G gravitaional constant
////               I use the row-major matrix which the fftw need.
////Parameter   :  rho = density of particle 
////               phi = potential of particle gravity
////            
////
////Return      : phi
////--------------------------------------------------------------------------------------------
void Potential( double *rho, double *phi, double G, int BC, int GN, double gs, int Nthread )
{

  if(BC == 0)                                       //period BC
  {
    //fft rho to rhok
    fftw_complex *rhok;
    rhok = (fftw_complex*) fftw_malloc( GN*GN*(GN/2+1) * sizeof(fftw_complex) );
    fftw_plan fft;
    fft = fftw_plan_dft_r2c_3d( GN, GN, GN, rho, rhok, FFTW_ESTIMATE);
    fftw_execute(fft);
    //calculat potential phi =-4*M_PI*G/(kx**2+ky**2+kz**2)
    //need normailze with 1/N*N*N
    double _n;
    _n = 1 / (GN*gs*GN*gs*GN*gs);                                //normalize factor
    fftw_complex *phik;
    phik = (fftw_complex*) fftw_malloc( GN*GN*(GN/2+1) * sizeof(fftw_complex) );

    for(int i = 0; i < GN; i++){
    for(int j = 0; j < GN; j++){
    for(int k = 0; k < (GN/2+1); k++)
    {
        double _k;
        double kx,ky,kz;
        if (i > GN/2)
        {   
           kx = GN-i;
        }
        else
        {
           kx = i;
        }
        if (j > GN/2)
        {
           ky = GN-j;
        }
        else
        {
           ky = j;
        }
        kz = k;
	//set up the phi in row major-form
	
        _k = -1/gs/gs/((kx*kx*2*M_PI/GN/gs*2*M_PI/GN/gs)+(ky*ky*2*M_PI/GN/gs*2*M_PI/GN/gs)+(kz*kz*2*M_PI/GN/gs*2*M_PI/GN/gs));
     //   printf("k = %2f %2f %2f ,_k = %5.5f\n ",kx,ky,kz,_k);
        phik[k+(GN/2+1)*(j+GN*i)][0] = 4*M_PI*G*_k*_n*rhok[k+(GN/2+1)*(j+GN*i)][0];  //real part
        phik[k+(GN/2+1)*(j+GN*i)][1] = 4*M_PI*G*_k*_n*rhok[k+(GN/2+1)*(j+GN*i)][1];  //imagine part
    }}}
   
    //DC = 0
    phik[0][0] = 0; 
    phik[0][1] = 0;
    //fft phik to phi
    fftw_plan ifft;
    ifft = fftw_plan_dft_c2r_3d( GN, GN, GN, phik, phi, FFTW_ESTIMATE);  
    fftw_execute(ifft);

    for(int i = 0; i < GN; i++){
    for(int j = 0; j < GN; j++){
    for(int k = 0; k < (GN/2+1); k++)
    {
      //printf("r =%5.5f\n", phik[k+(GN/2+1)*(j+GN*i)][0]);
      //printf("i =%5.5f\n", phik[k+(GN/2+1)*(j+GN*i)][1]);
    }}}

    for(int i = 0; i < GN; i++){
    for(int j = 0; j < GN; j++){
    for(int k = 0; k < GN; k++)
    { 
      //printf("phi =%5.5f\n", phi[k+(GN)*(j+GN*i)]);
    }}}
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    fftw_free(rhok);
    fftw_free(phik);
  }  
  if (BC == 1)           //isolated boundary
  {
     //zero padding M
     double *zM;         //zero padding M     
     zM = (double*) fftw_malloc( (2*GN)*(2*GN)*(2*GN) * sizeof(double) );
 
     for (int i = 0; i < 2*GN; i++)
     for (int j = 0; j < 2*GN; j++)
     for (int k = 0; k < 2*GN; k++)
     {
         zM[k+(2*GN)*(j+2*GN*i)] = 0.0;
     }
     for (int i = 0; i < GN; i++)
     for (int j = 0; j < GN; j++)
     for (int k = 0; k < GN; k++)
     {
         zM[k+(2*GN)*(j+2*GN*i)] = rho[k+GN*(j+GN*i)]*gs*gs*gs;
     }
     double *dgf;        //discrete Green's function
     dgf = (double*) fftw_malloc( (2*GN)*(2*GN)*(2*GN) * sizeof(double) );   // dgf = -1*/R

     for (int i = 0; i < 2*GN; i++)
     for (int j = 0; j < 2*GN; j++)
     for (int k = 0; k < 2*GN; k++)
     {
         double nx,ny,nz;
         if (i > GN )
            {
               nx = 2*GN - i;
            }
         else
            {
               nx = i;
            }
         if (i > GN )
            {
               ny = 2*GN - j;
            }
         else
            {
               ny = j;
            }

         if (i > GN )
            {
               nz = 2*GN - k;
            }
         else
            {
               nz = k;
            }
         double _R;
         _R = -1 / gs*sqrt(nx*nx + ny*ny + nz*nz);
         if (i == 0 && j == 0 && k == 0)
            {
               dgf[k+(2*GN)*(j+2*GN*i)] = 0.0;
            }
         else
            {
                dgf[k+(2*GN)*(j+2*GN*i)] = _R;
            }
      }   
      //  FFT
      fftw_complex *zMk, *dgfk;
      zMk = (fftw_complex*) fftw_malloc( (2*GN)*(2*GN)*(GN+1) * sizeof(fftw_complex) );
      dgfk = (fftw_complex*) fftw_malloc( (2*GN)*(2*GN)*(GN+1) * sizeof(fftw_complex) );
      fftw_plan fftM, fftR;
      fftM = fftw_plan_dft_r2c_3d( 2*GN, 2*GN, 2*GN, zM, zMk, FFTW_ESTIMATE );
      fftR = fftw_plan_dft_r2c_3d( 2*GN, 2*GN, 2*GN, dgf, dgfk, FFTW_ESTIMATE );
      fftw_execute(fftM);
      fftw_execute(fftR);
      // convolution to get phi  ( a+bi * c+di )
      fftw_complex *conk,*phik;
      conk = (fftw_complex*) fftw_malloc( (2*GN)*(2*GN)*(GN+1) * sizeof(fftw_complex) );
      phik = (fftw_complex*) fftw_malloc( (2*GN)*(2*GN)*(GN+1) * sizeof(fftw_complex) );

      for (int i = 0; i < 2*GN; i++)
      for (int j = 0; j < 2*GN; j++)
      for (int k = 0; k < GN+1; k++)
      {
          conk[k+(GN+1)*(j+2*GN*i)][0] = (zMk[k+(GN+1)*(j+2*GN*i)][0] * dgfk[k+(GN+1)*(j+2*GN*i)][0]) - (zMk[k+(GN+1)*(j+2*GN*i)][1] * dgfk[k+(GN+1)*(j+2*GN*i)][1]);  // real part 
          conk[k+(GN+1)*(j+2*GN*i)][1] = (zMk[k+(GN+1)*(j+2*GN*i)][0] * dgfk[k+(GN+1)*(j+2*GN*i)][1]) + (zMk[k+(GN+1)*(j+2*GN*i)][1] * dgfk[k+(GN+1)*(j+2*GN*i)][0]);  // imagine part
      }
      double _n;
      _n = 1/(2*GN*2*GN*2*GN );      //normailize factor
 
      for (int i = 0; i < 2*GN; i++)
      for (int j = 0; j < 2*GN; j++)
      for (int k = 0; k < GN+1; k++)
      {
          phik[k+(GN+1)*(j+2*GN*i)][0] = G*_n*conk[k+(GN+1)*(j+2*GN*i)][0];  //real part
          phik[k+(GN+1)*(j+2*GN*i)][1] = G*_n*conk[k+(GN+1)*(j+2*GN*i)][1];  //imagine part
      }
      double *_phi; 
      _phi = (double*) fftw_malloc( (2*GN)*(2*GN)*(2*GN) * sizeof(double) );
      fftw_plan ifft;
      ifft = fftw_plan_dft_c2r_3d( 2*GN, 2*GN, 2*GN, phik, _phi, FFTW_ESTIMATE );
      fftw_execute(ifft);
 
      for (int i = 0; i < GN; i++)
      for (int j = 0; j < GN; j++)
      for (int k = 0; k < GN; k++)
      {
         phi[k+GN*(j+GN*i)] = _phi[k+(2*GN)*(j+2*GN*i)];
      }
      fftw_destroy_plan(fftM);
      fftw_destroy_plan(fftR);
      fftw_destroy_plan(ifft);
      fftw_free(zM);
      fftw_free(zMk);
      fftw_free(dgf);
      fftw_free(dgfk);
      fftw_free(conk);
      fftw_free(phik);
      fftw_free(_phi);
  }
}


void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G )
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
double sp = 0.5;
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
//printf("x = %lf\n", x[0]);

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      ry = y[j]-y[n];
      rz = z[j]-z[n];
//      rex = pow(rx, 2)+pow(sp, 2);
//      rey = pow(ry, 2)+pow(sp, 2);
//      rez = pow(rz, 2)+pow(sp, 2);
      if (fabs(rx)>sp) afx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else afx[n] += 0;
      if (fabs(ry)>sp) afy[n] += G*M[j]*ry/fabs(pow(ry, 3));
	else afy[n] += 0;
      if (fabs(rz)>sp) afz[n] += G*M[j]*rz/fabs(pow(rz, 3));
	else afz[n] += 0;
      rvx = vx[j] - vx[n];
      rvy = vy[j] - vy[n];
      rvz = vz[j] - vz[n];
      if (fabs(rx)>sp) jfx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
	else jfx[n] += 0;
      if (fabs(ry)>sp) jfy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
	else jfy[n] += 0;
      if (fabs(rz)>sp) jfz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
	else jfz[n] += 0;
			   }
   else                    {
      afx[n] += 0;
      jfx[n] += 0;
      afy[n] += 0;
      jfy[n] += 0;
      afz[n] += 0;
      jfz[n] += 0;
                           }

   			   }
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
/*
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
*/
}
void hermiteDKD( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G )
{
float dt[N];         //timestep based on all particles' properties respectively.
double rx;
double ry;
double rz;
double rvx;
double rvy;
double rvz;
double ahx[N];
double ahy[N];
double ahz[N];
double jhx[N];
double jhy[N];
double jhz[N];
double afx[N];//f stand for final, the predicted acceleration at time t+ts
double afy[N];
double afz[N];
double jfx[N];
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
double hts = ts/2;
double sp = 0.1;

//parameter initialization
   for (int i=0; i<N; i++){
      aaf[i] = 0;
      aaf2[i] = 0;
      ajf[i] = 0;
      aa3[i] = 0;
      ahx[i] = 0;
      ahy[i] = 0;
      ahz[i] = 0;
      jhx[i] = 0;
      jhy[i] = 0;
      jhz[i] = 0;
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

   for (int n=0; n<N; n++){//first drift of 0.5ts     
      x[n] += ( hts*vx[n] + pow(hts, 2)*ax[n]/2 + pow(hts, 3)*jx[n]/6 );
      y[n] += ( hts*vy[n] + pow(hts, 2)*ay[n]/2 + pow(hts, 3)*jy[n]/6 );
      z[n] += ( hts*vz[n] + pow(hts, 2)*az[n]/2 + pow(hts, 3)*jz[n]/6 );
		      }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      ry = y[j]-y[n];
      rz = z[j]-z[n];
      if ( fabs(rx)>sp ) ahx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else ahx[n] += 0;
      if ( fabs(ry)>sp ) ahy[n] += G*M[j]*ry/fabs(pow(ry, 3));
        else ahy[n] += 0;
      if ( fabs(rz)>sp ) ahz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else ahz[n] += 0;
      rvx = vx[j] - vx[n];
      rvy = vy[j] - vy[n];
      rvz = vz[j] - vz[n];
      if ( fabs(rx)>sp ) jhx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
	else jhx[n] += 0;
      if ( fabs(ry)>sp ) jhy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else jhy[n] += 0;
      if ( fabs(rz)>sp ) jhz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else jhz[n] += 0;
		           }
   else                    {
      ahx[n] += 0;
      jhx[n] += 0;
      ahy[n] += 0;
      jhy[n] += 0;
      ahz[n] += 0;
      jhz[n] += 0;
                           }
		      }
		      }

//kick and second drift
   for (int n=0; n<N; n++){
      vx[n] += ( ts*ahx[n] + pow(ts, 2)*jhx[n]/2 );
      vy[n] += ( ts*ahy[n] + pow(ts, 2)*jhy[n]/2 );
      vz[n] += ( ts*ahz[n] + pow(ts, 2)*jhz[n]/2 );
      x[n] += ( hts*vx[n] + pow(hts, 2)*ahx[n]/2 + pow(hts, 3)*jhx[n]/6 );
      y[n] += ( hts*vy[n] + pow(hts, 2)*ahy[n]/2 + pow(hts, 3)*jhy[n]/6 );
      z[n] += ( hts*vz[n] + pow(hts, 2)*ahz[n]/2 + pow(hts, 3)*jhz[n]/6 );
		      }

//predict a and jerk at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      ry = y[j]-y[n];
      rz = z[j]-z[n];
      if ( fabs(rx)>sp ) afx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else afx[n] += 0;
      if ( fabs(ry)>sp ) afy[n] += G*M[j]*ry/fabs(pow(ry, 3));
        else afy[n] += 0;
      if ( fabs(rz)>sp ) afz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else afx[n] += 0;
      rvx = vx[j] - vx[n];
      rvy = vy[j] - vy[n];
      rvz = vz[j] - vz[n];
      if ( fabs(rx)>sp ) jfx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
        else jfx[n] += 0;
      if ( fabs(ry)>sp ) jfy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else jfy[n] += 0;
      if ( fabs(rz)>sp ) jfz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else jfz[n] += 0;
			   }
   else                    {
      afx[n] += 0;
      jfx[n] += 0;
      afy[n] += 0;
      jfy[n] += 0;
      afz[n] += 0;
      jfz[n] += 0;
                           }
   			   }
   			   }

//high order correction
   for(int n=0; n<N; n++){
      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*jx[n] + 2*jfx[n]) )/(pow(ts, 2)+1e-7);
      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*jy[n] + 2*jfy[n]) )/(pow(ts, 2)+1e-7);
      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*jz[n] + 2*jfz[n]) )/(pow(ts, 2)+1e-7);
      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(jx[n] + jfx[n]) )/(pow(ts, 3)+1e-7);
      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(jy[n] + jfy[n]) )/(pow(ts, 3)+1e-7);
      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(jz[n] + jfz[n]) )/(pow(ts, 3)+1e-7);
                         }

//final acceleration, position and velocity
  for (int n=0; n<N; n++){
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
/*
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
*/
}
void hermiteKDK( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G )
{
float dt[N];         //timestep based on all particles' properties respectively.
double rx;
double ry;
double rz;
double rvx;
double rvy;
double rvz;
double ahx[N];
double ahy[N];
double ahz[N];
double jhx[N];
double jhy[N];
double jhz[N];
double afx[N];//f stand for final, the predicted acceleration at time t+ts
double afy[N];
double afz[N];
double jfx[N];
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
double hts = ts/2;
double sp = 0.1;
   for (int i=0; i<N; i++){//parameters initialization
      aaf[i] = 0;
      aaf2[i] = 0;
      ajf[i] = 0;
      aa3[i] = 0;
      ahx[i] = 0;
      ahy[i] = 0;
      ahz[i] = 0;
      jhx[i] = 0;
      jhy[i] = 0;
      jhz[i] = 0;
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

//first kick of 0.5ts
   for (int n=0; n<N; n++){
      vx[n] += ( hts*ax[n] + pow(hts, 2)*jx[n]/2 );
      vy[n] += ( hts*ay[n] + pow(hts, 2)*jy[n]/2 );
      vz[n] += ( hts*az[n] + pow(hts, 2)*jz[n]/2 );
		      }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)       {
      rx = x[j]-x[n];
      ry = y[j]-y[n];
      rz = z[j]-z[n];
      if ( fabs(rx)>sp ) ahx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else ahx[n] += 0;
      if ( fabs(ry)>sp ) ahy[n] += G*M[j]*ry/fabs(pow(ry, 3));
        else ahy[n] += 0;
      if ( fabs(rz)>sp ) ahz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else ahz[n] += 0;
      rvx = vx[j] - vx[n];
      rvy = vy[j] - vy[n];
      rvz = vz[j] - vz[n];
      if ( fabs(rx)>sp ) jhx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
        else jhx[n] += 0;
      if ( fabs(ry)>sp ) jhy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else jhy[n] += 0;
      if ( fabs(rz)>sp ) jhz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else ahz[n] += 0;
			   }
   else                    {
      ahx[n] += 0;
      jhx[n] += 0;
      ahy[n] += 0;
      jhy[n] += 0;
      ahz[n] += 0;
      jhz[n] += 0;
                           }
   			   }
			   }

//drift and second kick
   for (int n=0; n<N; n++){
      x[n] += ( ts*vx[n] + pow(ts, 2)*ahx[n]/2 + pow(ts, 3)*jhx[n]/6 );
      y[n] += ( ts*vy[n] + pow(ts, 2)*ahy[n]/2 + pow(ts, 3)*jhy[n]/6 );
      z[n] += ( ts*vz[n] + pow(ts, 2)*ahz[n]/2 + pow(ts, 3)*jhz[n]/6 );
      vx[n] += ( hts*ahx[n] + pow(hts, 2)*jhx[n]/2 );
      vy[n] += ( hts*ahy[n] + pow(hts, 2)*jhy[n]/2 );
      vz[n] += ( hts*ahz[n] + pow(hts, 2)*jhz[n]/2 );
		      }

//predict a and jerk at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      ry = y[j]-y[n];
      rz = z[j]-z[n];
      if (fabs(rx)>sp) afx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else afx[n] += 0;
      if (fabs(ry)>sp) afy[n] += G*M[j]*ry/fabs(pow(ry, 3));
        else afy[n] += 0;
      if (fabs(rz)>sp) afz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else afz[n] += 0;
      rvx = vx[j] - vx[n];
      rvy = vy[j] - vy[n];
      rvz = vz[j] - vz[n];
      if (fabs(rx)>sp) jfx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
        else jfx[n] += 0;
      if (fabs(ry)>sp) jfy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else jfy[n] += 0;
      if (fabs(rz)>sp) jfz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else jfz[n] += 0;
     			   }
   else                    {
      afx[n] += 0;
      jfx[n] += 0;
      afy[n] += 0;
      jfy[n] += 0;
      afz[n] += 0;
      jfz[n] += 0;
                           }
			   }
      			   }

//high order correction
   for(int n=0; n<N; n++){
      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*jx[n] + 2*jfx[n]) )/(pow(ts, 2)+1e-7);
      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*jy[n] + 2*jfy[n]) )/(pow(ts, 2)+1e-7);
      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*jz[n] + 2*jfz[n]) )/(pow(ts, 2)+1e-7);
      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(jx[n] + jfx[n]) )/(pow(ts, 3)+1e-7);
      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(jy[n] + jfy[n]) )/(pow(ts, 3)+1e-7);
      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(jz[n] + jfz[n]) )/(pow(ts, 3)+1e-7);
                         }

//final acceleration, position and velocity
  for (int n=0; n<N; n++){
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
/*
   for (int n=0; n<N; n++){
      aaf[n]  = sqrt((pow(ax[n], 2) + pow(ay[n], 2) + pow(az[n], 2)));
      ajf[n] = sqrt((pow(jx[n], 2) + pow(jy[n], 2) + pow(jz[n], 2)));
      aa3[n]  = sqrt((pow(a3x[n], 2) + pow(a3y[n], 2) + pow(a3z[n], 2)));
      aaf2[n] = sqrt((pow(af2x[n], 2) + pow(af2y[n], 2) + pow(af2z[n], 2)));
      dt[n]   = sqrt(0.01*((aaf[n]*aaf2[n] + pow(ajf[n], 2))/(ajf[n]*aa3[n] + pow(aaf2[n], 2) + 1e-7)));
			  }
   ts = dt[0];
   for (int i=0; i<N; i++){
   if  ( dt[i] < ts )     {
       ts = dt[i];
       if ( ts < 1e-6 ) ts = 1e-7;
			  }
                          }
*/
}

void hermiteMPI( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G )
{
   int NRank, MyRank;

   MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );

//parameters needed
float dt[N];  //timestep based on all particles' properties respectively.
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
const int Tag = 123;
const int Count = 1;
double SendBuf;
double RecvBuf;
double sp = 0.5;

//parameter initialization
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

//x direction
   if ( MyRank == 0 ){
   for (int n=0; n<N; n++){//first position, velocity update (drift and kick)
      x[n] += ( ts*vx[n] + pow(ts, 2)*ax[n]/2 + pow(ts, 3)*jx[n]/6 );
      vx[n] += ( ts*ax[n] + pow(ts, 2)*jx[n]/2 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      if ( fabs(rx)>sp ) afx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else afx[n] += 0;
      rvx = vx[j] - vx[n];
      if ( fabs(rx)>sp ) jfx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
	else jfx[n] += 0;
                          }
   else                   {
      afx[n] += 0;
      jfx[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*jx[n] + 2*jfx[n]) )/(pow(ts, 2)+1e-7);
      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(jx[n] + jfx[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jx[n] = jfx[n];
      x[n] += ( pow(ts, 4)*a2x[n]/24 + pow(ts, 5)*a3x[n]/120 );
      vx[n] += ( pow(ts, 3)*a2x[n]/6 + pow(ts, 4)*a3x[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      SendBuf = x[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      y[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = x[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      z[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = vx[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vy[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = vx[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vz[n] = RecvBuf;
                          }
		     }

// y direction
   if ( MyRank == 1 ){
   for (int n=0; n<N; n++){//first position, velocity update (drift and kick)
      y[n] += ( ts*vy[n] + pow(ts, 2)*ay[n]/2 + pow(ts, 3)*jy[n]/6 );
      vy[n] += ( ts*ay[n] + pow(ts, 2)*jy[n]/2 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      ry = y[j]-y[n];
      if (fabs(ry)>sp) afy[n] += G*M[j]*ry/fabs(pow(ry, 3));
	else afy[n] += 0;
      rvy = vy[j] - vy[n];
      if (fabs(ry)>sp) jfy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
                          }
   else                   {
      afy[n] += 0;
      jfy[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*jy[n] + 2*jfy[n]) )/(pow(ts, 2)+1e-7);
      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(jy[n] + jfy[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jy[n] = jfy[n];
      y[n] += ( pow(ts, 4)*a2y[n]/24 + pow(ts, 5)*a3y[n]/120 );
      vy[n] += ( pow(ts, 3)*a2y[n]/6 + pow(ts, 4)*a3y[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      x[n] = RecvBuf;
      SendBuf = y[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      SendBuf = y[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      z[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vx[n] = RecvBuf;
      SendBuf = vy[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      SendBuf = vy[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vz[n] = RecvBuf;
                          }
                     }

//z direction
   if ( MyRank == 2 ){
   for (int n=0; n<N; n++){//first position, velocity update (drift and kick)
      z[n] += ( ts*vz[n] + pow(ts, 2)*az[n]/2 + pow(ts, 3)*jz[n]/6 );
      vz[n] += ( ts*az[n] + pow(ts, 2)*jz[n]/2 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rz = z[j]-z[n];
      if ( fabs(rz)>sp ) afz[n] += G*M[j]*rz/fabs(pow(rz, 3));
      rvz = vz[j] - vz[n];
      if ( fabs(rz)>sp ) jfz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
                          }
   else                   {
      afz[n] += 0;
      jfz[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*jz[n] + 2*jfz[n]) )/(pow(ts, 2)+1e-7);
      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(jz[n] + jfz[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jz[n] = jfz[n];
      z[n] += ( pow(ts, 4)*a2z[n]/24 + pow(ts, 5)*a3z[n]/120 );
      vz[n] += ( pow(ts, 3)*a2z[n]/6 + pow(ts, 4)*a3z[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      x[n] = RecvBuf;
      SendBuf = z[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      y[n] = RecvBuf;
      SendBuf = z[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vx[n] = RecvBuf;
      SendBuf = vz[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vy[n] = RecvBuf;
      SendBuf = vz[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
                          }
                     }
   MPI_Barrier(MPI_COMM_WORLD);
}

void herMPIDKD( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G )
{
   int NRank, MyRank;

   MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );

//parameters needed
float dt[N];         //timestep based on all particles' properties respectively.
double rx;
double ry;
double rz;
double rvx;
double rvy;
double rvz;
double SendBuf;
double RecvBuf;
const int Tag = 123;
double sp = 0.1;
double ahx[N];
double ahy[N];
double ahz[N];
double afx[N];//f stand for final, the predicted acceleration at time t+ts
double afy[N];
double afz[N];
double jhx[N];
double jhy[N];
double jhz[N];
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
double hts = 0.5*ts;
//parameter initialization
   for (int i=0; i<N; i++){ //initialize value to zero
      aaf[i] = 0;
      aaf2[i] = 0;
      ajf[i] = 0;
      aa3[i] = 0;
      afx[i] = 0;
      afy[i] = 0;
      afz[i] = 0;
      ahx[i] = 0;
      ahy[i] = 0;
      ahz[i] = 0;
      jhx[i] = 0;
      jhy[i] = 0;
      jhz[i] = 0;
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

//x direction
if ( MyRank == 0 ){
//first drift of 0.5ts
   for (int n=0; n<N; n++){
      x[n] += ( hts*vx[n] + pow(hts, 2)*ax[n]/2 + pow(hts, 3)*jx[n]/6 );
                          }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      if (fabs(rx)>sp) ahx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else ahx[n] += 0;
      rvx = vx[j] - vx[n];
      if (fabs(rx)>sp) jhx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
	else jhx[n] += 0;
                          }
   else                   {
      ahx[n] += 0;
      jhx[n] += 0;
                          }
                          }
                          }

//kick and second drift
   for (int n=0; n<N; n++){
   vx[n] += ( ts*ahx[n] + pow(ts, 2)*jhx[n]/2 );
   x[n] += ( hts*vx[n] + pow(hts, 2)*ahx[n]/2 + pow(hts, 3)*jhx[n]/6 );
                          }

//predict a and jerk at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      if (fabs(rx)>sp) afx[n] += G*M[j]*rx/fabs(pow(rx, 3));
	else afx[n] += 0;
      rvx = vx[j] - vx[n];
      if (fabs(rx)>sp) jfx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
	else jfx[n] += 0;
                          }
   else                   {
      afx[n] += 0;
      jfx[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*jx[n] + 2*jfx[n]) )/(pow(ts, 2)+1e-7);
      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(jx[n] + jfx[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jx[n] = jfx[n];
      x[n] += ( pow(ts, 4)*a2x[n]/24 + pow(ts, 5)*a3x[n]/120 );
      vx[n] += ( pow(ts, 3)*a2x[n]/6 + pow(ts, 4)*a3x[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      SendBuf = x[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      y[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = x[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      z[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = vx[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vy[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = vx[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vz[n] = RecvBuf;
                          }

                     }

// y direction
   if ( MyRank == 1 ){
//first drift of 0.5ts
   for (int n=0; n<N; n++){
      y[n] += ( hts*vy[n] + pow(hts, 2)*ay[n]/2 + pow(hts, 3)*jy[n]/6 );
                          }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      ry = y[j]-y[n];
      if ( fabs(ry)>sp ) ahy[n] += G*M[j]*ry/fabs(pow(ry, 3));
	else ahy[n] += 0;
      rvy = vy[j] - vy[n];
      if ( fabs(ry)>sp ) jhy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
	else jhy[n] += 0;
                          }
   else                   {
      ahy[n] += 0;
      jhy[n] += 0;
                          }
                          }
                          }

//kick and second drift
   for (int n=0; n<N; n++){
      vy[n] += ( ts*ahy[n] + pow(ts, 2)*jhy[n]/2 );
      y[n] += ( hts*vy[n] + pow(hts, 2)*ahy[n]/2 + pow(hts, 3)*jhy[n]/6 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      ry = y[j]-y[n];
      if ( fabs(ry)>sp ) afy[n] += G*M[j]*ry/fabs(pow(ry, 3));
	else afy[n] += 0;
      rvy = vy[j] - vy[n];
      if ( fabs(ry)>sp ) jfy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else jfy[n] += 0;
                          }
   else                   {
      afy[n] += 0;
      jfy[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*jy[n] + 2*jfy[n]) )/(pow(ts, 2)+1e-7);
      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(jy[n] + jfy[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jy[n] = jfy[n];
      y[n] += ( pow(ts, 4)*a2y[n]/24 + pow(ts, 5)*a3y[n]/120 );
      vy[n] += ( pow(ts, 3)*a2y[n]/6 + pow(ts, 4)*a3y[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      x[n] = RecvBuf;
      SendBuf = y[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      SendBuf = y[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      z[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vx[n] = RecvBuf;
      SendBuf = vy[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      SendBuf = vy[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vz[n] = RecvBuf;
                          }
                     }

//z direction
   if ( MyRank == 2 ){
//first drift of 0.5ts
   for (int n=0; n<N; n++){
      z[n] += ( hts*vz[n] + pow(hts, 2)*az[n]/2 + pow(hts, 3)*jz[n]/6 );
                          }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rz = z[j]-z[n];
      if ( fabs(rz)>sp ) ahz[n] += G*M[j]*rz/fabs(pow(rz, 3));
	else ahz[n] += 0;
      rvz = vz[j] - vz[n];
      if ( fabs(rz)>sp ) jhz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else jhz[n] += 0;
                          }
   else                   {
      ahz[n] += 0;
      jhz[n] += 0;
                          }
                          }
                          }

//kick and second drift
   for (int n=0; n<N; n++){
      vz[n] += ( ts*ahz[n] + pow(ts, 2)*jhz[n]/2 );
      z[n] += ( hts*vz[n] + pow(hts, 2)*ahz[n]/2 + pow(hts, 3)*jhz[n]/6 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rz = z[j]-z[n];
      if ( fabs(rz)>sp ) afz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else afz[n] += 0;
      rvz = vz[j] - vz[n];
      if ( fabs(rz)>sp ) jfz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else afz[n] += 0;
                          }
   else                   {
      afz[n] += 0;
      jfz[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*jz[n] + 2*jfz[n]) )/(pow(ts, 2)+1e-7);
      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(jz[n] + jfz[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jz[n] = jfz[n];
      z[n] += ( pow(ts, 4)*a2z[n]/24 + pow(ts, 5)*a3z[n]/120 );
      vz[n] += ( pow(ts, 3)*a2z[n]/6 + pow(ts, 4)*a3z[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      x[n] = RecvBuf;
      SendBuf = z[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      y[n] = RecvBuf;
      SendBuf = z[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vx[n] = RecvBuf;
      SendBuf = vz[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vy[n] = RecvBuf;
      SendBuf = vz[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
                          }
                     }

   MPI_Barrier(MPI_COMM_WORLD);
}

void herMPIKDK( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *jx, double *jy, double *jz, double ts, double G )
{
   int NRank, MyRank;

   MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );

//parameters needed
float dt[N];         //timestep based on all particles' properties respectively.
double rx;
double ry;
double rz;
double rvx;
double rvy;
double rvz;
double SendBuf;
double RecvBuf;
const int Tag = 123;
double sp = 0.5;
double ahx[N];
double ahy[N];
double ahz[N];
double afx[N];//f stand for final, the predicted acceleration at time t+ts
double afy[N];
double afz[N];
double jhx[N];
double jhy[N];
double jhz[N];
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
double hts = 0.5*ts;

//parameter initialization
   for (int i=0; i<N; i++){ //initialize value to zero
      aaf[i] = 0;
      aaf2[i] = 0;
      ajf[i] = 0;
      aa3[i] = 0;
      afx[i] = 0;
      afy[i] = 0;
      afz[i] = 0;
      ahx[i] = 0;
      ahy[i] = 0;
      ahz[i] = 0;
      jhx[i] = 0;
      jhy[i] = 0;
      jhz[i] = 0;
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

//x direction
if ( MyRank == 0 ){
//first kick of 0.5ts
   for (int n=0; n<N; n++){
   vx[n] += ( hts*ax[n] + pow(hts, 2)*jx[n]/2 );
                          }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      if (fabs(rx)>sp) ahx[n] += G*M[j]*rx/fabs(pow(rx, 3));
        else ahx[n] += 0;
      rvx = vx[j] - vx[n];
      if (fabs(rx)>sp) jhx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
        else jhx[n] += 0;
                          }
   else                   {
      ahx[n] += 0;
      jhx[n] += 0;
                          }
                          }
                          }

//drift and second kick
   for (int n=0; n<N; n++){
   x[n] += ( ts*vx[n] + pow(ts, 2)*ahx[n]/2 + pow(ts, 3)*jhx[n]/6 );
   vx[n] += ( hts*ahx[n] + pow(hts, 2)*jhx[n]/2 );
                          }

//predict a and jerk at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      if ( fabs(rx)>sp ) afx[n] += G*M[j]*rx/fabs(pow(rx, 3));
        else afx[n] += 0;
      rvx = vx[j] - vx[n];
      if ( fabs(rx)>sp ) jfx[n] += G*M[j]*(rvx/fabs(pow(rx, 3)) + 3*(rvx*rx)*rx/fabs(pow(rx, 5)));
        else afx[n] += 0;
                          }
   else                   {
      afx[n] += 0;
      jfx[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*jx[n] + 2*jfx[n]) )/(pow(ts, 2)+1e-7);
      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(jx[n] + jfx[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jx[n] = jfx[n];
      x[n] += ( pow(ts, 4)*a2x[n]/24 + pow(ts, 5)*a3x[n]/120 );
      vx[n] += ( pow(ts, 3)*a2x[n]/6 + pow(ts, 4)*a3x[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      SendBuf = x[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      y[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = x[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      z[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = vx[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vy[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      SendBuf = vx[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vz[n] = RecvBuf;
                          }

                  }

//y direction
if ( MyRank == 1 ){
//first kick of 0.5ts
   for (int n=0; n<N; n++){
   vy[n] += ( hts*ay[n] + pow(hts, 2)*jy[n]/2 );
                          }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      ry = y[j]-y[n];
      if (fabs(ry)>sp) ahy[n] += G*M[j]*ry/fabs(pow(ry, 3));
        else ahy[n] += 0;
      rvy = vy[j] - vy[n];
      if (fabs(ry)>sp) jhy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else jhy[n] += 0;
                          }
   else                   {
      ahy[n] += 0;
      jhy[n] += 0;
                          }
                          }
                          }

//drift and second kick
   for (int n=0; n<N; n++){
   y[n] += ( ts*vy[n] + pow(ts, 2)*ahy[n]/2 + pow(ts, 3)*jhy[n]/6 );
   vy[n] += ( hts*ahy[n] + pow(hts, 2)*jhy[n]/2 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      ry = y[j]-y[n];
      if ( fabs(ry)>sp ) afy[n] += G*M[j]*ry/fabs(pow(ry, 3));
        else afy[n] += 0;
      rvy = vy[j] - vy[n];
      if ( fabs(ry)>sp ) jfy[n] += G*M[j]*(rvy/fabs(pow(ry, 3)) + 3*(rvy*ry)*ry/fabs(pow(ry, 5)));
        else afy[n] += 0;
                          }
   else                   {
      afy[n] += 0;
      jfy[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*jy[n] + 2*jfy[n]) )/(pow(ts, 2)+1e-7);
      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(jy[n] + jfy[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jy[n] = jfy[n];
      y[n] += ( pow(ts, 4)*a2y[n]/24 + pow(ts, 5)*a3y[n]/120 );
      vy[n] += ( pow(ts, 3)*a2y[n]/6 + pow(ts, 4)*a3y[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      x[n] = RecvBuf;
      SendBuf = y[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      SendBuf = y[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      z[n] = RecvBuf;
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vx[n] = RecvBuf;
      SendBuf = vy[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      SendBuf = vy[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD );
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 2, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vz[n] = RecvBuf;
                          }
		     }

//z direction
if ( MyRank == 2 ){
//first kick of 0.5ts
   for (int n=0; n<N; n++){
   vz[n] += ( hts*az[n] + pow(hts, 2)*jz[n]/2 );
                          }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rz = z[j]-z[n];
      if (fabs(rz)>sp) ahz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else ahz[n] += 0;
      rvz = vz[j] - vz[n];
      if (fabs(rz)>sp) jhz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else jhz[n] += 0;
                          }
   else                   {
      ahz[n] += 0;
      jhz[n] += 0;
                          }
                          }
                          }

//drift and second kick
   for (int n=0; n<N; n++){
   z[n] += ( ts*vz[n] + pow(ts, 2)*ahz[n]/2 + pow(ts, 3)*jhz[n]/6 );
   vz[n] += ( hts*ahz[n] + pow(hts, 2)*jhz[n]/2 );
                          }

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rz = z[j]-z[n];
      if ( fabs(rz)>sp ) afz[n] += G*M[j]*rz/fabs(pow(rz, 3));
        else afz[n] += 0;
      rvz = vz[j] - vz[n];
      if ( fabs(rz)>sp ) jfz[n] += G*M[j]*(rvz/fabs(pow(rz, 3)) + 3*(rvz*rz)*rz/fabs(pow(rz, 5)));
        else afz[n] += 0;
                          }
   else                   {
      afz[n] += 0;
      jfz[n] += 0;
                          }
                          }
                          }

//high order correction
   for (int n=0; n<N; n++){
      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*jz[n] + 2*jfz[n]) )/(pow(ts, 2)+1e-7);
      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(jz[n] + jfz[n]) )/(pow(ts, 3)+1e-7);
                          }

//final position, velocity and new jerk
   for (int n=0; n<N; n++){
      jz[n] = jfz[n];
      z[n] += ( pow(ts, 4)*a2z[n]/24 + pow(ts, 5)*a3z[n]/120 );
      vz[n] += ( pow(ts, 3)*a2z[n]/6 + pow(ts, 4)*a3z[n]/24 );
                          }
//update to other ranks
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      x[n] = RecvBuf;
      SendBuf = z[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      y[n] = RecvBuf;
      SendBuf = z[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vx[n] = RecvBuf;
      SendBuf = vz[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 0, Tag, MPI_COMM_WORLD );
                          }
   for (int n=0; n<N; n++){
      MPI_Recv( &RecvBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      vy[n] = RecvBuf;
      SendBuf = vz[n];
      MPI_Send( &SendBuf, 1, MPI_DOUBLE, 1, Tag, MPI_COMM_WORLD );
                          }
                  }
} 
