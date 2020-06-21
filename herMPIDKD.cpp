#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <math.h>

int main( int argc, char *argv[] )
{
   int NRank, MyRank;

   MPI_Init( &argc, &argv );

   MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );

// this test assumes only three ranks
   if ( NRank != 3 )
   {
      fprintf( stderr, "ERROR: NRank (%d) != 3\n", NRank );
      MPI_Abort( MPI_COMM_WORLD, 1 );
   }

//constant and initial condition
   const int N = 2;

   const int Tag = 123;

   const int Count = 1;

   double SendBuf;

   double RecvBuf;

   double M[N], x[N], y[N], z[N];

   const double t_end = 0.05; // end time

   const double ts = 0.05; //time step size of each step

   const double G = 0.25/M_PI; //(m3 kg-1 s-2)

   double vx[N], vy[N], vz[N], ax[N], ay[N], az[N], jx[N], jy[N], jz[N]; // jerk for x, y, z. j stand for jerk

   double rx, ry, rz, rex, rey, rez;

   double sp = 0.01;

   double t = 0;

   for(int n = 0; n<N; n++){
      M[n] = 1.0;//*(double)rand()/RAND_MAX;// 10 is maxium mass

      x[n] = 1.0 + (double) n;//(double)rand()/RAND_MAX*(GN-1);

      y[n] = 1.0 + (double) n;//(double)rand()/RAND_MAX*(GN-1);

      z[n] = 1.0 + (double) n;//(double)rand()/RAND_MAX*(GN-1);

      vx[n] =0.0;// (double)rand()/RAND_MAX*(GN-1);

      vy[n] =0.0;// (double)rand()/RAND_MAX*(GN-1);

      vz[n] =0.0;// (double)rand()/RAND_MAX*(GN-1);
      
      jx[n] = 0;

      jy[n] = 0;

      jz[n] = 0;
  			    }

   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      ry = y[j]-y[n];
      rz = z[j]-z[n];
      rex = pow(rx, 2)+pow(sp, 2);
      rey = pow(ry, 2)+pow(sp, 2);
      rez = pow(rz, 2)+pow(sp, 2);
      ax[n] += G*M[j]*rx/sqrt(pow(rex, 3));
      ay[n] += G*M[j]*ry/sqrt(pow(rey, 3));
      az[n] += G*M[j]*rz/sqrt(pow(rez, 3));
			  }
   else			  {
      ax[n] += 0;
      ay[n] += 0;
      az[n] += 0;
	   		  }
			  }
			  }
//parameters needed
float dt[N];         //timestep based on all particles' properties respectively.
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

   while(t <= t_end)
   {
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
//first drift of 0.5ts
   for (int n=0; n<N; n++){
      x[n] += ( hts*vx[n] + pow(hts, 2)*ax[n]/2 + pow(hts, 3)*jx[n]/6 );
			  }

//now update a and jerk at t+hts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      rex = pow(rx, 2)+pow(sp, 2);
      ahx[n] += G*M[j]*rx/sqrt(pow(rex, 3));
      rvx = vx[j] - vx[n];
      jhx[n] += G*M[j]*(rvx/sqrt(pow(rex, 3)) + 3*(rvx*rx)*rx/sqrt(pow(rex, 5)));
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

//predict acceleration and its derivative at time t + ts
   for (int n=0; n<N; n++){
   for (int j=0; j<N; j++){
   if  (n != j)           {
      rx = x[j]-x[n];
      rex = pow(rx, 2)+pow(sp, 2);
      afx[n] += G*M[j]*rx/sqrt(pow(rex, 3));
      rvx = vx[j] - vx[n];
      jfx[n] += G*M[j]*(rvx/sqrt(pow(rex, 3)) + 3*(rvx*rx)*rx/sqrt(pow(rex, 5)));
			  }
   else			  {
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
      rey = pow(ry, 2)+pow(sp, 2);
      ahy[n] += G*M[j]*ry/sqrt(pow(rey, 3));
      rvy = vy[j] - vy[n];
      jhy[n] += G*M[j]*(rvy/sqrt(pow(rey, 3)) + 3*(rvy*ry)*ry/sqrt(pow(rey, 5)));
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
      rey = pow(rx, 2)+pow(sp, 2);
      afy[n] += G*M[j]*ry/sqrt(pow(rey, 3));
      rvy = vy[j] - vy[n];
      jfy[n] += G*M[j]*(rvy/sqrt(pow(rey, 3)) + 3*(rvy*ry)*ry/sqrt(pow(rey, 5)));
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
      rez = pow(rz, 2)+pow(sp, 2);
      ahz[n] += G*M[j]*rz/sqrt(pow(rez, 3));
      rvz = vz[j] - vz[n];
      jhz[n] += G*M[j]*(rvz/sqrt(pow(rez, 3)) + 3*(rvz*rz)*rz/sqrt(pow(rez, 5)));
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
      rez = pow(rz, 2)+pow(sp, 2);
      afz[n] += G*M[j]*rz/sqrt(pow(rez, 3));
      rvz = vz[j] - vz[n];
      jfz[n] += G*M[j]*(rvz/sqrt(pow(rez, 3)) + 3*(rvz*rz)*rz/sqrt(pow(rez, 5)));
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
                     }
   t += ts;
   MPI_Barrier(MPI_COMM_WORLD);
   printf("1:%lf, %lf, %lf\n", x[0], y[0], z[0]);
   printf("2:%lf, %lf, %lf\n", x[1], y[1], z[1]);
   }
   MPI_Finalize();

   return EXIT_SUCCESS;
}
