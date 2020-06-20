#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

using namespace std;

void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *dax, double *day, double *daz, double gs, const int GN  );

void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *dax, double *day, double *daz, double gs, const int GN  )
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

double dafx[N];

double dafy[N];

double dafz[N];

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

double adaf[N];

double aa3[N];

double rex = 0;

double rey = 0;

double rez = 0;

   for (int i=0; i<N; i++){

      aaf[i] = 0;

      aaf2[i] = 0;

      adaf[i] = 0;

      aa3[i] = 0;

      afx[i] = 0; //set initial value before calculation to avoid default random number

      afy[i] = 0;

      afz[i] = 0;

      dafx[i] = 0;

      dafy[i] = 0;

      dafz[i] = 0;

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

      x[n] += ( ts*vx[n] + pow(ts, 2)*ax[n]/2 + pow(ts, 3)*dax[n]/6 );

      y[n] += ( ts*vy[n] + pow(ts, 2)*ay[n]/2 + pow(ts, 3)*day[n]/6 );

      z[n] += ( ts*vz[n] + pow(ts, 2)*az[n]/2 + pow(ts, 3)*daz[n]/6 );

      vx[n] += ( ts*ax[n] + pow(ts, 2)*dax[n]/2 );

      vy[n] += ( ts*ay[n] + pow(ts, 2)*day[n]/2 );

      vz[n] += ( ts*az[n] + pow(ts, 2)*daz[n]/2 );
			  }

//predict acceleration and its derivative at time t + ts
//
   for (int n=0; i<N; i++){
//
   for (int j=0; j<N; j++){
//
   if  (n != j)           {
//
   rx = x[j]-x[n];

   ry = y[j]-y[n];

   rz = z[j]-z[n];

   rex = pow(rx, 2)+pow(sp, 2);
//
   rey = pow(ry, 2)+pow(sp, 2);

   rez = pow(rz, 2)+pow(sp, 2);

   afx[n] += G*M[j]*rx/sqrt(pow(rex, 3));
//
   afy[n] += G*M[j]*ry/sqrt(pow(rey, 3));

   afz[n] += G*M[j]*rz/sqrt(pow(rez, 3));

   rvx = vx[j] - vx[n];

   rvy = vy[j] - vy[n];

   rvz = vz[j] - vz[n];

   dafx[n] += G*M[j]*(rvx/sqrt(pow(rex, 3)) + 3*(rvx*rx)*rx/sqrt(pow(rex, 5)));

   dafy[n] += G*M[j]*(rvy/sqrt(pow(rey, 3)) + 3*(rvy*ry)*ry/sqrt(pow(rey, 5)));

   dafz[n] += G*M[j]*(rvz/sqrt(pow(rez, 3)) + 3*(rvz*rz)*rz/sqrt(pow(rez, 5)));
			   }
//high order correction
//
   for (int n=0; n<N; n++){

      a2x[n] = ( 6*(afx[n] - ax[n]) - ts*(4*dax[n] + 2*dafx[n]) )/(pow(ts, 2)+1e-7);

      a2y[n] = ( 6*(afy[n] - ay[n]) - ts*(4*day[n] + 2*dafy[n]) )/(pow(ts, 2)+1e-7);

      a2z[n] = ( 6*(afz[n] - az[n]) - ts*(4*daz[n] + 2*dafz[n]) )/(pow(ts, 2)+1e-7);

      a3x[n] = ( 12*(ax[n] - afx[n]) - 6*ts*(dax[n] + dafx[n]) )/(pow(ts, 3)+1e-7);

      a3y[n] = ( 12*(ay[n] - afy[n]) - 6*ts*(day[n] + dafy[n]) )/(pow(ts, 3)+1e-7);

      a3z[n] = ( 12*(az[n] - afz[n]) - 6*ts*(daz[n] + dafz[n]) )/(pow(ts, 3)+1e-7);
                          }
//final acceleration, position and velocity
//
   for (int n=0; n<N; n++){
//
      ax[n] += ( ts*dax[n] + pow(ts, 2)*a2x[n]/2 + pow(ts, 3)*a3x[n]/6 );

      ay[n] += ( ts*day[n] + pow(ts, 2)*a2y[n]/2 + pow(ts, 3)*a3y[n]/6 );

      az[n] += ( ts*daz[n] + pow(ts, 2)*a2z[n]/2 + pow(ts, 3)*a3z[n]/6 );

      dax[n] = dafx[n];

      day[n] = dafy[n];

      daz[n] = dafz[n];

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

      adaf[n] = sqrt((pow(dax[n], 2) + pow(day[n], 2) + pow(daz[n], 2)));

      aa3[n]  = sqrt((pow(a3x[n], 2) + pow(a3y[n], 2) + pow(a3z[n], 2)));

      aaf2[n] = sqrt((pow(af2x[n], 2) + pow(af2y[n], 2) + pow(af2z[n], 2)));

      dt[n]   = sqrt(0.01*((aaf[n]*aaf2[n] + pow(adaf[n], 2))/(adaf[n]*aa3[n] + pow(aaf2[n], 2) + 1e-7)));
			  }

   ts = dt[0];

   for (int i=0; i<N; i++){

   if  ( dt[i] < tmin )   {

       ts = dt[i];

       if ( tmin < 1e-6 ) tmin = 1e-7;

			  }

                          }
}
