#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

using namespace std;

int N;                 //number of bodies simulating
double position[N][3]; //positions of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double velocity[N][3]; //velocity of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double mass[N];        //mass of bodies in 3-dimension

void hermite( const int N, double *M, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *ax, double *ay, double *az, double *dax, double *day, double *daz, double gs, const int GN  );
/*
//initial acceleration calculate
   double acceleration[N][3];
   float G = 6.67e-11;
   float r[3];
   float sp = 0.01;//softening parameter (unknown what to set yet)
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
	 for (int k=0; k<3; k++){
	    r[k] = abs(position[i][k]-position[j][k]);
	    acceleration[i][k] += G*mass[j]*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5));
				}
			     }
			  }

//initial acceleration derivative of time calculate
   double da[N][3];
   float rv[3];//relative speed between bodies
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         for (int k=0; k<3; k++){
            r[k] = abs(position[i][k] - position[j][k]);
	    rv[k] = abs(velocity[i][k] - velocity[j][k]);
            da[i][k] += G*mass[j]*(rv[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)) + 3*(rv[k]*r[k])*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)));
                                }
                             }
                          }

//first timestep. Use the minima timestep as the global timestep
   float t = 0;//time passed in simulation
   float dt[N];//timestep based on all particles' properties respectively.
   float tmin;//save the minimum timestep in this variable. All particles follow this timestep.
   float timescale;//the total time of simulation we want(to be given)
   for (int i=0; i<N; i++){
      dt[i] = 0.01*( abs(pow((pow(acceleration[i][0], 2) + pow(acceleration[i][1], 2) + pow(acceleration[i][2], 2)), 0.5))/abs(pow((pow(da[i][0], 2) + pow(da[i][1], 2) + pow(da[i][2], 2)), 0.5)) );
//      t[i] += dt[i];
			  }
   for (int i=0; i<N; i++){
      tmin = dt[0];
      if ( dt[i] < tmin ) tmin = dt[i];
			  }
   t += tmin;
*/
/*
//from now on, we want the simulation to go on by itself until it reach the time scale for simulation to complete.
while ( t<timescale ){
//predict position and velocity (to be corrected)
   for (int i=0; i<N; i++){
      for (int k=0; k<3; k++){
	 position[i][k] += ( tmin*velocity[i][k] + pow(tmin, 2)*acceleration[i][k]/2 + pow(tmin, 3)*da[i][k]/6 );
	 velocity[i][k] += ( tmin*acceleration[i][k] +  pow(tmin, 2)*da[i][k]/2 )
			     }
			  }

//update acceleration and its derivative (to be corrected)
   double af[N][3];
   double daf[N][3];
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         for (int k=0; k<3; k++){
	    r[k] = abs(position[i][k]-position[j][k]);
            af[i][k] = G*mass[j]*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5));
	    rv[k] = abs(velocity[i][k] - velocity[j][k]);
            daf[i][k] = G*mass[j]*(rv[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)) + 3*(rv[k]*r[k])*r[k]/(pow((pow(r[k],2)+pow(sp, 2)), 1.5)));
                                }
                             }
                          }
   
//high order correction
   double a2[N][3];//second order derivative of acceleration
   double a3[N][3]; //third order derivative of acceleration
   for (int i=0; i<N; i++){
      for (int k=0; k<3; k++){
	 a2[i][k] = ( -6*(acceleration[i][k] - af[i][k]) - 2*tmin(2*da[i][k] + daf[i][k]) )/pow(tmin, 2);
	 a3[i][k] = ( 12*(acceleration[i][k] - af[i][k]) - 6*tmin(da[i][k] + daf[i][k]) )/pow(tmin, 3);
                             }
                          }                                }

//final acceleration, position and velocity
   for (int i=0; i<N; i++){
      for (int k=0; k<3; k++){
//	 acceleration[i][k] += ( tmin*da[i][k] + pow(dt[i], 2)*a2[i][k]/2 + pow(tmin, 3)*a3[i][k]/6 );
	 acceleration[i][k] = af[i][k];
	 da[i][k] = daf[i][k];
	 position[i][k] += ( pow(tmin, 4)*a2[i][k]/24 + pow(tmin, 5)*a3[i][k]/120 );
	 velocity[i][k] += ( pow(tmin, 3)*a2[i][k]/6 + pow(tmin, 4)*a3[i][k]/24 );
                             }
                          }   

//new timestep to next step
   double aaf[N];
   double aaf2[N];
   double adaf[N];
   double aa3[N];
   for (int i=0; i<N; i++){
      aaf[i]  = pow((pow(af[i][0], 2) + pow(af[i][1], 2) + pow(af[i][2], 2)), 0.5);
      adaf[i] = pow((pow(daf[i][0], 2) + pow(daf[i][1], 2) + pow(daf[i][2], 2)), 0.5);
      aa3[i]  = pow((pow(a3[i][0], 2) + pow(a3[i][1], 2) + pow(a3[i][2], 2)), 0.5);
      aaf2[i] = pow((pow((a2[i][0] + dt[i]*a3[i][0]), 2) + pow((a2[i][1] + dt[i]*a3[i][1]), 2) + pow((a2[i][2] + dt[i]*a3[i][2]), 2)), 0.5);
      dt[i]   = pow(((aaf[i]*aaf2[i] + pow(adaf[i], 2))/(adaf[i]*aa3[i] + pow(aaf2[i], 2))), 0.5);
//      t[i]   += dt[i];
			  }
   for (int i=0; i<N; i++){
      tmin = dt[0];
      if ( dt[i] < tmin ) tmin = dt[i];
                          }
   t += tmin;
*/

//if t < given time scale, return to beginning of while loop, start another iteration.
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
   if  (i != j)           {
//
   rx = fabs(x[n]-x[j]);

   ry = fabs(y[n]-y[j]);

   rz = fabs(z[n]-z[j]);

   rex = pow(rx, 2)+pow(sp, 2);
//
   rey = pow(ry, 2)+pow(sp, 2);

   rez = pow(rz, 2)+pow(sp, 2);

   afx[n] += G*M[j]*rx/sqrt(pow(rex, 3));
//
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
