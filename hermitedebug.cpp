#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

using namespace std;

int main( int argc, char *argv[] )
{
const int N = 2;                 //number of bodies simulating
double position[N][3] = {
                        {10, 0.0, 0.0},
                        {-10, 0.0, 0.0}
                        }; //positions of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double velocity[N][3] = {
                        {-0.1, 0.0, 0.0},
                        {0.1, 0.0, 0.0}
                        }; //velocity of bodies in 3-dimension, [N][0] for x, [N][1] for y, [N][2] for z
double mass[N] = {1e10, 1e10};        //mass of bodies in 3-dimension
printf("position: x1=%lf, y1=%lf, z1=%lf, x2=%lf, y2=%lf, z2=%lf\n", position[0][0], position[0][1], position[0][2], position[1][0], position[1][1], position[1][2]);
double acceleration[N][3];
float G = 6.67e-11;
float r[3];
float sp = 0.01;//softening parameter
   for (int i=0; i<N; i++){
   for (int j=0; j<N; j++){
   for (int k=0; k<3; k++){
   if (i != j)		  {
      r[k] = fabs(position[i][k]-position[j][k]);
      acceleration[i][k] += G*mass[j]*r[k]/sqrt(pow((pow(r[k],2)+pow(sp, 2)), 3));
                          }
			  }
                          }
      printf("acceleration:%lf, %lf, %lf\n",  acceleration[i][0],  acceleration[i][1],  acceleration[i][2]);
      printf("r:%f, %f, %f\n", r[0], r[1], r[2]);
                          }

//initial acceleration derivative of time calculate
   double da[N][3];
   float rv[3];//relative speed between bodies
   for (int i=0; i<N; i++){
   for (int j=0; j<N; j++){
   for (int k=0; k<3; k++){
   if (i != j)            {
      r[k] = fabs(position[i][k] - position[j][k]);
      rv[k] = fabs(velocity[i][k] - velocity[j][k]);
      da[i][k] += G*mass[j]*(rv[k]/sqrt(pow((pow(r[k],2)+pow(sp, 2)), 3)) + 3*(rv[k]*r[k])*r[k]/sqrt(pow((pow(r[k],2)+pow(sp, 2)), 5)));
                          }
			  }
                          }
//      printf("%lf, %lf, %lf\n",  da[i][0],  da[i][1],  da[i][2]);
//      printf("%f, %f, %f\n", rv[0], rv[1], rv[2]);
                          }

//first timestep. Use the minima timestep as the global timestep
float t = 0;         //time passed in simulation
float dt[N];         //timestep based on all particles' properties respectively.
float tmin;          //save the minimum timestep in this variable. All particles follow this timestep.
float timescale = 50;//the total time of simulation we want(to be given)
float aa[N];         //absolute value of acceleration
float ada[N];        //absolute value of acceleration derivative
   for (int i=0; i<N; i++){
      aa[i] = sqrt( pow(acceleration[i][0], 2)+ pow(acceleration[i][1], 2)+pow(acceleration[i][2], 2) );
      ada[i] = sqrt( pow(da[i][0], 2)+ pow(da[i][1], 2)+pow(da[i][2], 2) );
      dt[i] = 0.01*aa[i]/ada[i];
			  }
tmin = dt[0];
   for (int i=0; i<N; i++){
   if ( dt[i] < tmin ) tmin = dt[i];
                          }
   printf("%f\n", tmin);
   t += tmin;

//pre-set some parameter needed later
double af[N][3];
double daf[N][3];
double a2[N][3]; //second order derivative of acceleration
double a3[N][3]; //third order derivative of acceleration
double af2[N][3];   //second derivative of acceleration at t + dt
double aaf[N];   //first a of all below stands for absolute value
double aaf2[N];
double adaf[N];
double aa3[N];
double re[3];

while ( t<timescale )
{
//predict position and velocity (to be corrected)
   for (int i=0; i<N; i++){
   for (int k=0; k<3; k++){
      position[i][k] += ( tmin*velocity[i][k] + pow(tmin, 2)*acceleration[i][k]/2 + pow(tmin, 3)*da[i][k]/6 );
      velocity[i][k] += ( tmin*acceleration[i][k] +  pow(tmin, 2)*da[i][k]/2 );
                          }
//      printf("position[%d]: %lf, %lf, %lf\n", i, position[i][0], position[i][1], position[i][2]);
//      printf("velocity[%d]: %lf, %lf, %lf\n", i, velocity[i][0], velocity[i][1], velocity[i][2]);
                          }

//update acceleration and its derivative
   for (int i=0; i<N; i++){
   for (int j=0; j<N; j++){
   for (int k=0; k<3; k++){
   if  (i != j)           {
      r[k] = fabs(position[i][k]-position[j][k]);
      re[k] = pow(r[k], 2)+pow(sp, 2);
//      printf("re[%d]:%f\n", k, re[k]);
      af[i][k] += G*mass[j]*r[k]/sqrt(pow(re[k], 3));
//      printf("af[%d][%d]:%f\n", i, k, af[i][k]);
//      printf("分母:%f\n", sqrt(pow(re[k], 3)));
      rv[k] = fabs(velocity[i][k] - velocity[j][k]);
//      printf("rv[%d]:%f\n", k, rv[k]);
      daf[i][k] += G*mass[j]*(rv[k]/sqrt(pow(re[k], 3)) + 3*(rv[k]*r[k])*r[k]/sqrt(pow(re[k], 5)));
//      printf("daf[%d][%d]:%f\n", i, k, daf[i][k]);

                          }
			  }
                          }
                          }

//high order correction
   for (int i=0; i<N; i++){
   for (int k=0; k<3; k++){
      a2[i][k] = ( 6*(af[i][k] - acceleration[i][k]) - tmin*(4*da[i][k] + 2*daf[i][k]) )/pow(tmin, 2);
//      printf("deltaa[%d][%d]=%lf  ", i, k, (af[i][k]-acceleration[i][k]));
//      printf("deltada[%d][%d]=%lf  ", i, k, (4*da[i][k]+2*daf[i][k]));
//      printf("a2[%d][%d]=%lf  ", i, k, a2[i][k]);
      a3[i][k] = ( 12*(acceleration[i][k] - af[i][k]) - 6*tmin*(da[i][k] + daf[i][k]) )/pow(tmin, 3);
//      printf("a3[%d][%d]=%lf  ", i, k, a3[i][k]);
                          }
//      printf("\n");
                          }

//final acceleration, position and velocity
   for (int i=0; i<N; i++){
   for (int k=0; k<3; k++){
      acceleration[i][k] = af[i][k];
      da[i][k] = daf[i][k];
      position[i][k] += ( pow(tmin, 4)*a2[i][k]/24 + pow(tmin, 5)*a3[i][k]/120 );
      velocity[i][k] += ( pow(tmin, 3)*a2[i][k]/6 + pow(tmin, 4)*a3[i][k]/24 );
      af2[i][k] = a2[i][k] + tmin*a3[i][k];
                          }
                          }
//check the result
 printf("position: x1=%lf, y1=%lf, z1=%lf, x2=%lf, y2=%lf, z2=%lf\n", position[0][0], position[0][1], position[0][2], position[1][0], position[1][1], position[1][2]);

//new timestep to next step
//double aaf[N];
//double aaf2[N];
//double adaf[N];
//double aa3[N];
   for (int i=0; i<N; i++){
      aaf[i]  = sqrt((pow(af[i][0], 2) + pow(af[i][1], 2) + pow(af[i][2], 2)));
      printf("aaf[%d]:%lf\n", i, aaf[i]);
      adaf[i] = sqrt((pow(daf[i][0], 2) + pow(daf[i][1], 2) + pow(daf[i][2], 2)));
      printf("adaf[%d]:%lf\n", i, adaf[i]);
      aa3[i]  = sqrt((pow(a3[i][0], 2) + pow(a3[i][1], 2) + pow(a3[i][2], 2)));
      printf("aa3[%d]:%lf\n", i, aa3[i]);
      aaf2[i] = sqrt((pow(af2[i][0], 2) + pow(af2[i][1], 2) + pow(af2[i][2], 2)));
      printf("aaf2[%d]:%lf\n", i, aaf2[i]);
      dt[i]   = sqrt(0.01*((aaf[i]*aaf2[i] + pow(adaf[i], 2))/(adaf[i]*aa3[i] + pow(aaf2[i], 2) + 1e-7)));
      printf("dt[%d]:%f\n", i, dt[i]);
			  }
   tmin = dt[0];
   for (int i=0; i<N; i++){
   if  ( dt[i] < tmin )   {
       tmin = dt[i];
       if ( tmin < 1e-6 ) tmin = 1e-7;
			  }
                          }
   printf("%f\n", tmin);
   t += tmin;

}
   return EXIT_SUCCESS;
}
