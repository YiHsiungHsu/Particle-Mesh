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
void Potential( double *rho, double *phi, double G, int BC, int GN, double gs )
{
  if(BC == 0)                                       //period BC
  {
    //fft rho to rhok
    fftw_complex *rhok;
    rhok = (fftw_complex*) fftw_malloc( GN*GN*(GN/2+1) * sizeof(fftw_complex) );
    fftw_plan fft;
    fft = fftw_plan_dft_r2c_3d( GN, GN, GN, rho, rhok, FFTW_ESTIMATE);
    fftw_execute(fft);
    for(int i = 0; i < GN; i++){
    for(int j = 0; j < GN; j++){
    for(int k = 0; k < (GN/2+1); k++){
   // printf("mkr[%2d][%2d][%2d] = %5.5f\n",i,j,k, rhok[k+(GN)*(j+GN*i)][0]);
   // printf("mki = %5.5f\n", rhok[k+(GN)*(j+GN*i)][1]);
    }}}
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
