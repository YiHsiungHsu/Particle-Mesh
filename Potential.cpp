//---------------------------------------------------------------------------------------------
////Function    :  poisson  
////Description :  use the DFT for self-gravity to solve the poisson solver to get the potential
////
////Note        :  it can be optmize at convolution part
////
////
////
////
////Parameter   :  rho = density of particle 
////               phi = potential of particle gravity
////            
////
////Return      : phi
////--------------------------------------------------------------------------------------------
void Possion( double *rho, double *phi )
{
  if(BC == 0)                                       //period BC
  {
    //fft rho to rhok
    fftw_complex *rhok;
    rhok = (fftw_complex*) fftw_malloc( N*N*(N/2+1) * sizeof(fftw_complex) );
    fftw_plan fft;
    fft = fftw_plan_dft_r2c_3d( N, N, N, rho, rhok, FFTW_ESTIMATE);
    fftw_execute(fft);
    //calculat potential phi =-4*M_PI*G/(kx**2+ky**2+kz**2)
    //need normailze with 1/N*N*N
    double _n;
    _n = 1 / (N*N*N);                                //normalize factor
    fftw_complex *phik;
    phik = (fftw_complex*) fftw_malloc( N*N*(N/2+1) * sizeof(fftw_complex) );
    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    for(int k = 0; k < (N/2+1); k++)
    {
        double _k;
        double kx,ky,kz;
        if (i > N/2)
        {   
           kx = N-i;
        }
        else
        {
           kx = i;
        }
        if (j > N/2)
        {
           ky = N-j;
        }
        else
        {
           ky = j;
        }
        kz = k;
        _k = -1/((kx*kx)+(ky*ky)+(kz*kz));
        phik[k+(N/2+1)*(j+N*i)][0] = 4*M_PI*G*_k*_n*rhok[k+(N/2+1)*(j+N*i)][0];  //real part
        phik[k+(N/2+1)*(j+N*i)][1] = 4*M_PI*G*_k*_n*rhok[k+(N/2+1)*(j+N*i)][1];  //imagine part
    }
    //DC = 0
    phik[0][0] = 0; 
    phik[0][1] = 0;
    //fft phik to phi
    fftw_plan ifft;
    ifft = fftw_plan_dft_c2r_3d( N, N, N, phik, phi, FFTW_ESTIMATE);  
    fftw_execute(ifft);
    fftw_destroy_plan(fft)
    fftw_destroy_plan(ifft)
    fftw_free(rhok)
    fftw_free(phik)
  }  
  if (BC == 1)           //isolated boundary
  {
     //zero padding M
     double *zM;         //zero padding M     
     zM = (double*) fftw_malloc( (2*N)*(2*N)*(2*N) * sizeof(double) );
     for (int i = 0; i < 2*N; i++)
     for (int j = 0; j < 2*N; j++)
     for (int k = 0; k < 2*N; k++)
     {
         zM[k+(2*N)*(j+2*N*i)] = 0.0;
     }
     for (int i = 0; i < N; i++)
     for (int j = 0; j < N; j++)
     for (int k = 0; k < N; k++)
     {
         zM[k+(2*N)*(j+2*N*i)] = rho[k+N*(j+N*i)]*dx*dx*dx;
     }
     double *dgf;        //discrete Green's function
     dgf = (double*) fftw_malloc( (2*N)*(2*N)*(2*N) * sizeof(double) );   // dgf = -1*/R
     for (int i = 0; i < 2*N; i++)
     for (int j = 0; j < 2*N; j++)
     for (int k = 0; k < 2*N; k++)
     {
         double nx,ny,nz;
         if (i > N )
            {
               nx = 2*N - i;
            }
         else
            {
               nx = i;
            }
         if (i > N )
            {
               ny = 2*N - j;
            }
         else
            {
               ny = j;
            }

         if (i > N )
            {
               nz = 2*N - k;
            }
         else
            {
               nz = k;
            }
         double _R;
         _R = -1 / dx*sqrt(nx*nx + ny*ny + nz*nz);
         if (i == 0 && j == 0 && k == 0)
            {
               dgf[k+(2*N)*(j+2*N*i)] = 0.0;
            }
         else
            {
                dgf[k+(2*N)*(j+2*N*i)] = _R;
            }
      }   
      //  FFT
      fftw_complex zMk,dftk;
      zMK = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );
      dftk = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );
      fftw_plan fftM, fftR;
      fftM = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, zM, zMk, FFTW_ESTIMATE );
      fftR = fftw_plan_dft_r2c_3d( 2*N, 2*N, 2*N, dft, dftk, FFTW_ESTIMATE );
      fftw_execute(fftM);
      fftw_execute(fftR);
      // convolution to get phi  ( a+bi * c+di )
      fftw_complex conk,phik;
      conk = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );
      phik = (fftw_complex*) fftw_malloc( (2*N)*(2*N)*(N+1) * sizeof(fftw_complex) );
      for (int i = 0; i < 2*N; i++)
      for (int j = 0; j < 2*N; j++)
      for (int k = 0; k < N+1; k++)
      {
          conk[k+(N+1)*(j+2*N*i)][0] = (zMk[k+(N+1)*(j+2*N*i)][0] * dftk[k+(N+1)*(j+2*N*i)][0]) - (zMk[k+(N+1)*(j+2*N*i)][1] * dftk[k+(N+1)*(j+2*N*i)][1]);  // real part 
          conk[k+(N+1)*(j+2*N*i)][1] = (zMk[k+(N+1)*(j+2*N*i)][0] * dftk[k+(N+1)*(j+2*N*i)][1]) + (zMk[k+(N+1)*(j+2*N*i)][1] * dftk[k+(N+1)*(j+2*N*i)][0]);  // imagine part
      }
      double _n;
      _n = 1/(2*N*2*N*2*N );      //normailize factor
      for (int i = 0; i < 2*N; i++)
      for (int j = 0; j < 2*N; j++)
      for (int k = 0; k < N+1; k++)
      {
          phik[k+(N+1)*(j+2*N*i)][0] = G*_n*conk[k+(N+1)*(j+2*N*i)][0];  //real part
          phik[k+(N+1)*(j+2*N*i)][1] = G*_n*conk[k+(N+1)*(j+2*N*i)][1];  //imagine part
      }
      double *_phi; 
      _phi = (double*) fftw_malloc( (2*N)(2*N)(2*N) * sizeof(double) );
      fftw_plan ifft;
      ifft = fftw_plan_dft_c2r_3d( 2*N, 2*N, 2*N, phik, _phi, FFTW_ESTIMATE );
      fftw_execute(ifft);
      for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
      for (int k = 0; k < N; k++)
      {
         phi[k+N*(j+N*i)] = _phi[k+(2*N)*(j+2*N*i)];
      }
      fftw_destroy_plan(fftM);
      fftw_destroy_plan(fftR);
      fftw_destroy_plan(ifft);
      fftw_free(zM);
      fftw_free(zMk);
      fftw_free(dft);
      fftw_free(dftk);
      fftw_free(conk);
      fftw_free(phik);
      fftw_free(_phi);
}


