#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <string>
#include <omp.h>

using namespace  std;
using namespace arma;

ofstream ofile;

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define NUM_THREADS 4

// long periode (> 2 x 10^18) random number generator of L'Ecuyer and Bays-Durham shuffle and added safeguards, from lib.cpp
double ran2(long *idum)
{
   int            j;
   long           k;
   static long    idum2 = 123456789;
   static long    iy=0;
   static long    iv[NTAB];
   double         temp;

   if(*idum <= 0) {
      if(-(*idum) < 1) *idum = 1;
      else             *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ1;
	 *idum = IA1*(*idum - k*IQ1) - k*IR1;
	 if(*idum < 0) *idum +=  IM1;
	 if(j < NTAB)  iv[j]  = *idum;
      }
      iy=iv[0];
   }
   k     = (*idum)/IQ1;
   *idum = IA1*(*idum - k*IQ1) - k*IR1;
   if(*idum < 0) *idum += IM1;
   k     = idum2/IQ2;
   idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
   if(idum2 < 0) idum2 += IM2;
   j     = iy/NDIV;
   iy    = iv[j] - idum2;
   iv[j] = *idum;
   if(iy < 1) iy += IMM1;
   if((temp = AM*iy) > RNMX) return RNMX;
   else return temp;
}


//compute quantities analytically for 2x2 lattice
void Analytical_2times2(double T, double& Z, double& E_expect, double& M_expect, double& M_abs_expect, double& C_V, double& chi){
   double J = 1.0; double k = 1.0;
   double beta = 1.0/(k*T);

   Z = 2.0*exp(-8.0*J*beta) + 2.0*exp(8.0*J*beta) + 12.0;
   E_expect = (16.0*J/Z)*(exp(-8.0*beta*J)-exp(8.0*beta*J));
   double E_squared_expect = (128.0*J*J/Z)*(exp(-8.0*beta*J)+exp(8.0*beta*J));
   C_V = (1.0/(k*T*T))*(E_squared_expect-E_expect*E_expect);
   M_expect = 0.;
   double M_squared_expect = (32.0/Z)*(1.0 + exp(8.0*beta*J));
   chi = beta*(M_squared_expect-M_expect*M_expect);
   M_abs_expect = (8.0/Z)*(2.0 + exp(8.0*beta*J));
}


// inline function for periodic boundary conditions
inline int PB(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}

// function to initialize energy, spin matrix and magnetization
void Initialize(int L, mat &Lattice,  double& E, double& M, int ordered_start){
  long idum = -1;
  // setup square LxL spin matrix (ordered or not) and initial magnetization
  for(int x =0; x < L; x++) {
    for (int y= 0; y < L; y++){
      if (ordered_start == 1){
         Lattice(x,y) = 1.0;// spin orientation for the ground state
      }
      else {
         
         if (ran2(&idum) > 0.5){
            Lattice(x,y) = 1.0;
            }
         else
            {Lattice(x,y) = -1.0;
            }
      } 
      M +=  (double) Lattice(x,y); //sum over all spins for a given configuration
    }
  }
  // setup initial energy
  for(int x=0; x < L; x++) {
    for (int y=0; y < L; y++){
      E -=  (double) Lattice(x,y)*(Lattice(PB(x,L,-1),y) + Lattice(x,PB(y,L,-1)));
    }
  }
}




// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void Metropolis(int L, int MC_cycles, double T, vec &Quantities, int ordered_start){

  long idum;
  idum = -1;  // random starting point
  
  mat Lattice = zeros<mat>(L,L); //lattice of spins

  double E = 0.;  //energy
  double M = 0.;  // magnetization
  Initialize(L, Lattice, E, M, ordered_start);
  
  vec Energy_Difference = zeros<mat>(17); // array for possible energy changes
  for( int de =-8; de <= 8; de+=4){ 
      Energy_Difference(de+8) = exp(-de/T);
   }
   int cycles, x, y;
   
   // Start Monte Carlo cycles  
   for (cycles = 1; cycles <= MC_cycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(x=0; x < L; x++) {
      for (y=0; y < L; y++){
	int ix = (int) (ran2(&idum)*(double)L);
        int iy = (int) (ran2(&idum)*(double)L);
	int deltaE =  2*Lattice(ix,iy)*(Lattice(ix,PB(iy,L,-1)) + Lattice(PB(ix,L,-1),iy) + Lattice(ix,PB(iy,L,1)) + Lattice(PB(ix,L,1),iy));
	if (  ran2(&idum) <= Energy_Difference(deltaE+8) ) {
	  Lattice(ix,iy) *= -1.0;  // flip one spin and accept the new configuration
	  M += (double) 2*Lattice(ix,iy);
	  E += (double) deltaE;
	}
      }
    }
    // Update sums: E, E^2, M, M^2, |M|
    Quantities(0) += E;    
    Quantities(1) += E*E;
    Quantities(2) += M;    
    Quantities(3) += M*M; 
    Quantities(4) += fabs(M);
  }
   
} 


// The Monte Carlo part with the Metropolis algo with sweeps over the lattice for different numbers of MC cycles/ part 4c)
void Metropolis_likely_state(int L, int MC_cycles_max, double T, int ordered_start){

  int N = L*L;
  long idum;
  idum = -1;  // random starting point
  
  mat Lattice = zeros<mat>(L,L); //lattice of spins
  double E = 0.;  //energy
  double M = 0.;  // magnetization
  Initialize(L, Lattice, E, M, ordered_start);
  vec Energy_Difference = zeros<mat>(17); // array for possible energy changes
  for( int de =-8; de <= 8; de+=4){ 
      Energy_Difference(de+8) = exp(-de/T);
   }

  int accepted_conf = 0;
  double E_expect = 0.;
  double E2_expect = 0;
  double Mabs_expect = 0;
 
  // Start Monte Carlo cycles
  for (int cycles = 1; cycles <= MC_cycles_max; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int x =0; x < L; x++) {
      for (int y= 0; y < L; y++){
        //position at random location in the lattice
	int ix = (int) (ran2(&idum)*(double)L);
        int iy = (int) (ran2(&idum)*(double)L);
	int deltaE =  2*Lattice(ix,iy)*(Lattice(ix,PB(iy,L,-1)) + Lattice(PB(ix,L,-1),iy) + Lattice(ix,PB(iy,L,1)) + Lattice(PB(ix,L,1),iy));
	if (  ran2(&idum) <= Energy_Difference(deltaE+8) ) {
	  Lattice(ix,iy) *= -1.0;  // flip one spin and accept the new configuration
          accepted_conf += 1;
	  M += (double) 2*Lattice(ix,iy);
	  E += (double) deltaE;
	}
      }
    }
    // Update sums: E,|M|
    E_expect += E;    
    E2_expect += E*E;
    Mabs_expect += fabs(M);
    
    double factor = 1.0/((double) (cycles));
    double E_expectation = E_expect*factor; //<E>
    double Mabs_expectation = Mabs_expect*factor; //<|M|>
    double E2_expectation = E2_expect*factor; //<E^2>
    double E_variance = E2_expectation- E_expectation*E_expectation;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << cycles;
    ofile << setw(15) << setprecision(8) << E_expectation;
    ofile << setw(15) << setprecision(8) << Mabs_expectation;
    ofile << setw(15) << setprecision(8) << accepted_conf;
    ofile << setw(15) << setprecision(8) << E_variance;
    ofile << setw(15) << setprecision(8) << E << endl; //this is then properly cut in Python program ot start after equilibrium point
  }
} 

void Write_to_file(int L, int MC_cycles, double T, vec Quantities){
  
  int N = L*L; //number of spins
  double factor = 1.0/((double) (MC_cycles));  
  double E_expectation = Quantities(0)*factor; //<E>
  double E2_expectation = Quantities(1)*factor; //<E^2>
  double M_expectation = Quantities(2)*factor; // <M>
  double M2_expectation = Quantities(3)*factor; // <M^2>
  double Mabs_expectation = Quantities(4)*factor; //<|M|>

  // all output values are per spin (divided by N)
  double E_variance = (E2_expectation- E_expectation*E_expectation)/N;
  double M_variance = (M2_expectation - M_expectation*M_expectation)/N; //use Mabs_expectation^2 here for 4e) and M_expectation^2 for 4b)
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << T;
  ofile << setw(15) << setprecision(8) << E_expectation/N;
  ofile << setw(15) << setprecision(8) << E_variance/(T*T); //C_V
  ofile << setw(15) << setprecision(8) << M_expectation/N;
  ofile << setw(15) << setprecision(8) << M_variance/T; //chi
  ofile << setw(15) << setprecision(8) << Mabs_expectation/N << endl;
  
} 


// Main program
int main()
{

  // Print the analytical values for 2x2 lattice
  double T = 1.0;
  double Z_analytical, E_analytical, M_analytical, M_abs_expected, C_V_analytical, chi_analytical;
  Analytical_2times2(T, Z_analytical, E_analytical, M_analytical, M_abs_expected, C_V_analytical, chi_analytical);
  cout << "Analytical values for 2x2 lattice and T = 1 KT/J:" << endl;
  cout << "<E> per spin: " << E_analytical/4.0 << endl;
  cout << "<M> per spin: " << M_analytical/4.0 << endl;
  cout << "<|M|> per spin: " << M_abs_expected/4.0 << endl;
  cout << "C_V per spin: " << C_V_analytical/4.0 << endl;
  cout << "chi per spin: " << chi_analytical/4.0 << endl;

  int L, MC_cycles;
  double T_initial, T_final, dT;
  bool likely_state_computations = false; //<-----------------true to perform part 4c) 

  L = 40; //number of spins in one direction
  MC_cycles = 100000; //number of MC cycles    
  T_initial = 2.1; //initial temperature
  T_final = 2.4; //final temperature
  dT = 0.02; //temperature step
  
  string fileout = "output_data40.tex"; //<------------------------------set the filename here
  ofile.open(fileout);
  
  omp_set_num_threads(NUM_THREADS);
 
  double T_points = (T_final-T_initial)/dT + 1;
  T_points  = (int) T_points;
  vec Temperature = zeros<mat>(T_points+1); // array for possible energy changes
  for( int i=0; i <= T_points; i++){ 
      Temperature(i) = T_initial+i*dT;
   }

  if (likely_state_computations == false){
  ofile << "        " << "  T  " << "           <E>  " << "         C_V  " << "        <M>  " << "            chi  " << "          <|M|>  " <<endl;
  //Loop over temperature, perform computations and write results to file
  //# pragma omp parallel for
  for (int g=0; g <= T_points; g++){
    vec Quantities = zeros<mat>(5);
    Metropolis(L, MC_cycles, Temperature(g), Quantities, 1);
    Write_to_file(L, MC_cycles, Temperature(g), Quantities);
    }

  }
  if (likely_state_computations == true){
  //Likely state computations (4c)
     double T_likely_1 = 1.0;
     double T_likely_2 = 2.4;
     int ordered_start = 0; //1 for ordered starting configuration, 0 for random
     Metropolis_likely_state(20, 1000000, T_likely_1, ordered_start);
  }

  ofile.close();  
  return 0;
}

