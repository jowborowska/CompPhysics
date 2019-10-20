#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include<time.h>
#include <omp.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
#define NUM_THREADS 5

using namespace std;



//Function to integrate, in Cartesian coordinates - used in Gauss-Legendre Quadrature
double Integrand_Cartesian(double x1, double y1, double z1, double x2, double y2, double z2){
       double r1_minus_r2 = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
       double r1 = sqrt(x1*x1+y1*y1+z1*z1);
       double r2 = sqrt(x2*x2+y2*y2+z2*z2);
       double alpha = 2.0;
       if (r1_minus_r2 < ZERO) //account for potential numerical problems
          return 0;
       else
          return exp(-2*alpha*(r1+r2))/r1_minus_r2;
}


//Function to integrate, in Spherical coordinates, used by Gaussian-Laguerre (reduced) and Monte Carlo
double Integrand_Spherical(bool reduced, double r1, double r2, double theta1, double theta2, double phi1, double phi2){
        double alpha = 2.0;
	double cos_beta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
        double r_12_squared = r1*r1+r2*r2-2.0*r1*r2*cos_beta;
	if(r_12_squared < ZERO) //account for potential numerical problems
	    return 0;
        if (reduced == true and r_12_squared > ZERO)//here reduce exp, absorbed in weights for Gauss-Laguerre
            return exp(-(2.0*alpha-1.0)*(r1+r2))*r1*r1*r2*r2*sin(theta1)*sin(theta2)/sqrt(r_12_squared); 
        else //use in Monte Carlo with importance sampling
            return (1./16.)*r1*r1*r2*r2*sin(theta1)*sin(theta2)/sqrt(r_12_squared); 
}       


//  Calculate weights and abcissas for Legendre polynomials
void Legendre(double x1, double x2, double x[], double w[], int n){
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2; // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;  // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {  // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

      do {
         p1 =1.0;
	 p2 =0.0;


	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp; // Newton's method
      } while(fabs(z - z1) > ZERO);

        
      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
}



//  Calculate weights and abcissas for Laguerre polynomials
double gammln( double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void Laguerre(double *x, double *w, int n, double alf)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;
        
	for (i=1;i<=n;i++) { 
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in Laguerre" << endl;
		x[i]=z;
       
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}

//Gauss-Legendre Quadrature
double Gauss_Legendre(int N, double a, double b,  double &CPU_time){

	double *x = new double [N]; //vector containing abcissas
	double *w = new double [N]; //vector containing weights

	Legendre(a,b,x,w,N);

	double integral = 0.0;

        clock_t start, end;
        start=clock();
	for (int g = 0;  g < N; g++){
	for (int h = 0;  h < N; h++){
	for (int i = 0;  i < N; i++){
	for (int j = 0;  j < N; j++){
	for (int k = 0;  k < N; k++){
	for (int l = 0;  l < N; l++){
		integral += w[g]*w[h]*w[i]*w[j]*w[k]*w[l]*Integrand_Cartesian(x[g],x[h],x[i],x[j],x[k],x[l]);
	}}}}}}
        end=clock();
        CPU_time = ((double)end-(double)start)/CLOCKS_PER_SEC;

        // free the space in memory
	delete [] x;
	delete [] w;

        return integral;
}

//Improved Gauss-Quadrature, using Laguerre poly for r1,r2 and Legendre poly for angles
double Gauss_improved(int N,  double &CPU_time){
        double pi = 3.14159265359;

	double *x_r = new double [N+1];
	double *w_r = new double [N+1];
	double *x_theta = new double [N];
	double *w_theta = new double [N];
	double *x_phi = new double [N];
	double *w_phi = new double [N];
        
        Laguerre(x_r,w_r,N,0.0); //Laguerre for radial part, alf = 0
	Legendre(0,pi,x_theta,w_theta,N); //Legendre for theta, limit 0 to pi
	Legendre(0,2.0*pi,x_phi,w_phi,N); //Legendre for phi, limit 0 to 2pi
	double integral = 0.0;

	clock_t start, end;
        start=clock();

        for (int g = 1;  g < N+1;  g++){  //r1
	for (int h = 0;  h <  N;  h++){  //theta1
	for (int i = 0;  i <  N;  i++){  //phi1
	for (int j = 1;  j < N+1;  j++){  //r2
	for (int k = 0;  k <  N;  k++){  //theta2
	for (int l = 0;  l <  N;  l++){  //phi2
              double Integrand = Integrand_Spherical(true, x_r[g],x_r[j],x_theta[h],x_theta[k],x_phi[i],x_phi[l]);
              integral += w_r[g]*w_r[j]*w_theta[h]*w_theta[k]*w_phi[i]*w_phi[l]*Integrand;
	}}}}}}

	end=clock();
        CPU_time = ((double)end-(double)start)/CLOCKS_PER_SEC;

	//free space in memory
	delete [] x_r;
	delete [] w_r;
	delete [] x_theta;
	delete [] w_theta;
	delete [] x_phi;
	delete [] w_phi;

        return integral;
}


// Minimal random number generator, returns a uniform deviate between 0.0 and 1.0
double ran0(long *idum)
{
   long     k;
   double   ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}

//Monte Carlo integration, uniform distribution
void Monte_Carlo_uniform(int N, double a, double b, double  &integral, double  &std, double &CPU_time){

        long idum = -1;
	double cart[6]; //Cartesian coordinates x1,y1,z1, x2,y2,z2
        double Function;
	double Integral_MC = 0.0;
	double sum_variance = 0.0;
        double variance;
	double volume = pow((b-a),6);
	clock_t start, end;
        start=clock();

	for (int i = 1; i <= N; i++){
                for (int j = 0; j < 6; j++) {
		cart[j] = a + (b-a)*ran0(&idum);
		}
   
	    Function = Integrand_Cartesian(cart[0], cart[1], cart[2],cart[3], cart[4], cart[5]);
            Integral_MC += Function;
            sum_variance += Function*Function;
	}

        Integral_MC = Integral_MC/((double) N );
	integral = volume*Integral_MC;
        end=clock();
        CPU_time = ((double)end-(double)start)/CLOCKS_PER_SEC;

        sum_variance = sum_variance/((double) N );
	variance = sum_variance-Integral_MC*Integral_MC;
        std = volume*sqrt(variance/((double) N ));
	
}

//Monte Carlo integration, importance sampling
void Monte_Carlo_Importance_Sampl(int N, double  &integral, double  &std, double &CPU_time){
	
        long idum = -1;
        double pi = 3.14159265359;
        int i;
	double r1, theta1, phi1, r2, theta2, phi2;
        double Function;
	double Integral_MC = 0.0;
	double sum_variance = 0.0;
        double variance;
	double volume = 4*pow(pi,4); 
	
        double start_time = omp_get_wtime();
        #pragma omp parallel for reduction(+:Integral_MC) reduction(+:sum_variance)  private (idum,i,r1,theta1,phi1,r2,theta2,phi2) 
       
	for (i = 1; i <= N; i++){
		r1 = -0.25*log(1.0-ran0(&idum));
		theta1 = pi*ran0(&idum);
                phi1 = 2.0*pi*ran0(&idum);
                r2 = -0.25*log(1.0-ran0(&idum));
		theta2 = pi*ran0(&idum);
		phi2 = 2.0*pi*ran0(&idum);
                
		Function = Integrand_Spherical(false, r1, r2, theta1, theta2, phi1, phi2);
                Integral_MC += Function;
                sum_variance += Function*Function;
                
	}
        Integral_MC = Integral_MC/((double) N );
        integral = volume*Integral_MC;
        
        CPU_time = omp_get_wtime() - start_time;
        
        sum_variance = sum_variance/((double) N );
	variance = sum_variance-Integral_MC*Integral_MC; 
        std = volume*sqrt(variance/((double) N ));
	
}


// Main function
int main(){ 
        omp_set_num_threads(NUM_THREADS); //NUM_THREADS = 5
        double pi = 3.14159265359;
        double Exact_solution = 5.0*pi*pi/(16.0*16.0);
	int N, N_MC;
	double a,b;
	N = 10; //Used in Gaussian Quadrature
        N_MC = 100000000; //Used in Monte Carlo
	a = -3;
	b =  3;
        double Integral_Gauss_Legendre, CPU_time_gauleg;
        double Integral_Gauss_Laguerre, CPU_time_gaulag;
        double Integral_MC_uniform,std_MC_uniform, CPU_time_MC_uniform;
        double Integral_MC_improved, std_MC_improved, CPU_time_MC_improved;
	Integral_Gauss_Legendre = Gauss_Legendre(N,a,b, CPU_time_gauleg);
        Integral_Gauss_Laguerre = Gauss_improved(N, CPU_time_gaulag);
        Monte_Carlo_uniform(N_MC,a,b, Integral_MC_uniform, std_MC_uniform, CPU_time_MC_uniform);
        Monte_Carlo_Importance_Sampl(N_MC, Integral_MC_improved, std_MC_improved, CPU_time_MC_improved);

        cout << "Exact closed form solution: " << Exact_solution << endl;
	cout << "Gauss-Legendre: " << Integral_Gauss_Legendre << " Relative error: " << (Exact_solution-Integral_Gauss_Legendre)/Exact_solution <<" CPU time (sec): "<< CPU_time_gauleg <<endl;
        cout << "Gauss-Laguerre-Legendre: " << Integral_Gauss_Laguerre << " Relative error: " << (Exact_solution-Integral_Gauss_Laguerre)/Exact_solution <<" CPU time (sec): "<< CPU_time_gaulag <<endl;
        cout << "---Monte Carlo, brute force---" << endl; 
        cout<<" CPU time (sec): "<< CPU_time_MC_uniform <<endl;
        cout << " Integral: " << Integral_MC_uniform << " Standard deviation: " <<std_MC_uniform <<endl;
        cout << "---Monte Carlo, importance sampling---" << endl; 
        cout<<" CPU time (sec): "<< CPU_time_MC_improved <<endl;
        cout << " Integral: " << Integral_MC_improved << " Standard deviation: " <<std_MC_improved <<endl;

	return 0;
}
