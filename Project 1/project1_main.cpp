#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include <armadillo>

using namespace arma;
using namespace std;

//object for output file
ofstream ofile;

//source term
inline double f(double x){
   return 100.0*exp(-10.0*x);
}

//closed-form solution
inline double exact_solution(double x) {
   return 1.0-(1.0-exp(-10))*x-exp(-10*x);
}

int main(int argc, char *argv[]){

   int n;
   string which_algorithm;
      string filename = "solution.out";
      // We read the number of grid points n and whether to use general, special-case algorithm or LU-decomposition
      if( argc <= 1 ){
          cout << "Bad usage. Input also the dimension n in the same line and 'general', 'special' or 'LU' to choose the algorithm. " << endl;
          exit(1);
      }
      else{
        n = atoi(argv[1]);
        which_algorithm = argv[2];
      }
      cout << "Using " << which_algorithm << " algorithm with " << n << " grid points." << endl;
   string fileout = filename;
   string argument = to_string(n);
   fileout.append(argument);

   // stepsize
   double h = 1.0/(n+1);
   double hh = h*h;

   // Specify elements common for first two algorithms (general and special)----------------------


   // Set up arrays
   double *b_tilde = new double [n+1];
   double *v = new double [n+2];
   double *x = new double[n+2];
   // Dirichlet boundary conditions
   v[0] = v[n+1] = 0.0;
   //Grid points - interval
   x[0] = 0.0;
   x[n+1] = 1.0;
   for (int i = 0; i <= n; i++){
      x[i]= i*h;
      b_tilde[i] = hh*f(i*h);
      }


    if (which_algorithm=="general"){
   //--------------GENERAL AGORITHM--------------------------------

   clock_t start1, finish1;

   // Set up arrays
   double *a = new double [n+1];
   double *b = new double [n+1];
   double *c = new double [n+1];
   a[0] = a[n] = 0;
   b[0] = 0; b[n] = 2.0;
   c[0] = c[n] = 0;
   for (int i = 1; i < n; i++){
      a[i] = -1.0;
      b[i] = 2.0;
      c[i] = -1.0;
      }


    //create modified coefficients
    double *c_prim = new double [n+1];
    double *b_tilde_prim = new double [n+1];
    c_prim[0] = b_tilde_prim[0] = 0.0;

    start1 = clock();

    // Forward substitution
    for (int i = 1; i <= n; i++){
       double factor = b[i]-a[i-1]*c_prim[i-1]; //m-factor
       c_prim[i] = c[i]/factor;
       b_tilde_prim[i] = (b_tilde[i] - a[i-1]*b_tilde_prim[i-1])/factor;
    }

    // Backward substitution
    for (int i = n; i > 0; i--){
       v[i] = b_tilde_prim[i]-c_prim[i]*v[i+1];
       }

    finish1 = clock();
    double time1;
    time1 = (double(finish1 - start1)/CLOCKS_PER_SEC);
    cout << "Time elapsed for general algorithm: " << time1 << " s" << endl;
    delete [] a; delete [] b; delete [] c; delete [] c_prim; delete [] b_tilde_prim;

    } else if (which_algorithm=="special") {

    // SPECIAL CASE ALGORITHM---------------------------------------------------

    clock_t start2, finish2;

    //create modified coefficients
    double *b_tilde_prim_special = new double [n+1];
    b_tilde_prim_special[0] = 0.0;
    double *factor_special = new double [n+1];
    factor_special[0] = 0.0;
    for (int i = 1; i <= n; i++){
       factor_special[i] = 2.0 + (-i+1.0)/i;
       }
    double *factor_special2 = new double [n+1];
    for (int i = n; i > 0; i--){
       factor_special2[i] = -i/(i+1.0);
       }
    factor_special2[0] = 0.0;

    start2 = clock();
    //Forward substitution
    for (int i = 1; i <= n; i++){
       b_tilde_prim_special[i] = (b_tilde[i] + b_tilde_prim_special[i-1])/factor_special[i];

    }

    // Backward substitution
    for (int i = n; i > 0; i--){
       v[i] = b_tilde_prim_special[i]-factor_special2[i]*v[i+1];
    }

    finish2 = clock();
    double time2;
    time2 = (double(finish2 - start2)/CLOCKS_PER_SEC);
    cout << "Time elapsed for special-case algorithm: " << time2 << " s" << endl;
    delete [] b_tilde_prim_special;
    }
    else if  (which_algorithm=="LU") {

    // LU-DECOMPOSITION---------------

    clock_t start3, finish3;
    // Set up necessary vectors and arrays
    //n = n-1;
    mat A = zeros<mat>(n,n);
    vec b_tilde_vector(n);
    vec x_vector(n);
    A(0,0) = 2.0;
    A(0,1) = -1;
    x_vector(0) = h;
    b_tilde_vector(0) =  hh*f(x_vector(0));
    x_vector(n-1) = x_vector(0)+(n-1)*h;
    b_tilde_vector(n-1) = hh*f(x_vector(n-1));
    for (int i = 1; i < n-1; i++){
        x_vector(i) = x_vector(i-1)+h;
        b_tilde_vector(i) = hh*f(x_vector(i));
        A(i,i-1)  = -1.0;
        A(i,i)    = 2.0;
        A(i,i+1)  = -1.0;
        }
    A(n-1,n-1) = 2.0;
    A(n-2,n-1) = -1.0;
    A(n-1,n-2) = -1.0;

    start3 = clock();
    //perform LU-decomposition
    mat L, U;
    lu(L,U,A);

    // solve Av = LUv = b_tilde in two steps
    vec v_vector0 = solve(L,b_tilde_vector);
    vec v_vector = solve(U,v_vector0);
    finish3 = clock();
    double time3;
    time3 = (double(finish3 - start3)/CLOCKS_PER_SEC);
    cout << "Time elapsed for the algorithm with LU-decomposition: " << time3 << " s" << endl;

    ofile.open(fileout);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    //      ofile << "       x-values:             approximation:          exact solution:" <<endl;
    for (int i = 0; i < n;i++) {
       ofile << setw(15) << setprecision(8) << x_vector(i);
       ofile << setw(15) << setprecision(8) << v_vector(i);
       ofile << setw(15) << setprecision(8) << exact_solution(x_vector(i)) <<endl;
     }
     ofile.close();

    }

    if (which_algorithm != "LU"){
    //Compute relative error
    double *Relative_error = new double [n+1];
    double max_Relative_error;
    max_Relative_error = -20.0;
    for (int i = 1; i < n+1;i++) {
       Relative_error[i] = log10(fabs((v[i] - exact_solution(x[i]))/exact_solution(x[i])));
       if (Relative_error[i] > max_Relative_error) { max_Relative_error = Relative_error[i];}
       }

    cout << "Maximal relative error for the step-size log(h) = " << log10(h) << " is " << max_Relative_error << "." << endl;

    ofile.open(fileout);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    //      ofile << "       x-values:             approximation:          exact solution:       relative error:" << endl;
    for (int i = 1; i < n+1;i++) {
         ofile << setw(15) << setprecision(8) << x[i];
         ofile << setw(15) << setprecision(8) << v[i];
         ofile << setw(15) << setprecision(8) << exact_solution(x[i]);
         ofile << setw(15) << setprecision(8) << Relative_error[i] << endl;
      }
    ofile.close();
    delete [] x; delete [] v; delete [] b_tilde;
    }


    return 0;
    }













