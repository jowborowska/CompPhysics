#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector> 
#include <time.h> 

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

using namespace  std;
using namespace  arma;

//Object for output file
ofstream ofile;

//initialize Hamiltonian matrix, rho-vector and Eigenvectors-matrix
void initialize(int N, double h, mat& Hamiltonian_matrix, vec& rho, mat& v,int num_electrons,double omega_r){
    
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);
    
    for (int i=0; i<N ;i++){
        rho(i)=(i+1)*h;
        }

    for (int i=0;i<N;i++){
        double V = rho(i)*rho(i);
        for (int j=0;j<N;j++){
            
            if(i==j && num_electrons==0){ //the buckling beam problem, 1D
                Hamiltonian_matrix(i,j)=d;
                v(i,j)=1;
            }
            else if (i==j && num_electrons==1){ //one electron, qunatum dots in 3D
                Hamiltonian_matrix(i,j)=d+V;
                v(i,j)=1;
            }
            else if (i==j && num_electrons==2){ //two electrons, quantum dots in 3D
                Hamiltonian_matrix(i,j)=d+omega_r*omega_r*V + 1.0/rho(i);
                v(i,j)=1;
            }
            else if (i==j+1 or i==j-1){
                Hamiltonian_matrix(i,j)= a;
            }
            else{
                Hamiltonian_matrix(i,j)=0;
                v(i,j)=0;
            }}}}

//Use Armadillo to diagonalize the matrix and write exact and analytical solution to an output file
void Diagonalize_with_Armadillo(int N, double h, mat& Hamiltonian_matrix, mat& Eigenvectors, int num_electrons){

     double d = 2.0/(h*h);
     double a = -1.0/(h*h);

     //Diagonalize and obtain eigenvalues
     clock_t start, end;
     start=clock();

     vec Eigenvalues(N);

     eig_sym(Eigenvalues, Eigenvectors, Hamiltonian_matrix);	//eigen decomposition of symmetric/hermitian matrix
 
     end=clock();
     cout<<"Armadillo-diagonalization CPU time [s] : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
     
     //Calculate analytical solutions for buckling beam and write them togehter with numerical ones to a file
     if (num_electrons == 0){
     double pi = acos(-1.0);
     string filename;
     filename = "diagonalize_arma.out";
     ofile.open(filename);
     ofile << setiosflags(ios::showpoint | ios::uppercase);
     //      ofile << "       Analytical:             Numerical:          Abs_Difference:" << endl;
     for(int i = 0; i < N; i++) {
        int j = i+1;
        double analytical = d+2.0*a*cos(j*pi/(N+1));
        ofile << setw(15) << setprecision(8) << analytical;
        ofile << setw(15) << setprecision(8) << Eigenvalues[i];
        ofile << setw(15) << setprecision(8) << fabs(Eigenvalues[i] - analytical) << endl;
      }}}



//Order eigenvalues
vector<double> Get_eigenvalues(mat a,int N){
    vector<double>Eigenvalues;
    for(int i=0;i<N;i++){
        Eigenvalues.push_back(a(i,i));
    }
    sort (Eigenvalues.begin(), Eigenvalues.begin()+N);
    return Eigenvalues;
}



//Order matrix of eigenvectors (so that they correspond to ordered eigenvalues)
mat Get_eigenvectors(mat a, mat v, int N){
    vector<double>eigenvalues = Get_eigenvalues(a,N);
    mat Eigenvectors(N,N);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(a(j,j)==eigenvalues[i]){
                for(int k=0;k<N;k++){
                      Eigenvectors(i,k)=v(k,j);
                }}}}
    return Eigenvectors;
}

//Find non-diagonal matrix element with maximal absolute value
void Find_max_element(mat Hamiltonian_matrix,int& p,int& q,double& a_pq,int N){
    for (int i=0; i<N; i++){
         for (int j=0; j<N; j++){
            if(i!=j && abs(Hamiltonian_matrix(i,j))>= abs(a_pq)){
                a_pq = Hamiltonian_matrix(i,j);
                p = i;
                q = j;
            }}}}

//Perform Jacobi rotations to diagonalize the matrix
int Jacobi_algorithm(int N, double delta, mat& a, mat& v) {
    
    double a_ip=0, a_iq=0, v_ip=0, v_iq=0;
    double tau=0, t_min=0, s=0, c=0;//tau, tan(theta), sin(theta), cos(theta)
    int iterations=1; //count number of iterations
    int p=N-2, q=N-1; //all off-diagonal elements the same at the start - pick a(N-2, N-1) as first maximum
                                
    clock_t start, end;

    double a_pp=a(p,p);
    double a_qq=a(q,q);
    double a_pq=a(p,q);

    start=clock();

    while(abs(a_pq)> delta){
        if(iterations>1){
            a_pq=0;
            Find_max_element(a,p,q,a_pq,N);
        }

        //Calculate sin(theta) and cos(theta)
        a_qq=a(q,q);
        a_pp=a(p,p);
        tau=(a_qq-a_pp)/(2.0*a_pq);
        if(tau>=0){
            t_min = 1.0/(tau+sqrt(1.0+tau*tau));
        }
        else{
            t_min = 1.0/(tau-sqrt(1.0+tau*tau));
        }
        c = 1.0/sqrt(1.0+t_min*t_min);
        s = c*t_min;

        //Overwrite matrix elements in A-matrix and Eigenvalue-vectors in v-matrix
        for(int i=0;i<N;i++){
            if(i!=p && i!=q){
                a_ip=a(i,p);
                a_iq=a(i,q);
                a(i,p)=a_ip*c-a_iq*s;
                a(p,i)=a_ip*c-a_iq*s;
                a(i,q)=a_iq*c+a_ip*s;
                a(q,i)=a_iq*c+a_ip*s;
            }
            v_ip=v(i,p);
            v_iq=v(i,q);
            v(i,p)=c*v_ip-s*v_iq;
            v(i,q)=c*v_iq+s*v_ip;
        }
        a(p,p)=a_pp*c*c-2*a_pq*c*s+a_qq*s*s;
        a(q,q)=a_pp*s*s+2*a_pq*c*s+a_qq*c*c;
        a(p,q)=0;
        a(q,p)=0;

        iterations++;
    }

    end=clock();


    cout<<"Diagonalization took "<< iterations <<" iterations."<<endl;
    cout<<"Jacobi-algorithm CPU time [s] : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;

    return 0;
}




//Write eigenvalues and eigenvectors to an output file (after performed Jacobi algorithm)
void Write_to_file(mat a_mat, mat v, int N, int num_electrons, double h){
     vec Eigenvalues = Get_eigenvalues(a_mat, N);
     mat Eigenvectors = Get_eigenvectors(a_mat, v, N);
     
     double d = 2.0/(h*h);
     double a = -1.0/(h*h);
     
     if (num_electrons == 0){ //for buckling beam, write also analytical, numerical and error into a file
     double pi = acos(-1.0);
     string filename;
     filename = "jacobi_eigenvalues_beam_analytical.out";
     ofile.open(filename);
     ofile << setiosflags(ios::showpoint | ios::uppercase);
     //      ofile << "       Analytical:             Numerical:          Abs_Difference:" << endl;
     for(int i = 0; i < N; i++) {
        int j = i+1;
        double analytical = d+2.0*a*cos(j*pi/(N+1));
        ofile << setw(15) << setprecision(8) << analytical;
        ofile << setw(15) << setprecision(8) << Eigenvalues[i];
        ofile << setw(15) << setprecision(8) << fabs(Eigenvalues[i] - analytical) << endl;
      }}
     string filename1;
     filename1 = "jacobi_eigen.out";
     ofile.open(filename1);
     ofile << setiosflags(ios::showpoint | ios::uppercase);
     //    ofile << "       Eigenvalue:             Eigenvector: " << endl;
     for(int i = 0; i < N; i++) {
        ofile << setw(15) << setprecision(8) << Eigenvalues(i) << endl;
        ofile << setw(15) << setprecision(8) << "  " << endl;
        for(int j = 0; j < N; j++){
           ofile << setw(15) << setprecision(8) << Eigenvectors(i,j) << endl;
      }
           ofile << setw(15) << setprecision(8) << "  " << endl; 
      }}


//Unit tests

    TEST_CASE("Testing maximal off-diagonal element"){
    int N=5;
    double rho_min = 0;
    double rho_max = 10.0;
    double h = (rho_max-rho_min)/N;
    cout << "Unit-test: ";
    mat a = zeros<mat>(N,N);
    mat v = zeros<mat>(N,N);
    vec rho(N);
    initialize(N,h,a,rho,v,0,0);
    int p=0;
    int q=0;
    double a_pq=0;

    //Find maximal element of matrix a
    Find_max_element(a,p,q,a_pq,N);
    
    REQUIRE(p==4);
    REQUIRE(q==3);
    REQUIRE(a_pq==Approx(-0.25));
} 


TEST_CASE("Testing agreement of eigenvalues with analytical solution"){
    int N = 5;
    double tolerance = pow(10.0, -8.0);
    double rho_min = 0;
    double rho_max = 1.0;
    double h = (rho_max-rho_min)/N;
    N = N-1;
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);
    
    mat a_matrix = zeros<mat>(N,N);
    mat v = zeros<mat>(N,N);
    vec rho(N);
    initialize(N,h,a_matrix,rho,v,0,0);
    //perform Jacobi algorithm
    Jacobi_algorithm(N,tolerance,a_matrix,v);
    
    //get eigenvalues-vector
    vector<double>Eigenvalues=Get_eigenvalues(a_matrix,N);
    
    //Require agreement with analytical solution
    double pi = acos(-1.0);
    REQUIRE(Eigenvalues[0]==Approx(d+2.0*a*cos(1*pi/(N+1))));
    REQUIRE(Eigenvalues[1]==Approx(d+2.0*a*cos(2*pi/(N+1))));
    REQUIRE(Eigenvalues[2]==Approx(d+2.0*a*cos(3*pi/(N+1))));
}

//Main program
int main(int argc, char* argv[]){
     //Result of unit-tests
     int test_result = Catch::Session().run( argc, argv );
     //Use Armadillo or Jacobi algorithm <----------------------------set "NO" to use Jacobi
     string use_arma = "NO";

     int N;
     double rho_min, rho_max, h;
     rho_min = 0.0;
     rho_max = 10; //set rho_max = 10. for case with two electrons
     N = 100;
     h = (rho_max-rho_min)/N;
     N = N-1; //u_0 and u_N are fixed by boundary conditions, proceed with points u_1 to u_N-1
     mat Hamiltonian_matrix = zeros<mat>(N,N);
     mat v =  zeros<mat>(N,N);
     vec rho(N);
     int num_electrons = 2; //0 - buckling beam, 1 - one electron, 2 - two electrons
     double omega_r0 = 0.0, omega_r1 = 0.01, omega_r2 = 0.5, omega_r3 = 1.0, omega_r4 = 5.0;
     double delta_convergence = pow(10.0, -8.0); //tolerance of convergence
     initialize(N, h, Hamiltonian_matrix, rho, v, num_electrons, omega_r4);
     if(use_arma != "NO"){
        Diagonalize_with_Armadillo(N, h, Hamiltonian_matrix, v, num_electrons);
     }
     else{
     Jacobi_algorithm(N, delta_convergence,Hamiltonian_matrix,v);
     Write_to_file(Hamiltonian_matrix, v, N, num_electrons,h);
     vec Eigenvalues = Get_eigenvalues(Hamiltonian_matrix, N);
     cout.precision(6);
     cout << Eigenvalues(0) <<"  " << Eigenvalues(1) << "  " << Eigenvalues(2) << "  " <<  Eigenvalues(3) << endl;
     }
   return test_result;
} 
