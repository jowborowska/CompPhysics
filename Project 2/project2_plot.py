from numpy import *
import matplotlib.pyplot as plt

#Convergence and efficiency of Jacobi's method - the buckling beam problem
no_transformations = array([122,610,1457,4161,10928,17323,25007,39353,57024,70229,85608,110360])
dim_N_Jacobi = array([10,20,30,50,80,100,120,150,180,200,220,250])
CPU_time_Jacobi = array([0.001756,0.038143,0.138336, 0.468779,1.19434,2.42911,4.80059,11.0612,23.116,35.8579,53.8603,88.4976])
dim_N_Armadillo = array([10,20,30,50,80,100,120,150,180,200,220,250,400,500,600,800,1000,1500,2000,2200,2500])
CPU_time_Armadillo = array([0.000812, 0.000904,0.004055,0.016012,0.019633,0.021261,0.038311, 0.096616,  0.055834,0.060839, 0.061215,0.061215,0.174384,0.284758,0.407042,0.789388,1.56462,5.20555,12.6725,14.5773,21.184])

#One electron in 3D harmonic oscillator potential
lambda_analytical = ['3  7  11  15']
q_max = [5,5,7,7,7,9,9,10,10,10]
N = [100,200,100,200,300,200,300,200,300,400]
lambda_numerical = ['2.99922  6.99609  10.9907  14.9881','2.9998  6.99903  10.9978  15.0014','2.99847  6.99234  10.9813  14.9653','2.99962  6.99809  10.9953  14.9913','2.99983  6.99915  10.9979  14.9962','2.99937  6.99683  10.9923  14.9857','2.99972  6.99859  10.9966  14.9936','2.99922  6.99609  10.9905  14.9823', '2.99965  6.99826  10.9958  14.9921','2.9998  6.99902  10.9976  14.9956']

#Two interacting electrons in 3D harmonic oscillator potentail
omega_r = array([0.01,0.5,1.,5.])
lambda_0 = array([0.31162 ,2.22932 ,4.0546,17.3664]) 

plt.figure()
plt.title('Convergence of Jacobi algorithm', fontsize=14)
plt.plot(dim_N_Jacobi, 1.8*dim_N_Jacobi**2, color='darkcyan')
plt.scatter(dim_N_Jacobi, no_transformations,color='salmon')
plt.yscale('log')
plt.xlabel('Dimension of matrix, N', fontsize=14)
plt.ylabel('Number of Jacobi rotations',fontsize=14)
plt.legend(['Approximation,$1.8N^2$','Exact'], fontsize=14)
plt.savefig('convergence_plot.png')
plt.show()

plt.figure()
plt.title('Jacobi algorithm', fontsize='17')
plt.plot(dim_N_Jacobi, CPU_time_Jacobi, color='purple')
plt.scatter(dim_N_Jacobi, CPU_time_Jacobi, color='purple')
plt.xlabel('N',fontsize=17)
plt.ylabel('CPU time [s]', fontsize=17)
plt.savefig('efficiency_plot_jacobi.png')
plt.show()

plt.figure()
plt.title(' $\mathsf{eig\_sym}$ from Armadillo', fontsize='17')
plt.plot(dim_N_Armadillo, CPU_time_Armadillo, color='violet')
plt.scatter(dim_N_Armadillo, CPU_time_Armadillo, color='violet')
plt.xlabel('N',fontsize=17)
plt.ylabel('CPU time [s]', fontsize=17)
plt.savefig('efficiency_plot_arma.png')
plt.show()

plt.figure()
plt.scatter(omega_r,lambda_0, color = 'salmon')
plt.plot(omega_r,lambda_0, color='salmon')
plt.xlabel('$\omega_r$', fontsize=15)
plt.ylabel('$\lambda_0$',fontsize=15)
plt.savefig('omega_plot.png')
plt.show()
