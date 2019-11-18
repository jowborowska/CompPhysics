from numpy import *
import matplotlib.pyplot as plt
from scipy.interpolate import spline

# Read file determining the most likely state vs number of MC cycles, 4c)
def read_file_1(filename, n):
   infile = open(filename, 'r')
   MC_cycles = zeros(n)
   E_expectation = zeros(n)
   M_expectation = zeros(n)
   accepted = zeros(n)
   E_variance = zeros(n)
   E = zeros(n) #at current MC cycle
   i = 0
   for line in infile:
      values = line.split()
      MC_cycles[i] = float(values[0])
      E_expectation[i] = float(values[1])
      M_expectation[i] = float(values[2])
      accepted[i] = float(values[3])
      E_variance[i] = float(values[4])
      E[i] = float(values[5])
      i += 1
   infile.close()
   return MC_cycles, E_expectation, M_expectation, accepted, E_variance, E

# Read file for 4e)
def read_file_2(filename, n):
   infile = open(filename, 'r')
   T = zeros(n)
   E_expectation = zeros(n)
   C_V = zeros(n)
   M_expectation = zeros(n)
   chi = zeros(n)
   Mabs_expectation = zeros(n)
   i = 0
   for line in infile:
      if i > 0:
         values = line.split()
         T[i-1] = float(values[0])
         E_expectation[i-1] = float(values[1])
         C_V[i-1] = float(values[2])
         M_expectation[i-1] = float(values[3])
         chi[i-1] = float(values[4])
         Mabs_expectation[i-1] = float(values[5])
      i += 1
   infile.close()
   return T, E_expectation, C_V, M_expectation, chi, Mabs_expectation

#smooth the curves
def smooth(T, E, C_V, M, chi, Mabs):
   T_interpolated = linspace(2.1, 2.4, 100)
   E_interpolated = spline( T, E,T_interpolated)
   C_V_interpolated = spline( T, C_V,T_interpolated)
   M_interpolated = spline( T, M,T_interpolated)
   chi_interpolated = spline(T, chi,T_interpolated)
   Mabs_interpolated = spline(T, Mabs,T_interpolated)
   return T_interpolated, E_interpolated, C_V_interpolated, M_interpolated, chi_interpolated, Mabs_interpolated



T100, E_expectation100, C_V100, M_expectation100, chi100, Mabs_expectation100 = read_file_2("L100_final_dT002.tex", 16)
T80, E_expectation80, C_V80, M_expectation80, chi80, Mabs_expectation80 = read_file_2("L80_final_dT002.tex", 16)
T60, E_expectation60, C_V60, M_expectation60, chi60, Mabs_expectation60 = read_file_2("L60_final_dT002.tex", 16)
T40, E_expectation40, C_V40, M_expectation40, chi40, Mabs_expectation40 = read_file_2("L40_final_dT002.tex", 16)

'''
#these when we use values recorded after 10^4 MC cycles
T100, E_expectation100, C_V100, M_expectation100, chi100, Mabs_expectation100 = read_file_2("new_try_dT002_100.tex", 16)
T80, E_expectation80, C_V80, M_expectation80, chi80, Mabs_expectation80 = read_file_2("new_try_dT002_80.tex", 16)
T60, E_expectation60, C_V60, M_expectation60, chi60, Mabs_expectation60 = read_file_2("new_try_dT002_60.tex", 16)
T40, E_expectation40, C_V40, M_expectation40, chi40, Mabs_expectation40 = read_file_2("new_try_dT002_40.tex", 16)
'''


T100, E_expectation100, C_V100, M_expectation100, chi100, Mabs_expectation100 = smooth(T100, E_expectation100, C_V100, M_expectation100, chi100, Mabs_expectation100)   
T80, E_expectation80, C_V80, M_expectation80, chi80, Mabs_expectation80 = smooth(T80, E_expectation80, C_V80, M_expectation80, chi80, Mabs_expectation80)   
T60, E_expectation60, C_V60, M_expectation60, chi60, Mabs_expectation60 = smooth(T60, E_expectation60, C_V60, M_expectation60, chi60, Mabs_expectation60)   
T40, E_expectation40, C_V40, M_expectation40, chi40, Mabs_expectation40 = smooth(T40, E_expectation40, C_V40, M_expectation40, chi40, Mabs_expectation40)   

#Find critical temperatures for different L at the chi plot
T_C_40 = T40[argmax(chi40)]
T_C_60 = T60[argmax(chi60)]
T_C_80 = T80[argmax(chi80)]
T_C_100 = T100[argmax(chi100)]
print T_C_40, T_C_60, T_C_80, T_C_100 #not-smoothed: 2.32 2.3 2.3 2.28, smoothed: 2.315151515151515 2.306060606060606 2.309090909090909 2.2818181818181817


plt.figure()
plt.plot(T40, E_expectation40)
plt.plot(T60, E_expectation60)
plt.plot(T80, E_expectation80)
plt.plot(T100, E_expectation100)
plt.legend(['L=40','L=60','L=80','L=100'], fontsize='15')
plt.xlabel('kT/J', fontsize='15')
plt.ylabel('<E>', fontsize='14')
plt.show()   

plt.figure()
plt.plot(T40, Mabs_expectation40)
plt.plot(T60, Mabs_expectation60)
plt.plot(T80, Mabs_expectation80)
plt.plot(T100, Mabs_expectation100)
plt.legend(['L=40','L=60','L=80','L=100'], fontsize='15')
plt.xlabel('kT/J', fontsize='15')
plt.ylabel('<|M|>', fontsize='14')
plt.show()   

plt.figure()
plt.plot(T40, chi40)
plt.scatter(T_C_40, chi40[argmax(chi40)])
plt.plot(T60, chi60)
plt.scatter(T_C_60, chi60[argmax(chi60)])
plt.plot(T80, chi80)
plt.scatter(T_C_80, chi80[argmax(chi80)])
plt.plot(T100, chi100)
plt.scatter(T_C_100, chi100[argmax(chi100)])
plt.legend(['L=40','L=60','L=80','L=100'], fontsize='15')
plt.xlabel('kT/J', fontsize='15')
plt.ylabel('$\chi$', fontsize='14')
plt.show()   

#Find critical temperatures for different L at the C_V plot
T_C_40 = T40[argmax(C_V40)]
T_C_60 = T60[argmax(C_V60)]
T_C_80 = T80[argmax(C_V80)]
T_C_100 = T100[argmax(C_V100)]
print T_C_40, T_C_60, T_C_80, T_C_100 #smoothed: 2.2818181818181817 2.293939393939394 2.278787878787879 2.275757575757576


plt.figure()
plt.plot(T40, C_V40)
plt.scatter(T_C_40, C_V40[argmax(C_V40)])
plt.plot(T60, C_V60)
plt.scatter(T_C_60, C_V60[argmax(C_V60)])
plt.plot(T80, C_V80)
plt.scatter(T_C_80, C_V80[argmax(C_V80)])
plt.plot(T100, C_V100)
plt.scatter(T_C_100, C_V100[argmax(C_V100)])
plt.legend(['L=40','L=60','L=80','L=100'], fontsize='15')
plt.xlabel('kT/J', fontsize='15')
plt.ylabel('$C_V$', fontsize='14')
plt.show()   


MC, E_expect1ord, M_expect1ord, acc1ord, E_variance1ord, E1ord = read_file_1('likely_T1_ordered_new.tex', int(1e6))
MC, E_expect1nord, M_expect1nord, acc1nord, E_variance1nord, E1nord = read_file_1('likely_T1_notordered_new.tex', int(1e6))
MC, E_expect2ord, M_expect2ord, acc2ord, E_variance2ord, E2ord = read_file_1('likely_T2_ordered_new.tex', int(1e6))
MC, E_expect2nord, M_expect2nord, acc2nord, E_variance2nord, E2nord = read_file_1('likely_T2_notordered_new.tex', int(1e6))


plt.figure()
plt.plot(MC, M_expect1ord, color='darkcyan')
plt.plot(MC, M_expect1nord, color='navy')
plt.xscale("log")
plt.xlabel('MC cycles', fontsize='15')
plt.ylabel('<|M|>', fontsize='15')
plt.legend(['Ordered', 'Random'],fontsize='15' )
plt.show()   

plt.figure()
plt.plot(MC, E_expect1ord, color='darkcyan')
plt.plot(MC, E_expect1nord, color='navy')
plt.xscale("log")
plt.xlabel('MC cycles', fontsize='15')
plt.ylabel('<E>', fontsize='15')
plt.legend(['Ordered', 'Random'],fontsize='15' )
plt.show()   

plt.figure()
plt.plot(MC, M_expect2ord, color='darkcyan')
plt.plot(MC, M_expect2nord, color='navy')
plt.xscale("log")
plt.xlabel('MC cycles', fontsize='15')
plt.ylabel('<|M|>', fontsize='15')
plt.legend(['Ordered', 'Random'],fontsize='15' )
plt.show()   

plt.figure()
plt.plot(MC, E_expect2ord, color='darkcyan')
plt.plot(MC, E_expect2nord, color='navy')
plt.xscale("log")
plt.xlabel('MC cycles', fontsize='15')
plt.ylabel('<E>', fontsize='15')
plt.legend(['Ordered', 'Random'],fontsize='15' )
plt.show()   

plt.figure()
plt.plot(MC, acc1ord, color='darkcyan')
plt.plot(MC, acc1nord, color='navy')
plt.xscale("log")
plt.xlabel('MC cycles', fontsize='15')
plt.ylabel('no. of accepted conf.', fontsize='15')
plt.legend(['Ordered', 'Random'],fontsize='15' )
plt.show()   

plt.figure()
plt.plot(MC, acc2ord, color='darkcyan')
plt.plot(MC, acc2nord, color='navy')
plt.xscale("log")
plt.xlabel('MC cycles', fontsize='15')
plt.ylabel('no. of accepted conf.', fontsize='15')
plt.legend(['Ordered', 'Random'],fontsize='15' )
plt.show()   



E_probability_T1 = E1ord[10000:1000000]
print E_variance1ord[-1] #9.2953742
E1_variance = (sum(E_probability_T1**2))/(1e6-1e4) - (sum(E_probability_T1)/(1e6-1e4))**2
print E1_variance #9.29840234562289

E_probability_T2 = E2ord[10000:1000000]
print E_variance2ord[-1] #3250.0166
E2_variance = (sum(E_probability_T2**2))/(1e6-1e4) - (sum(E_probability_T2)/(1e6-1e4))**2
print E2_variance #3246.3098298036784


plt.hist(E_probability_T1, bins=30, density=True, color='firebrick')
plt.xlabel('E', fontsize='15')
plt.ylabel('P(E)', fontsize='15')
plt.show()

plt.hist(E_probability_T2, bins=30, density=True, color='salmon')
plt.xlabel('E', fontsize='15')
plt.ylabel('P(E)', fontsize='15')
plt.show()



