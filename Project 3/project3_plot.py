#Plots the single-particle hydrogenic wave function

from numpy import *
import matplotlib.pyplot as plt

def single_wave_function(limit, N):
   alpha = 2.
   domain = linspace(0., limit, N)
   return domain, exp(-alpha*domain)

r_array, wave_array = single_wave_function(4.,1000)

plt.figure()
plt.plot(r_array, wave_array, color='salmon')
plt.xlabel(r'$r_i$', fontsize = 16)
plt.ylabel(r'$\psi_{1s}$', fontsize = 16)
plt.show()
