from numpy import *
import matplotlib.pyplot as plt


def read_and_plot(filename, n, chosen_color1, chosen_color2):
   infile = open(filename, 'r')
   x_values = zeros(n)
   exact_solution = zeros(n)
   approximation = zeros(n)
   i = 0
   for line in infile:
      values = line.split()
      x_values[i] = float(values[0])
      approximation[i] = float(values[1])
      exact_solution[i] = float(values[2])
      i += 1
   infile.close()
   plt.plot(x_values, exact_solution, color=chosen_color1)
   plt.plot(x_values, approximation, color=chosen_color2, linestyle='dashed')
   #plt.title('n=%.0f grid points' % (n), fontsize='15')
   plt.legend(['Exact solution','Approximation'], fontsize='17')
   plt.xlabel('x', fontsize='18')
   plt.savefig('n%.0f.png' % (n))
   plt.show()

# Create plots of exact and numerical solutions for different n values
read_and_plot('solution.out10', 10, 'teal', 'navy')
read_and_plot('solution.out100', 100, 'salmon', 'darkred')
read_and_plot('solution.out1000', 1000, 'darkgreen', 'lime')

#Create arrays (read off from the Terminal) and plot maximal relative error as a function of stepsize
max_relative_error = array([-1.1797,-3.08804,-5.08005,-7.07929,-8.84297, -6.07547, -5.52523])
log10h = array([-1.04139 ,-2.00432, -3.00043,-4.00004, -5., -6., -7.])
plt.plot(log10h, max_relative_error, color='tomato')
plt.scatter(log10h, max_relative_error, color='tomato')
plt.title('Relative error as a function of stepsize', fontsize='15')
plt.xlabel('log$_{10}$(h)', fontsize='15')
plt.ylabel('$\\varepsilon_{max}$',  fontsize='16')
plt.axis([-0.8,-7.1,-9.,-1.])
plt.savefig('error.png')
plt.show()
