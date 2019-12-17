from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def read_energy_momentum(filename, time_step, total_time):
   infile = open(filename, 'r')
   Energy_list = []
   Momentum_list = []
   kinetic_list = []
   potential_list = []
   for line in infile:
      values = line.split()
      E = float(values[0])
      L = float(values[1])
      Ek =  float(values[2])
      Ep =  float(values[3])
      Energy_list.append(E)
      Momentum_list.append(L)
      kinetic_list.append(Ek)
      potential_list.append(Ep)
   Energy_array = array(Energy_list)
   Momentum_array = array(Momentum_list)
   kinetic_array = array(kinetic_list)
   potential_array = array(potential_list)
   infile.close()
   #create time array
   infile = open(filename, 'r')
   number_of_lines = len(infile.readlines())
   time = linspace(time_step, total_time, number_of_lines)
   infile.close()
   return Energy_array, Momentum_array, time, kinetic_array, potential_array

def read_positions(filename, bodies_number):
   infile = open(filename, 'r')
   number_of_lines = len(infile.readlines())
   position_matrix = zeros((number_of_lines/bodies_number, bodies_number, 3))
   j = 0
   infile.close()
   infile = open(filename, 'r')
   for line in infile:
      values = line.split()
      body = int(values[0])
      x = float(values[1])
      y = float(values[2])
      z = float(values[3])
      position_matrix[j,body-1,:] = array([x,y,z])
      if body == bodies_number:
         j += 1
   infile.close()
   return position_matrix

def read_positions_precession(filename, bodies_number):
   infile = open(filename, 'r')
   number_of_lines = len(infile.readlines())
   position_matrix = zeros((number_of_lines/bodies_number, bodies_number, 3))
   vec_length_matrix = zeros((number_of_lines/bodies_number, bodies_number))
   j = 0
   infile.close()
   infile = open(filename, 'r')
   for line in infile:
      values = line.split()
      body = int(values[0])
      x = float(values[1])
      y = float(values[2])
      z = float(values[3])
      vector_length = float(values[4]) #added for Mercury precession
      position_matrix[j,body-1,:] = array([x,y,z])
      vec_length_matrix[j,body-1] = array([vector_length])
      if body == bodies_number:
         j += 1
   infile.close()
   return position_matrix, vec_length_matrix

#create the function to find latest index corresponding to perihelion (argmin would find the earliest one)
def find_min_distance(distance):
   minimum = 0.3075
   for i in range(len(distance)):
      if distance[i] == minimum:
         index_min = i
   return index_min

positions_Mercury_nocorrection, length_no = read_positions_precession("positions_Mercury_nocorrection.xyz",2) #100years, 1e8 timesteps
positions_Mercury_correction, length_with = read_positions_precession("positions_Mercury_correction.xyz",2) #100years, 1e8 timesteps, with relativistic effect corrected

Mercury_position_no = positions_Mercury_nocorrection[:,1,:]
Mercury_position_with = positions_Mercury_correction[:,1,:]

distance_no = length_no[:,1]
distance_with =length_with[:,1]

index_min_no = find_min_distance(distance_no)
index_min_with = find_min_distance(distance_with)

tan_perihelion_no = Mercury_position_no[index_min_no,1]/Mercury_position_no[index_min_no,0]
tan_perihelion_with = Mercury_position_with[index_min_with,1]/Mercury_position_with[index_min_with,0]

print "Net angle: ",  (arctan(tan_perihelion_with)*360/(2.*pi)-arctan(tan_perihelion_no)*360/(2.*pi))*60*60 #39.591952009684064 arcsec


plt.figure()

plt.plot(Mercury_position_no[:,0], Mercury_position_no[:,1], color='darkcyan', label='Mercury, no correction')
plt.plot(Mercury_position_with[:,0], Mercury_position_with[:,1], color='maroon', label='Mercury, with correction')
plt.scatter(Mercury_position_with[index_min_no,0],Mercury_position_with[index_min_no,1], label='perihelion final' )
plt.scatter(0.3075, 0, label='initial perihelion')
plt.legend(fontsize=14)
plt.xlabel('x [AU]', fontsize=14)
plt.ylabel('y [AU]', fontsize=14)
plt.show()



#whole model, 150years, 1.5e6 time-steps
positions_model = read_positions("positions_model.xyz",10)

Sun_position_model = positions_model[:,0,:]
Mercury_position_model = positions_model[:,1,:]
Venus_position_model = positions_model[:,2,:]
Earth_position_model = positions_model[:,3,:]
Mars_position_model = positions_model[:,4,:]
Jupiter_position_model = positions_model[:,5,:]
Saturn_position_model = positions_model[:,6,:]
Uranus_position_model = positions_model[:,7,:]
Neptune_position_model = positions_model[:,8,:]
Pluto_position_model = positions_model[:,9,:]

plt.figure()
plt.scatter(Sun_position_model[:,0], Sun_position_model[:,1], s=20, color='gold', label='Sun')
plt.plot(Mercury_position_model[:,0], Mercury_position_model[:,1], color='maroon', label='Mercury')
plt.plot(Venus_position_model[:,0], Venus_position_model[:,1], color='hotpink', label='Venus')
plt.plot(Earth_position_model[:,0], Earth_position_model[:,1], color='lightseagreen', label='Earth')
plt.plot(Mars_position_model[:,0], Mars_position_model[:,1], color='salmon', label='Mars')
plt.plot(Jupiter_position_model[:,0], Jupiter_position_model[:,1], color='orchid', label='Jupiter')
plt.plot(Saturn_position_model[:,0], Saturn_position_model[:,1], color='peru', label='Saturn')
plt.plot(Uranus_position_model[:,0], Uranus_position_model[:,1], color='seagreen', label='Uranus')
plt.plot(Neptune_position_model[:,0], Neptune_position_model[:,1], color='cornflowerblue', label='Neptune')
plt.plot(Pluto_position_model[:,0], Pluto_position_model[:,1], color='purple', label='Pluto')
plt.legend(fontsize=16, loc='right', bbox_to_anchor=(1.1, 0.5))
plt.xlabel('x [AU]', fontsize=15)
plt.ylabel('y [AU]', fontsize=15)
plt.show()

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(Sun_position_model[::100,0], Sun_position_model[::100,1], Sun_position_model[::100,2], color='gold', label='Sun')
ax.plot3D(Mercury_position_model[::100,0], Mercury_position_model[::100,1], Mercury_position_model[::100,2], color='maroon', label='Mercury')
ax.plot3D(Venus_position_model[::100,0], Venus_position_model[::100,1],Venus_position_model[::100,2], color='hotpink', label='Venus')
ax.plot3D(Earth_position_model[::100,0], Earth_position_model[::100,1], Earth_position_model[::100,2], color='lightseagreen', label='Earth')
ax.plot3D(Mars_position_model[::100,0], Mars_position_model[::100,1],Mars_position_model[::100,2], color='salmon', label='Mars')
ax.plot3D(Jupiter_position_model[::100,0], Jupiter_position_model[::100,1],Jupiter_position_model[::100,2], color='orchid', label='Jupiter')
ax.plot3D(Saturn_position_model[::100,0], Saturn_position_model[::100,1],Saturn_position_model[::100,2], color='peru', label='Saturn')
ax.plot3D(Uranus_position_model[::100,0], Uranus_position_model[::100,1], Uranus_position_model[::100,2], color='seagreen', label='Uranus')
ax.plot3D(Neptune_position_model[::100,0], Neptune_position_model[::100,1],  Neptune_position_model[::100,2],color='cornflowerblue', label='Neptune')
ax.plot3D(Pluto_position_model[::100,0], Pluto_position_model[::100,1],  Pluto_position_model[::100,2],color='purple', label='Pluto')
#ax.legend()
ax.set_xlabel('x [AU]')
ax.set_ylabel('y [AU]')
ax.set_zlabel('z [AU]')
plt.show()


#15years, 1.5e6 datapoins, Sun, Earth, Jupiter move about the center of the mass
positions_Jupiter_Sun_moves = read_positions("positions_Jupiter_Sun_moves_report.xyz", 3)

Sun_position_Jupiter_Sun_moves = positions_Jupiter_Sun_moves[:,0,:]
Earth_position_Jupiter_Sun_moves = positions_Jupiter_Sun_moves[:,1,:]
Jupiter_position_Jupiter_Sun_moves = positions_Jupiter_Sun_moves[:,2,:]

plt.figure()
plt.scatter(Sun_position_Jupiter_Sun_moves[:,0], Sun_position_Jupiter_Sun_moves[:,1], s=60, color='gold', label='Sun')
plt.plot(Earth_position_Jupiter_Sun_moves[:,0], Earth_position_Jupiter_Sun_moves[:,1], color='lightseagreen', label='Earth')
plt.plot(Jupiter_position_Jupiter_Sun_moves[:,0], Jupiter_position_Jupiter_Sun_moves[:,1], color='salmon', label='Jupiter')
plt.legend(fontsize=14)
plt.xlabel('x [AU]', fontsize=14)
plt.ylabel('y [AU]', fontsize=14)
plt.show()


#15years, 1.5e6 datapoins, Sun kept at rest as the origin (so we move Earth and Jupiter coordinates from NASA with respect, and keep 0 position, 0 velocity for Sun)
positions_Jupiter = read_positions("positions_Jupiter_report.xyz", 3) #Jupiter with its regular mass
positions_Jupiter_10mass = read_positions("positions_Jupiter_10mass_report.xyz", 3) #Jupiter with 10 times its mass
positions_Jupiter_100mass = read_positions("positions_Jupiter_100mass_report.xyz", 3) #Jupiter with 100 times its mass
positions_Jupiter_1000mass = read_positions("positions_Jupiter_1000mass_report.xyz", 3) #Jupiter with 1000 times its mass

Sun_position_Jupiter = positions_Jupiter[:,0,:]
Earth_position_Jupiter = positions_Jupiter[:,1,:]
Jupiter_position_Jupiter = positions_Jupiter[:,2,:]

Sun_position_Jupiter_10mass = positions_Jupiter_10mass[:,0,:]
Earth_position_Jupiter_10mass = positions_Jupiter_10mass[:,1,:]
Jupiter_position_Jupiter_10mass = positions_Jupiter_10mass[:,2,:]

Sun_position_Jupiter_100mass = positions_Jupiter_100mass[:,0,:]
Earth_position_Jupiter_100mass = positions_Jupiter_100mass[:,1,:]
Jupiter_position_Jupiter_100mass = positions_Jupiter_100mass[:,2,:]

Sun_position_Jupiter_1000mass = positions_Jupiter_1000mass[:,0,:]
Earth_position_Jupiter_1000mass = positions_Jupiter_1000mass[:,1,:]
Jupiter_position_Jupiter_1000mass = positions_Jupiter_1000mass[:,2,:]

#plt.figure()
plt.subplots(2,2,sharex=True, sharey=True)
plt.subplot(2,2,1)
plt.scatter(Sun_position_Jupiter[0,0], Sun_position_Jupiter[0,1], s=60, color='gold', label='Sun')
plt.plot(Earth_position_Jupiter[:,0], Earth_position_Jupiter[:,1], color='lightseagreen', label='Earth')
plt.plot(Jupiter_position_Jupiter[:,0], Jupiter_position_Jupiter[:,1], color='salmon', label='Jupiter')
plt.legend(fontsize=15, loc='upper right')
plt.xlabel('x [AU]', fontsize=15)
plt.ylabel('y [AU]', fontsize=15)
plt.subplot(2,2,2)
plt.scatter(Sun_position_Jupiter_10mass[0,0], Sun_position_Jupiter_10mass[0,1], s=60, color='gold')
plt.plot(Earth_position_Jupiter_10mass[:,0], Earth_position_Jupiter_10mass[:,1], color='lightseagreen')
plt.plot(Jupiter_position_Jupiter_10mass[:,0], Jupiter_position_Jupiter_10mass[:,1], color='salmon', label='10 times heavier Jupiter')
plt.legend(fontsize=15, loc='upper right')
plt.xlabel('x [AU]', fontsize=15)
plt.ylabel('y [AU]', fontsize=15)
plt.subplot(2,2,3)
plt.scatter(Sun_position_Jupiter_100mass[0,0], Sun_position_Jupiter_100mass[0,1], s=60, color='gold')
plt.plot(Earth_position_Jupiter_100mass[:,0], Earth_position_Jupiter_100mass[:,1], color='lightseagreen')
plt.plot(Jupiter_position_Jupiter_100mass[:,0], Jupiter_position_Jupiter_100mass[:,1], color='salmon', label='100 times heavier Jupiter')
plt.legend(fontsize=15, loc='upper right')
plt.xlabel('x [AU]', fontsize=15)
plt.ylabel('y [AU]', fontsize=15)
plt.subplot(2,2,4)
plt.scatter(Sun_position_Jupiter_1000mass[0,0], Sun_position_Jupiter_1000mass[0,1], s=60, color='gold')
plt.plot(Earth_position_Jupiter_1000mass[:,0], Earth_position_Jupiter_1000mass[:,1], color='lightseagreen')
plt.plot(Jupiter_position_Jupiter_1000mass[:,0], Jupiter_position_Jupiter_1000mass[:,1], color='salmon', label='1000 times heavier Jupiter')
plt.legend(fontsize=15, loc='upper right')
plt.xlabel('x [AU]', fontsize=15)
plt.ylabel('y [AU]', fontsize=15)
plt.show()


#5years, 5e5 data-points, Sun at rest, Earth: v0=7, x0 = 1AU 
positions_beta_2 = read_positions("positions_beta_2_report.xyz",2) #beta = 2
positions_beta_22 = read_positions("positions_beta_22_report.xyz",2) #beta = 2.2
positions_beta_25 = read_positions("positions_beta_25_report.xyz",2) #beta = 2.5
positions_beta_28 = read_positions("positions_beta_28_report.xyz",2) #beta = 2.8
positions_beta_29 = read_positions("positions_beta_29_report.xyz",2) #beta = 2.9
positions_beta_3 = read_positions("positions_beta_3_report.xyz",2) #beta = 3

Sun_position = positions_beta_2[:,0,:]
Earth_position_beta_2 = positions_beta_2[:,1,:]
Earth_position_beta_22 = positions_beta_22[:,1,:]
Earth_position_beta_25 = positions_beta_25[:,1,:]
Earth_position_beta_28 = positions_beta_28[:,1,:]
Earth_position_beta_29 = positions_beta_29[:,1,:]
Earth_position_beta_3 = positions_beta_3[:,1,:]

plt.figure()
plt.scatter(Sun_position[:,0], Sun_position[:,1], s=60, color='gold', label='Sun')
plt.plot(Earth_position_beta_2[:,0], Earth_position_beta_2[:,1], color='lightseagreen', label='Earth, $\\beta = 2$')
plt.plot(Earth_position_beta_22[:,0], Earth_position_beta_22[:,1], color='hotpink', label='Earth, $\\beta = 2.2$')
plt.plot(Earth_position_beta_25[:,0], Earth_position_beta_25[:,1], color='purple', label='Earth, $\\beta = 2.5$')
plt.plot(Earth_position_beta_28[:,0], Earth_position_beta_28[:,1], color='navy', label='Earth, $\\beta = 2.8$')
plt.plot(Earth_position_beta_29[:,0], Earth_position_beta_29[:,1], color='orchid', label='Earth, $\\beta = 2.9$')
plt.plot(Earth_position_beta_3[:,0], Earth_position_beta_3[:,1], color='green', label='Earth, $\\beta = 3$')
plt.legend(fontsize=14, loc='upper left',  bbox_to_anchor=(0.0, 1.1))
plt.xlabel('x [AU]', fontsize=14)
plt.ylabel('y [AU]', fontsize=14)
plt.show()

#these are produced over 100 years, 1e6 data-points, Verlet  velocity algorithm, Sun at rest 
positions_escape = read_positions("positions_escape_vel_report.xyz",2) #v = escape velocity
positions_vel_8 = read_positions("positions_vel_8_report.xyz",2) #v = 8
positions_vel_88 = read_positions("positions_vel_88_report.xyz",2) #v = 8.8
positions_vel_89 = read_positions("positions_vel_89_report.xyz",2) #v = 8.9
positions_vel_9 = read_positions("positions_vel_9_report.xyz",2) #v = 9
positions_vel_7 = read_positions("positions_vel_7_report.xyz",2) #v = 7

Sun_position = positions_escape[:,0,:] 
Earth_position_v_escape = positions_escape[:,1,:]
Earth_position_v_8 = positions_vel_8[:,1,:]
Earth_position_v_88 = positions_vel_88[:,1,:]
Earth_position_v_89 = positions_vel_89[:,1,:]
Earth_position_v_9 = positions_vel_9[:,1,:]
Earth_position_v_7 = positions_vel_7[:,1,:]

plt.figure()
plt.scatter(Sun_position[:,0], Sun_position[:,1], s=20, color='gold', label='Sun')
plt.plot(Earth_position_v_7[int(9e5):int(1e6-1),0], Earth_position_v_7[int(9e5):int(1e6-1),1], color='green', label='Earth, $v_0 = 7$')
plt.plot(Earth_position_v_8[int(9e5):int(1e6-1),0], Earth_position_v_8[int(9e5):int(1e6-1),1], color='purple', label='Earth, $v_0 = 8$')
plt.plot(Earth_position_v_88[:,0], Earth_position_v_88[:,1], color='steelblue', label='Earth, $v_0 = 8.8$')
plt.plot(Earth_position_v_89[:,0], Earth_position_v_89[:,1], color='cyan', label='Earth, $v_0 = 8.9$')
plt.plot(Earth_position_v_escape[:,0], Earth_position_v_escape[:,1], color='firebrick', label='Earth, $v_0 = 2 \pi \sqrt{2} \\approx 8.886 $ ')
plt.plot(Earth_position_v_9[:,0], Earth_position_v_9[:,1], color='hotpink', label='Earth, $v_0 = 9$')
plt.legend(fontsize=14, loc='lower right', bbox_to_anchor=(1.1, 0.5))
plt.xlabel('x [AU]', fontsize=14)
plt.ylabel('y [AU]', fontsize=14)
plt.show()


positions = read_positions("positions_report.xyz",2) #v = 2pi
Sun_position = positions[:,0,:] 
Earth_position = positions[:,1,:]

def distance(x,y):
   return sqrt((1.-x)**2 + y**2)

print Earth_position[-1] #read last position
#distamces between initial and final position
N_array = array([1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7])
Verlet_array = array([8.26822e-5, 3.30733e-6, 8.26834e-7, 3.30736e-8, 8.26821e-9, 3.31384e-10, 8.2052e-11, 4.10248e-12, 5.5352e-13])
Euler_array = array([distance(1.01557,0.358194), distance(1.01299, 0.0740556), distance(1.00719,0.0371252), distance(1.00155,0.00743846), distance(1.00078, 0.00372), distance(1.00016, 7.44121e-4), distance(1.00008,3.72068e-4), distance(1.00002, 7.44148e-5), distance(1.00001, 3.72075e-5)])

plt.figure()
plt.plot(N_array, Verlet_array,color='steelblue', label='velocity Verlet')
plt.scatter(N_array, Verlet_array,color='steelblue')
plt.plot(N_array, Euler_array, color='hotpink', label='Euler')
plt.scatter(N_array, Euler_array, color='hotpink')
plt.legend(fontsize=15)
plt.xlabel('N', fontsize=15)
plt.ylabel('$\Delta r \: [\mathrm{AU}]$', fontsize=15)
plt.xscale("log")
plt.yscale("log")
plt.show()

plt.figure()
plt.scatter(Sun_position[:,0], Sun_position[:,1], s=60, color='gold', label='Sun')
plt.plot(Earth_position[:,0], Earth_position[:,1], color='lightseagreen', label='Earth')
plt.legend(fontsize=14)
plt.xlabel('x [AU]', fontsize=14)
plt.ylabel('y [AU]', fontsize=14)
plt.show()

#keep Sun at rest, Earth starts at 1AU, with speed 2pi
E_Verlet, L_Verlet, time_Verlet, kinetic_Verlet, potential_Verlet = read_energy_momentum("energy_momentum_Verlet_report.txt", 1e-5, 10.)
E_Euler, L_Euler, time_Euler, kinetic_Euler, potential_Euler = read_energy_momentum("energy_momentum_Euler_report.txt", 1e-5, 10.)

plt.figure()
plt.plot(time_Verlet, E_Verlet, color='salmon', label='$\Delta E_{Verlet}$')
plt.plot(time_Verlet, L_Verlet, color='darkcyan', label='$\Delta L_{Verlet}$')
plt.plot(time_Euler[::10], L_Euler[::10], color='navy', label='$\Delta L_{Euler}$')
plt.plot(time_Euler[::10], E_Euler[::10], color='maroon', label='$\Delta E_{Euler}$')
plt.yscale("log")
plt.legend(fontsize=14)
plt.xlabel('t [yr]', fontsize=14)
plt.show()

plt.figure()
plt.subplot(2,1,1)
plt.plot(time_Verlet, kinetic_Verlet, color='darkcyan', label='$E_{k,Verlet}$',)
plt.plot(time_Euler[::10], kinetic_Euler[::10], color='maroon', label='$E_{k,Euler}$')
plt.legend(fontsize=15)
#plt.xlabel('t [yr]', fontsize=14)
plt.ylabel('$[\mathrm{AU}^2 M_{\odot}/\mathrm{yr}^2]$', fontsize=15)
plt.subplot(2,1,2)
plt.plot(time_Verlet, potential_Verlet, color='navy', label='$E_{p,Verlet}$')
plt.plot(time_Euler, potential_Euler, color='firebrick', label='$E_{p,Euler}$')
plt.legend(fontsize=15)
plt.ylabel('$[\mathrm{AU}^2 M_{\odot}/\mathrm{yr}^2]$', fontsize=15)
plt.xlabel('t [yr]', fontsize=15)
plt.show()

