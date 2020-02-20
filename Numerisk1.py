# Remember to import nescessary libraries
import numpy as np
import matplotlib.pyplot as plt

# Physical constants
theta_0 = 0.2
g = 9.8
l = 1.0
phi = 0

# Timing constants
T = 10.0       # [s], Plot up to t=10s. 
dt = 0.0001    # [s], timestep
times = np.linspace(0,T,int(T/dt))

# Function definitions
theta_t = lambda t : theta_0*np.cos(np.sqrt(g/l)*t + phi)

#Script

theta_val = np.zeros(int(T/dt))

for i in range(0,len(times)):
    theta_val[i] = theta_t(times[i])

#Ploting
plt.figure(0)
plt.title("Analytisk løsning")
plt.plot(times,theta_val)
plt.xlabel("Time(s)")
plt.ylabel("Position(rad)")

# Physical constants
theta_0 = 0.2
omega_0 = 0
g = 9.8
m = 5
l = 1.0
phi = 0

# Timing constants
t_i = 0.0
T = T       # [s], Plot up to t=10s. 
dt = dt    # [s], timestep
times = np.linspace(0,T,int(T/dt))

# Function definitions

theta_t = lambda t: theta_0*np.cos(np.sqrt(g/l)*t + phi)


def euler_method(theta_0,omega_0,dt):
    theta=np.zeros(int(T/dt)+1)
    omega=np.zeros(int(T/dt)+1)
    t=np.zeros(int(T/dt)+1)
    theta[0]=theta_0
    omega[0]=omega_0
    t[0]=t_i

    for n in range(int(T/dt)):
        theta_new=theta[n]+omega[n]*dt
        omega_new=omega[n]-g/l*theta[n]*dt
        
        t[n+1]=t[n]+dt
        theta[n+1]=theta_new
        omega[n+1]=omega_new
    return theta,omega,t

theta, omega, tim = euler_method(theta_0,omega_0,dt)


plt.figure(1)
plt.title("Displacement")
plt.plot(times,theta_val,"--r",label="Analytical")
plt.plot(tim,theta,"--b",label="Numerical solution (Euler)")
plt.xlabel("Time(s)")
plt.ylabel("Position(rad)")
plt.legend(loc="upper right")


# 2.

def energy_calculation(theta_0, omega_0, dt):

    theta, omega, tim = euler_method(theta_0,omega_0,dt)

    energy_func = lambda m,l,omega,theta: (1/2)*m*(l**2)*(omega**2) + (1/2)*m*g*l*(theta**2)
    
    t = np.linspace(t_i,T,int(T/dt))
    energy = np.zeros(int(T/dt))
    
    for i in range(len(t)):
        energy[i] = energy_func(m,l,omega[i],theta[i])
    
    
    E_total = energy

    return t, E_total

# 3.

dt1 = 0.001
dt2 = 0.004
dt3 = 0.007

# Physical constants
theta_0 = 0.2
omega_0 = 0
g = 9.8
m = 5
l = 1.0
phi = 0

plt.figure(2)
plt.title("Total Energy")
plt.plot(energy_calculation(theta_0,omega_0,dt1)[0],energy_calculation(theta_0,omega_0,dt1)[1], label="dt = " + str(dt1))
plt.plot(energy_calculation(theta_0,omega_0,dt2)[0],energy_calculation(theta_0,omega_0,dt2)[1], label="dt = " + str(dt2))
plt.plot(energy_calculation(theta_0,omega_0,dt3)[0],energy_calculation(theta_0,omega_0,dt3)[1], label="dt = " + str(dt3))
plt.xlabel("Time(s)")
plt.ylabel("Energy(J)")
plt.legend(loc="upper right")
plt.show()