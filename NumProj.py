import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Physical parameters

l = 1.0
g = 9.8
theta_0 = 0.2
omega_0 = 0.0
q = 1.0
omega_D = 3.13
F_D = 0.2

# Timing constants

t_f = 20
t_i = 0
dt = 0.01

# Functions

#       I think these are right, might want to do some checking though

d_theta = lambda theta, omega, t : omega

d_omega = lambda theta, omega, t : F_D*np.sin(omega_D*t)-((g/l)*theta)-(q*omega)



def RK4_step(f, k, theta, w, dt, t):
    """
    Calculates one step of the RK4-algorithm.
    
    theta: float
    previous value of theta
           
    w: float
    previous value of w (omega, angular velocity)
    
    dt: float
    timestep
    
    return: two floats 
    """
    k1 = k(theta,w,t)
    f1 = f(theta,w,t)
    k2 = k(theta + (dt/2)*f1,w,t)
    f2 = f(theta,w + (dt/2)*k1,t)
    k3 = k(theta + (dt/2)*f2,w,t)
    f3 = f(theta,w + (dt/2)*k2,t)
    k4 = k(theta + dt*f3,w,t)
    f4 = f(theta,w + dt*k3,t)
    return theta + (dt/6)*(f1 + (2*f2) + (2*f3) + f4), w + (dt/6)*(k1 + (2*k2) + (2*k3) + k4)

def RK4_method(k, f, theta__0, omega__0, dt):
    """
    Computes theta and w (omega).  
    
    Parameters
    -----------
    k: RHS of equation
    f: RHS of equation
    theta0: initial value of theta
    w0: initail value of omega
    dt: timestep
    
    return theta, w, t
    """
    t = np.linspace(0,t_f,int(t_f/dt))
    theta = np.zeros(len(t))
    omega = np.zeros(len(t))
    
    
    theta[0],omega[0] = theta__0,omega__0
    
    for i in range(1,len(t)):
        theta[i],omega[i] = RK4_step(k, f, theta[i-1], omega[i-1], dt, t[i])
    
    return theta, omega, t

""" Used this to figure out what the RK4-function did wrong and suceeded!
def euler_cromer_approx(theta_0, w_0, dt):
    
    #Calculates angular displacement and angular velocity 
    #using the Euler-Cromer method 
    
    N = int(t_f/dt)
    theta = np.zeros(N)
    w = np.zeros(N)
    t = np.linspace(0, t_f, N)
    theta[0] = theta_0
    w[0] = w_0
    for i in range(1,N):
        w[i] = w[i-1] + (F_D*np.sin(omega_D*t[i])-((g/l)*theta[i-1])-(q*w[i-1]))*dt
        theta[i] = theta[i-1] + w[i]*dt
    return theta, w, t
"""



theta_RK4,omega_RK4,t_RK4 = RK4_method(d_theta, d_omega, theta_0, omega_0, dt)
#theta_ec,omega_ec,t_ec = euler_cromer_approx(theta_0, omega_0, dt)

"""
plt.figure("Euler-Cromer")
plt.title("Calculation using Euler-Cromer")
plt.plot(t_ec,theta_ec,label="Displacement (rad)")
plt.legend(loc="upper right")
plt.show()
"""


plt.figure("RK4")
plt.title("Calculation using Runge-Kutta 4")
plt.plot(t_RK4,theta_RK4,label="Displacement (rad)")
plt.legend(loc="upper right")
plt.show()

"""
plt.figure("diff")
plt.title("Difference in Calculation")
plt.plot(t_RK4,theta_ec - theta_RK4,label="Displacement (rad)")
plt.legend(loc="upper right")
plt.show()

"""