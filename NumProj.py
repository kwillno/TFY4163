import numpy as np
import matplotli.pyplot as plt
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
t_i = =
dt = 0.01



def RK4_step(k, f, theta, w, dt):
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
    k1 = k(theta)
    f1 = f(w)
    k2 = k(theta + (dt/2)*f1)
    f2 = f(w + (dt/2)*k1)
    k3 = k(theta + (dt/2)*f2)
    f3 = f(w + (dt/2)*k2)
    k4 = k(theta + dt*f3)
    f4 = f(w + dt*k3)
    return theta + (dt/6)*(f1 + (2*f2) + (2*f3) + f4), w + (dt/6)*(k1 + (2*k2) + (2*k3) + k4)

def RK4_method(k, f, theta_0, omega_0, dt):
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
    t = np.linspace(0,T,T/dt)
    theta = np.zeros(len(t))
    omega = np.zeros(len(t))
    
    
    theta[0],omega[0] = theta_0,omega_0
    
    for i in range(1,len(t)):
        theta[i],omega[i] = RK4_step(k, f, theta[i-1], omega[i-1], dt)
    
    return theta, omega, t

theta_RK4,omega_RK4,t = RK4_method(k, f, theta_0, omega_0, dt)

plt.figure("RK4")
plt.title("RK4 Method")
plt.plot(t,theta_RK4,label="Displacement (rad)")
plt.legend(loc="upper right")
# plt.show()