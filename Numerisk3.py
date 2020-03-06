import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def euler_cromer_approx(theta_0, w_0, dt):
    """
    Calculates angular displacement and angular velocity 
    using the Euler-Cromer method 
    """
    N = int(T/dt)
    theta = np.zeros(N+1)
    w = np.zeros(N+1)
    t = np.linspace(0, T, N + 1)
    theta[0] = theta_0
    w[0] = w_0
    for i in range(N):
        w[i+1] = w[i] - g/l*theta[i]*dt
        theta[i+1] = theta[i] + w[i+1]*dt
    return theta, w, t


# RK5(4) method

def equation(t, vals):
    """
    Calculates the value of RHS of the differentail equations given an array (vals) which 
    contains the values of the parameters

    Parameters
    -----------
    t: float
    time
    
    vals: array
    values of theta and omega, [theta, omega] 
    
    Returns
    -------
    [dtheta, dw]: list with values of the RHS of the equations.
    """
    dw = -g/l*vals[0]  # Regner ut endring i w. 
    dtheta = vals[1]
    return [dtheta, dw]

def RK45_method(RHS, theta_0, w_0, t_1, dt):
    """
    Calculates the angular dispacement and angular velocity.
    
    Parameters:
    ------------
    RHS: right hand side of the differentail equations
    theta_0: initial value of angular displacement
    w_0: initial value of the angular velocity
    t_1: time to calculate up to
    dt: timestep 
    
    returns:
    ---------
    
    theta: array of theta values
    w: array of omega values
    t: array of time values
    """
    
    init_values = [theta_0, w_0]
    t_span = [0, t_1+dt]
    t = np.arange(0, t_1 + dt, dt)
    theta12 = integrate.solve_ivp(RHS, t_span, init_values, method = 'RK45', t_eval = t)
    theta = theta12.y[0, :]
    w = theta12.y[1, :]
    t = theta12.t
    return theta, w, t

theta_analytic_func = lambda theta_0,t : theta_0*np.cos(np.sqrt(g/l)*t)

# Initial Parameters

l = 1.0
m = 5.0
g = 9.8
theta_0 = 0.2
omega_0 = 0.0

# Timing parameters

dt = 0.1
T = 10
t = np.linspace(0,T,dt)


theta_ec,omega_ec,t = euler_cromer_approx(theta_0,omega_0,dt)

theta_analytic = theta_analytic_func(theta_0,t)

theta_RK45,omega_RK45,times = RK45_method(equation, theta_0, omega_0, T, dt)

plt.figure("Theta")
plt.title("Angular displacement")
plt.plot(t,theta_analytic,label="Analytic")
plt.plot(t,theta_ec,label="Euler Cromer")
plt.plot(t,theta_RK45,label="Runge-Kutta 5(4)")
plt.legend()
# plt.show()

plt.figure("Theta err from analytic")
plt.title("Error in Angular displacement")
# plt.plot(t,theta_analytic-theta_analytic,label="Analytic")
plt.plot(t,theta_ec-theta_analytic,label="Euler Cromer")
plt.plot(t,theta_RK45-theta_analytic,label="Runge-Kutta 5(4)")
plt.legend()
# plt.show()

k = lambda theta : -g/l*theta

f = lambda omega : omega

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
    k: RHS of equation for theta
    f: RHS of equation for omega
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

# Initial physical parameters

l = 1.0
m = 5.0
g = 9.8
theta_0 = 0.2
omega_0 = 0.0

# Timing parameters

T = 20 
dt = 0.01

theta_RK4,omega_RK4,t = RK4_method(k, f, theta_0, omega_0, dt)

plt.figure("RK4")
plt.title("RK4 Method")
plt.plot(t,theta_RK4,label="Displacement (rad)")
plt.legend(loc="upper right")
# plt.show()

E_kin = lambda omega : (1/2)*m*l**2*omega**2
E_pot = lambda theta : (1/2)*m*g*l*theta**2

E_tot = lambda theta,omega: E_kin(omega) + E_pot(theta)

theta_RK4, omega_RK4,t = RK4_method(k, f, theta_0, omega_0, dt)

E_t = E_tot(theta_RK4,omega_RK4)
E_k = E_kin(omega_RK4)
E_p = E_pot(theta_RK4)

plt.figure("Energy RK4")
plt.title("Total energy calculation using RK4 Method")
plt.plot(t,E_t,label="Total energy (J)")
plt.legend()
# plt.show()

plt.figure("Energy RK4")
plt.title("Total energy calculation using RK4 Method")
plt.plot(t,E_p,label="Potential energy (J)")
plt.plot(t,E_k,label="Kinetic energy (J)")
plt.legend(loc="upper right")
# plt.show()

dt_i = 0.0001
dt_f = 0.1
dt_delta = 0.01
T = 20

dt = np.linspace(dt_i,dt_f,T/dt_delta)
E_total = np.zeros(len(dt))

for i in range(0,len(dt)):
    theta_RK4, omega_RK4,t = RK4_method(k, f, theta_0, omega_0, dt[i])
    
    E_k_i = E_kin(omega_RK4[0])
    E_p_i = E_pot(theta_RK4[0])
    
    E_k_f = E_kin(omega_RK4[-1])
    E_p_f = E_pot(theta_RK4[-1])
    
    E_total[i] = abs((E_k_f + E_p_f) - (E_k_i + E_p_i))
    
plt.figure("Energy RK4 sfa dt")
plt.title("Total energy (RK4) with variable dt")
plt.plot(dt,E_total,label="Energy (J)")
plt.legend()
plt.show()