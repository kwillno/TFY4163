# Husk å definere nødvendige biblioteker
import numpy as np
import matplotlib.pyplot as plt


# Definer konstanter og initialbetingelser
l = 1
m = 5
g = 9.81
theta_0 = 0.2
omega_0 = 0.0

# Tidskonstanter
T_i = 0
T_f = 10
T = T_f - T_i
samples = 10000
dt = T/samples

def euler_cromer_approx(theta_0, omega_0, dt, T_i):
    """
    Calculates angular displacement and angular velocity using the Euler-Cromer method
    
    theta_0: initial angular displacement 
    omega_0: initial angular velocity
    dt: timestep
    """
    theta=np.zeros(samples)
    omega=np.zeros(samples)
    t=np.zeros(samples)
    theta[0]=theta_0
    omega[0]=omega_0
    t[0]=T_i

    for n in range(1,samples):
        """
        Beregner ny verdi av omega ved hjelp av gammel theta og omega
        """
        omega[n]=omega[n-1]-g/l*theta[n-1]*dt
        
        """
        Beregner ny verdi av omega ved hjelp av nye omega og gammel theta
        """
        theta[n]=theta[n-1]+omega[n]*dt
        
        """
        Får et tid-array med samme dimensjon som tidligere arrays for plotting
        """
        t[n]=t[n-1]+dt
        
    
    """
    theta: array with values of angular displacement
    w: array with values of angular velocity
    t: array with time-values
    """
    return theta,omega,t

# Tilordner arrays til output-ene fra approksimeringsfunksjonen
thetaArr,omegaArr,timeArr = euler_cromer_approx(theta_0,omega_0,dt,T_i)

plt.figure("Test")
plt.plot(timeArr,thetaArr,label="Theta")
plt.plot(timeArr,omegaArr,label="Omega")
plt.legend()
plt.grid()
# plt.show()


# Oppgave 2
def energy_calculation(theta_0, omega_0, dt):
    """
    Regner ut total energi for systemet numerisk ved euler-cromer methoden.
    Inputs:
    theta_0: startverdi av theta (vinkelutslag)
    omega_0: startverdi av omega (vinkelhastighet)
    
    """
    samples = int(T/dt) # Finner sample-rate for valgt dt
    
    # Lager arrays med verdier vha Euler-cromer approksimasjon
    thetaArr, omegaArr, timeArr = euler_cromer_approx(theta_0,omega_0,dt,T_i)
    
    # Funksjonsuttrykk for total-energi
    energy_func = lambda m,l,omega,theta: (1/2)*m*(l**2)*(omega**2) + (1/2)*m*g*l*(theta**2)
    
    # Får å få et samme-dimensjonert tid-array
    t = np.linspace(T_i,T,samples)
    energy = np.zeros(samples)
    
    for i in range(len(t)):
        """
        Beregning av total-energi for hvert t-element i t-array
        """
        energy[i] = energy_func(m,l,omegaArr[i],thetaArr[i])
    
    
    E_total = energy

    return t, E_total

# Verdier av dt vi ønsker å plotte for
dt1 = 0.001
dt2 = 0.004
dt3 = 0.007

plt.figure(2)
plt.title("Total Energi")
plt.plot(energy_calculation(theta_0,omega_0,dt1)[0],energy_calculation(theta_0,omega_0,dt1)[1], label="dt = 0.001")
plt.plot(energy_calculation(theta_0,omega_0,dt2)[0],energy_calculation(theta_0,omega_0,dt2)[1], label="dt = 0.004")
plt.plot(energy_calculation(theta_0,omega_0,dt3)[0],energy_calculation(theta_0,omega_0,dt3)[1], label="dt = 0.007")
plt.xlabel("Time(s)")
plt.ylabel("Energy(J)")
plt.legend()
plt.grid()
# plt.show()

# Oppgave 3
def energy_diff(time,E_tot):
    """
    Plotter en periode av Total Energi graf og finner differanse i Total Energi i start og slutt av periode.
    Inputs:
    
    time: Tid-array av samme dimensjon som E_tot
    E_tot: Array med total energi.
    """
    E_tot_0 = E_tot[0]
    
    # Uttryk for tid ved slutten av periode
    T_p = 2*np.pi*np.sqrt(l/g)
    
    i = 0
    while time[i] < T_p:
        i += 1
    
    time_p = np.zeros(i)
    E_tot_p = np.zeros(i)
    for j in range(i):
        time_p[j] = time[j]
        E_tot_p[j] = E_tot[j]
    
    plt.figure("En periode Total Energi")
    plt.title("En periode Total Energi")
    plt.plot(time_p,E_tot_p,label="dt = 0.001")
    plt.grid()
    # plt.show()
    
    print("Delta E etter en periode: " + str(abs(E_tot_p[0]-E_tot_p[-1])))
        
energy_diff(energy_calculation(theta_0,omega_0,dt3)[0],energy_calculation(theta_0,omega_0,dt3)[1])


# Oppgave 1
def euler_cromer_method(theta_0, omega_0, dt, T_i):
    """
    Calculates angular displacement and angular velocity using the Euler-Cromer method
    
    theta_0: initial angular displacement 
    omega_0: initial angular velocity
    dt: timestep
    """
    theta=np.zeros(samples)
    omega=np.zeros(samples)
    t=np.zeros(samples)
    theta[0]=theta_0
    omega[0]=omega_0
    t[0]=T_i

    for n in range(1,samples):
        """
        Beregner ny verdi av omega ved hjelp av gammel theta og omega
        """
        omega[n]=omega[n-1]-g/l*np.sin(theta[n-1])*dt
        
        """
        Beregner ny verdi av omega ved hjelp av nye omega og gammel theta
        """
        theta[n]=theta[n-1]+omega[n]*dt
        
        """
        Får et tid-array med samme dimensjon som tidligere arrays for plotting
        """
        t[n]=t[n-1]+dt
        
    
    """
    theta: array with values of angular displacement
    w: array with values of angular velocity
    t: array with time-values
    """
    return theta,omega,t

# Tilordner arrays til output-ene fra approksimeringsfunksjonen
thetaArr,omegaArr,timeArr = euler_cromer_method(theta_0,omega_0,dt,T_i)

plt.figure("New Test")
plt.plot(timeArr,thetaArr,label="Theta")
plt.plot(timeArr,omegaArr,label="Omega")
plt.legend()
plt.grid()
# plt.show()

# Oppgave 2

# Nye startverdier

theta_0 = np.radians(15)
omega_0 = 0
dt = 0.001

thetaArr_a,omegaArr_a,timeArr_a = euler_cromer_approx(theta_0,omega_0,dt,T_i)
thetaArr_m,omegaArr_m,timeArr_m = euler_cromer_method(theta_0,omega_0,dt,T_i)

plt.figure("Diff 15deg")
plt.title("Diff 15deg")
plt.plot(timeArr_a,thetaArr_a,label="Approx")
plt.plot(timeArr_m,thetaArr_m,label="Method")
plt.legend()
plt.grid()
# plt.show()

theta_0 = np.radians(40)
omega_0 = 0
dt = 0.001

thetaArr_a,omegaArr_a,timeArr_a = euler_cromer_approx(theta_0,omega_0,dt,T_i)
thetaArr_m,omegaArr_m,timeArr_m = euler_cromer_method(theta_0,omega_0,dt,T_i)

plt.figure("Diff 40deg")
plt.title("Diff 40deg")
plt.plot(timeArr_a,thetaArr_a,label="Approx")
plt.plot(timeArr_m,thetaArr_m,label="Method")
plt.legend()
plt.grid()
# plt.show()

def equation(t, vals):
    """
    Funksjonen tar inn t, som er en float, og vals som er en array med to elementer
    i dette tilfellet, [val1, val2], og returnerer verdien av høyre side
    av ligningene.
    
    t: float, tid
    vals: array med verdi for theta og omega
    
    Return: array med høyre-side av ligningene, som [ligning1, ligning2]
    """
    theta,omega = vals[0],vals[1]
    
    # Uttrykk for dw
    dw = -g/l*np.sin(theta)
    
    # Uttryk for dtheta
    dtheta = omega
    
    return [dtheta, dw]

"""
Note:
scipy.integrate.solve_ivp requires that fun returns
an object which is of type array_like. 
An ordinary list is of this type (aswell as e.g. integers, floats and numpy arrays)
so we can return a list, and do not have to first convert into an np.array. 
"""

from scipy.integrate import solve_ivp as DestoryerOfOrdinaryDifferentialEquations

def RK45_method(RHS, theta_0, omega_0, t_1, dt):
    """
    Calculates theta and omega using the scipy.integrate.solve_ivp function (RK54)
    
    RHS: right hand side of differential equations
    t_1: time-value to calculate up to (e.g. 10 seconds)
    dt: timestep
    
    return:
    -------
    theta: array of theta values
    w: array of omega values
    t: timevalues
    """
    initialValues = [theta_0,omega_0]

    timeSpan = [0,t_1 + dt]

    timeArr = np.arange(0,t_1+dt,dt)
    
    solution = DestoryerOfOrdinaryDifferentialEquations(RHS,timeSpan,initialValues,method="RK45",t_eval=timeArr)
    
    theta = solution.y[0, : ]
    omega = solution.y[1, : ]
    times = solution.t
    
    
    return theta, omega, times


theta,omega,times = RK45_method(equation, 0.2, 0, 10, 0.01)
# PLOT VINKELUTSLAGET


plt.figure("Scipy")
plt.plot(times,theta,label="Theta")
plt.plot(times,omega,label="omega")
plt.show()