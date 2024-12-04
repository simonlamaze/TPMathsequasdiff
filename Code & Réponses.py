#TP Equations Différentielles
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib notebook
plt.rcParams["figure.figsize"] = (14,7)
#Simon Lamaze - Corto Beck

# Question 1a - On dérive la norme au carré et hop c'est constant = 1

# Question 1b - Juste on écrit le produit vectoriel 

# Question 2

E=2
delta=1 
T=50
u = lambda t: (1-np.cos(2*np.pi*t/T))*np.cos(E*t+ np.sin(np.pi *t/T)/(np.pi/T))
psi_0=np.array([0,0,-1])
dt=1e-3
Omega_x=np.array([[0,0,0],[0,0,-1],[0,1,0]])
Omega_z=np.array([[0,-1,0],[1,0,0],[0,0,0]])

@np.vectorize
def norme(psi):
    return np.sqrt((psi**2).sum())

def euler_explicit(psi_0,t_0,t_f,E,delta,dt):
    
    temps = np.linspace(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    for i in range(1,len(Psi)):
        Psi[ i]= Psi[i-1]+ dt*(E*Omega_z + delta*u(dt*(i-1))*OMega_x)@Psi[i-1]
    plt.plot(t, norme(Psi))
    return Psi 

def euler_implicit(psi_0,t_0,t_f,E,delta,dt):
    ...
    return ...

def euler_projete(psi_0,t_0,t_f,E,delta,dt):
    ...
    return ...

