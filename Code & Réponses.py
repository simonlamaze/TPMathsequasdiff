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

def projete(psi):
    return psi/norme(psi)

def euler_explicit(psi_0,t_0,t_f,E,delta,dt):
    T = t_f - t_0
    
    temps = np.linspace(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    for i in range(1,len(Psi)):
        Psi[ i]= Psi[i-1]+ dt*(E*Omega_z + delta*u(dt*(i-1))*OMega_x)@Psi[i-1]
    plt.plot(t, norme(Psi))
    return Psi 

def euler_implicit(psi_0,t_0,t_f,E,delta,dt):
    T = t_f - t_0
    temps = np.linspace(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    epsilon = 0.00001
    for i in range(1,len(Psi)):
        P =Psi[i-1]
        A= E*Omega_z + delta*u(dt*(i))*OMega_x
        psi_a = P
        psi_b= P + dt*A@psi_a
        while norme(psi_a -psi_b)>=epsilon : 
            psi_a =  psi_b
            psi_b = P + dt*A@psi_a
        Psi[i]= psi_b
    plt.plot(t, norme(Psi))
    return Psi

def euler_projete(psi_0,t_0,t_f,E,delta,dt):
    T = t_f - t_0
    temps = np.linspace(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    epsilon = 0.00001
    for i in range(1,len(Psi)):
        P =Psi[i-1]
        A= E*Omega_z + delta*u(dt*(i))*OMega_x
        psi_a = P
        psi_b= P + dt*A@psi_a
        while norme(psi_a -psi_b)>=epsilon : 
            psi_a =  psi_b
            psi_b = P + dt*A@psi_a
        Psi[i]= projete(psi_b)

    fig = plt.figure()

    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')

    #    defining all 3 axis
    z = Psi[:,0]
    x = Psi[:,1]
    y = Psi[:,2]


    ax.plot3D(x, y, z, 'pink')
    ax.set_title('Trajectoire 3D')
    plt.show()

    return Psi

# Question 3

alpha=0.5
deltamin=0.4
deltamax=3.5
erreurs= np.array(10)
for i in range (10):
    deltai=random.randint(deltamin, deltamax)
    Ei=random.randint(E-alpha, E+alpha)
    Psi = euler_projete(psi_0,0,T,E,deltai,dt)
    erreurs[i]= norme (Psi[len(Psi)]-np.array([0,0,1]))
print erreurs.mean()


