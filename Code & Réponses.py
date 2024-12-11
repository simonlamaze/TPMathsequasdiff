#TP Equations Différentielles
import numpy as np 
import matplotlib.pyplot as plt
import random as rd
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib notebook
plt.rcParams["figure.figsize"] = (14,7)
#Simon Lamaze - Corto Beck

# Question 1a - On dérive la norme au carré et hop c'est constant = 1

# Question 1b - Juste on écrit le produit vectoriel 

# Question 2



Omega_x=np.array([[0,0,0],[0,0,-1],[0,1,0]])
Omega_z=np.array([[0,-1,0],[1,0,0],[0,0,0]])

@np.vectorize
def norme(psi):
    A = psi**2
    B= A.sum()
    return np.sqrt(B)


def projete(psi):
    return psi/norme(psi)

def euler_explicit(psi_0,t_0,t_f,E,delta,dt):
    
    T = t_f - t_0
    u = lambda t: (1-np.cos(2*np.pi*t/T))*np.cos(E*t+ np.sin(np.pi *t/T)/(np.pi/T))
    temps = np.arange(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    for i in range(1,len(Psi)):
        Psi[ i]= Psi[i-1]+ dt*(E*Omega_z + delta*u(dt*(i-1))*Omega_x)@Psi[i-1]
    plt.plot(temps, norme(Psi))
    return Psi 

def euler_implicit(psi_0,t_0,t_f,E,delta,dt):
    T = t_f - t_0
    u = lambda t: (1-np.cos(2*np.pi*t/T))*np.cos(E*t+ np.sin(np.pi *t/T)/(np.pi/T))
    A = lambda t: np.eye((3,3)) + dt(u(dt*i)*delta*Omega_x +E*Omega_z)
    temps = np.arange(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    epsilon = 0.00001
    for i in range(1,len(Psi)):
        P =Psi[i-1]
        psi_a = P
        psi_b= A(dt*(i-1))@psi_a
        while norme(psi_a -psi_b)>=epsilon : 
            psi_a =  psi_b
            psi_b = A(dt*i)@psi_a
        Psi[i]= psi_b
    plt.plot(temps, norme(Psi))
    return Psi

def euler_projete(psi_0,t_0,t_f,E,delta,dt):
    T = t_f - t_0
    u = lambda t: (1-np.cos(2*np.pi*t/T))*np.cos(E*t+ np.sin(np.pi *t/T)/(np.pi/T))
    A = lambda t: np.eye(3) + dt*(u(dt*i)*delta*Omega_x +E*Omega_z)
    temps = np.arange(0,T, dt)
    Psi=np.zeros((len(temps),3))
    Psi[0]=psi_0
    epsilon = 0.00001
    for i in range(1,len(Psi)):
        P =Psi[i-1]
        psi_a = P
        psi_b= A(dt*(i-1))@psi_a
        while norme(psi_a -psi_b)>=epsilon : 
            psi_a =  psi_b
            psi_b = A(dt*i)@psi_a
        Psi[i]= projete(psi_b)

    plt.plot(temps, norme(Psi))

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
def main(deltamin, deltamax ,alpha):
    psi_0=np.array([0,0,-1])
    T=50
    dt=1e-3 
    E= 2
    erreurs= np.array(10)
    for i in range (10):
        deltai=rd.uniform(deltamin, deltamax)
        Ei=rd.uniform(E-alpha, E+alpha)
        Psi = euler_projete(psi_0,0,T,E,deltai,dt)
        erreurs[i]= norme (Psi[len(Psi)]-np.array([0,0,1]))
    return erreurs.mean()

#Question 4 

#a) Thm d'unicité des solutions
#b) On dérive
#c)

def proj_so3 (M): # M une matrice
     U , sigma , V = np.linalg.svd(M)
     return U@(V.t)

def euler_so3 ( psi_0,t_0,t_f,E,delta,dt):
    T = t_f - t_0
    u = lambda t: (1-np.cos(2*np.pi*t/T))*np.cos(E*t+ np.sin(np.pi *t/T)/(np.pi/T))
    A = lambda t: np.eye(3) + dt*(u(dt*i)*delta*Omega_x +E*Omega_z)
    temps = np.arange(0,T, dt)

    U=np.zeros((len(temps),3))

    U[0]=np.eye(3)
    epsilon = 0.00001
    for i in range(1,len(U)):
        u =U[i-1]
        u_a = u
        u_b= A(dt*(i-1))@u_a
        while norme((u_a -u_b).flatten())>=epsilon : 
            u_a =  u_b
            u_b = A(dt*i)@u_a
        U[i]= proj_so3(u_b)

    plt.plot(temps, norme(U@psi_0)) # méthode implicite avec la projection

 #d)  Les preuves vont bien (j'y puis pas encore arrivé)  + structure de group de SO3