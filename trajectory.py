import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import time
import math
import matplotlib.animation as animation 

# Set common figure parameters
fig_size = 10
newparams = {'axes.grid': True,
             'lines.linewidth': 1.5, 'lines.markersize': 10,
             'font.size': 14}
plt.rcParams.update(newparams)
def Veff(rho, l):
    """ Evaluates the effective potential in the Schwarzschild geometry. """
    return -1/(2*rho) + l**2/(2*rho**2) - l**2/(2*rho**3)

def Vclassical(rho, l):
    """ Evaluates the classical effective potential. """
    return -1/(2*rho) + l**2/(2*rho**2)

def K(Z):
    """ Evaluates the kinetic part of the energy density. 
    Z = (x, y, u, v). 
    """
    rho2 = Z[0, :]**2 + Z[1, :]**2
    return .5 * (Z[3, :]*Z[1, :] + Z[2, :]*Z[0, :])**2/rho2
    
def getB(x, y, u, v):
    """ Computes the constant B. """
    l2 = (x*v - y*u)**2
    # Note that l2 has the "dimension" R_s 
    return 3*l2

def getA():
    """ Computes the constant A. """
    return .5

def RHS(Z, A, B):
    """ Returns the time derivatives of Z = X, Y, U and V. """
    
    rho = np.sqrt(Z[0]**2 + Z[1]**2)
    correction = 1 + B/rho**2
    dUdtau = -A*Z[0]/rho**3*correction
    dVdtau = -A*Z[1]/rho**3*correction
    
    return np.array([Z[2], Z[3], dUdtau, dVdtau])

def rk4step(f, y, h, A, B):
    """Calculate next step of an IVP of an ODE with a RHS described by the RHS function with RK4.
    Parameters:
        f: function. Right hand side of the ODE.
        y: float. Current step (position).
        h: float. Step-length.
    Returns:
        Next step.
    """
    
    s1 = f(y, A, B)
    s2 = f(y + h*s1/2.0, A, B)
    s3 = f(y + h*s2/2.0, A, B)
    s4 = f(y + h*s3, A, B)
    
    return y + h/6.0*(s1 + 2.0*s2 + 2.0*s3 + s4)

def getOrbit(n, T_max, Z0):
    """ Computes the orbit of the particle using the fourth order
    Runge-Kutta method.
    Parameters:
        n     : int. Number of iterations
        T_max : float. Stop time T=tau/t0
        Z0    : numpy-array, len(4), float.
                Position and velocities
                Z0 = [x, y, u, v]
    Returns:
        numpy-array, size(4, n). Z = [x[], y[], u[], v[]]
        Position and velocities for the orbit.
    """
    B = getB(*Z0)
    print("GR correction constant: %.2e"%(B))
    A = getA()
    
    h = T_max/float(n)
    Z = np.zeros((4, n + 1))
    Z[:, 0] = Z0

    tic = time.time()
    for i in range(0, n):
        Z[:, i + 1] = rk4step(RHS, Z[:, i], h, A, B)
        
    print("%.5f s, run time with %i steps."% (time.time() - tic, n))
    
    return Z
'''
import matplotlib.animation as animation 
#from matplotlib import animation, rc
from IPython.display import HTML

import matplotlib.pyplot as plt
plt.rcParams["animation.html"] = "jshtml"
from matplotlib import rc, animation
import matplotlib
matplotlib.rcParams['animation.embed_limit'] = 2**128

rc('animation', html='html5')
'''
def plotOrbit(Z, lim_fact):
    """ Plots the orbit, energy and effective potential.
    Parameters:
        Z: numpy-array, size(4, n). Z = [x[], y[], u[], v[]]
           Position and velocities for the orbit.
        lim_fact: float. Axis limits are given lim_fact
           multiplied by the start position.
    """

    xdata, ydata = [], [] 
    def init(): 
        # creating an empty plot/frame 
        line.set_data([], []) 
        return line, 
    
     
    def animate(i): 
        # appending new points to x, y axes points list 
        
        xdata.append(Z[0,i]) 
        ydata.append(Z[1,i]) 
        
        # set/update the x and y axes data 
        line.set_data(xdata, ydata) 
      
        # return line object 
        return line,   

    
    
    fig=plt.figure(figsize=(2*fig_size, fig_size))
        
    rho = (Z[0, :]**2 + Z[1, :]**2)**.5
    l = Z[0, :]*Z[3, :] - Z[1, :]*Z[2, :]
    #print(l)
    
    # Trajectory
    ax = plt.subplot(2, 2, (1, 3))
    plt.title("Orbit")
    #ax.plot(Z[0, :], Z[1, :], label="Orbit")
    ax.plot(Z[0, 0], Z[1, 0], "o", label="Start")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    ax.plot(0, 0, ".", label="Origin")
    circle = plt.Circle((0, 0), radius=1, color="gray")
    ax.add_artist(circle)
    ax_lim = max(Z[0, 0], Z[1, 0])
    plt.xlim([-lim_fact*ax_lim, lim_fact*ax_lim])
    plt.ylim([-lim_fact*ax_lim, lim_fact*ax_lim])
    plt.legend()
    
    
    line, = ax.plot([], [], lw=1) 
    
    anim = animation.FuncAnimation(fig, animate, init_func=init,repeat=False, 
                               frames=5000, interval=10, blit=True)
    #anim.save('scaterring.mp4', writer="ffmpeg")
    plt.show()
    
    
def plot2Orbit(Z1, Z2, lim_fact):
    """ Plots the orbit, energy and effective potential.
    Parameters:
        Z: numpy-array, size(4, n). Z = [x[], y[], u[], v[]]
           Position and velocities for the orbit.
        lim_fact: float. Axis limits are given lim_fact
           multiplied by the start position.
    """

    xdata1, ydata1 = [], [] 
    xdata2, ydata2 = [], [] 
    # center of mass
    xdata3, ydata3 = [], [] 
    
    def init(): 
        # creating an empty plot/frame 
        line1.set_data([], []) 
        line2.set_data([], []) 
        line3.set_data([], []) 
        return line1,line2,line3, 
    
     
    def animate(i): 
        # appending new points to x, y axes points list 
        #distance between particles
        distance  = math.sqrt((Z1[0,i]-Z2[0,i])**2 +(Z1[1,i]-Z2[1,i])**2)
        #angular momentum x*v-y*u for two particles
        
        rho1_i = (Z1[0, :]**2 + Z1[1, :]**2)**.5
        rho2_i = (Z2[0, :]**2 + Z2[1, :]**2)**.5
        l1_i = Z1[0, :]*Z1[3, :] - Z1[1, :]*Z1[2, :]
        l2_i = Z2[0, :]*Z2[3, :] - Z2[1, :]*Z2[2, :]
        Kinetic1=K(Z1)
        Kinetic2=K(Z2)
        print("angular momentum l1: ",l1_i[i]," and l2:",l2_i[i])
        
        #print(rho1_i[i],",",Kinetic1[i],",",rho1_i[i],",",l1_i[i],",",rho2_i[i],",",Kinetic2[i],",",rho2_i[i],",",l2_i[i])
        #print(i,",",np.sqrt(Z1[3, :]**2+Z1[2, :]**2)[i],",",Kinetic1[i],",",i,",",l1_i[i],",",i,",",np.sqrt(Z2[3, :]**2+Z2[2, :]**2)[i],",",Kinetic2[i],",",i,",",l2_i[i],",",i,",",(Kinetic1[i]+Kinetic2[i])/2)
        
        
        l1 = (Z1[0,i]*Z1[3,i] - Z1[1,i]*Z1[2,i])
        l2 = (Z2[0,i]*Z2[3,i] - Z2[1,i]*Z2[2,i])
        
        #print("Step: ",i," Angular momentum l1: ",l1," and l2:",l2)
        
        #print(distance)
        
        if distance < 0.01:
            print("Collision!")
        xdata1.append(Z1[0,i]) 
        ydata1.append(Z1[1,i]) 
        
        xdata2.append(Z2[0,i]) 
        ydata2.append(Z2[1,i]) 
        
        xdata3.append((Z1[0,i]+Z2[0,i])/2)
        ydata3.append((Z1[1,i]+Z2[1,i])/2)
        
        # set/update the x and y axes data 
        line1.set_data(xdata1, ydata1) 
        line2.set_data(xdata2, ydata2) 
        line3.set_data(xdata3, ydata3) 
      
        # return line object 
        return line1,line2,line3,   

    
    
    fig=plt.figure(figsize=(2*fig_size, fig_size))
    '''    
    rho1 = (Z1[0, :]**2 + Z1[1, :]**2)**.5
    l1 = Z1[0, :]*Z1[3, :] - Z1[1, :]*Z1[2, :]
    
    rho2 = (Z2[0, :]**2 + Z2[1, :]**2)**.5
    l2 = Z2[0, :]*Z2[3, :] - Z2[1, :]*Z2[2, :]
    print("Angular momentum l1: ",l1," and l2:",l2)
    '''
    #print(l)
    
    # Trajectory
    ax = plt.subplot(2, 2, (1, 3))
    plt.title("Orbit")
    #ax.plot(Z[0, :], Z[1, :], label="Orbit")
    ax.plot(Z1[0, 0], Z1[1, 0], "o", label="Start Z1")
    ax.plot(Z2[0, 0], Z2[1, 0], "o", label="Start Z2")
    #ax.plot(xdata3,ydata3, "o", label="Center of Mass")
    
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    ax.plot(0, 0, ".", label="Origin")
    circle = plt.Circle((0, 0), radius=1, color="gray")
    ax.add_artist(circle)
    ax_lim = max(Z1[0, 0], Z1[1, 0])
    plt.xlim([-lim_fact*ax_lim, lim_fact*ax_lim])
    plt.ylim([-lim_fact*ax_lim, lim_fact*ax_lim])
    plt.legend()
    
    
    line1, = ax.plot([], [], lw=1) 
    line2, = ax.plot([], [], lw=1) 
    line3, = ax.plot([], [], lw=1) 
    '''
    before_collision=0
    for i in range(5000):
        distance  = math.sqrt((Z1[0,i]-Z2[0,i])**2 +(Z1[1,i]-Z2[1,i])**2)
        if distance < 0.01:
            print("Collision at index:",i)
            before_collision=i
            break
    '''
    anim = animation.FuncAnimation(fig, animate, init_func=init,repeat=False, 
                               frames=5000, interval=1, blit=True)
    #anim.save('radial3.mp4', writer="ffmpeg")
    plt.show()
    
    rho1 = (Z1[0, :]**2 + Z1[1, :]**2)**.5
    l1 = Z1[0, :]*Z1[3, :] - Z1[1, :]*Z1[2, :]
    rho2 = (Z2[0, :]**2 + Z2[1, :]**2)**.5
    l2 = Z2[0, :]*Z2[3, :] - Z2[1, :]*Z2[2, :]
    
# Effective potential particle 1
    plt.subplot(2, 2, 1)
    plt.title("Effective potential particle 1")
    r1 = np.linspace(1.5, 20, 5000)
    plt.plot(r1, Veff(r1, l1[0]), label="Eff pot.")
    e1 = Veff(rho1[0], l1) + K(Z1[:, 0:2])[0]
    plt.plot([r1[0], r1[-1]], [e1, e1], label="Energy density")
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"$V$")
    
    
    # Energy per unit rest mass particle 1
    plt.subplot(2, 2, 2)
    plt.title("Kinetic and potential energy particle 1")
    V1 = Veff(rho1, l1)
    Ek1 = K(Z1)
    E1 = V1 + Ek1
    plt.plot(V1)
    plt.plot(Ek1, label=r"Kinetic $\frac{1}{2}\left(\frac{dr}{d\tau}\right)^2$")
    plt.plot(E1, label="Total En.")
    plt.plot(V1, label="Pot. En.")
    plt.xlabel("Step")
    plt.ylabel("Energy density")
    plt.legend()
    print("Relative change in\n E: %.2e\n l: %.2e"%((E1[0] - E1[-1])/E1[0], (l1[0] - l1[-1])/l1[-1]))
    #plt.savefig('precessing.png')
    plt.show()


# Effective potential particle 2
    plt.subplot(2, 2, 1)
    plt.title("Effective potential particle 2")
    r2 = np.linspace(1.5, 20, 5000)
    plt.plot(r2, Veff(r2, l2[0]), label="Eff pot.")
    e2 = Veff(rho2[0], l2) + K(Z2[:, 0:2])[0]
    plt.plot([r2[0], r2[-1]], [e2, e2], label="Energy density")
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"$V$")
    
    
    # Energy per unit rest mass particle 1
    plt.subplot(2, 2, 2)
    plt.title("Kinetic and potential energy particle 2")
    V2 = Veff(rho2, l2)
    Ek2 = K(Z2)
    E2 = V2 + Ek2
    plt.plot(V2)
    plt.plot(Ek2, label=r"Kinetic $\frac{1}{2}\left(\frac{dr}{d\tau}\right)^2$")
    plt.plot(E2, label="Total En.")
    plt.plot(V2, label="Pot. En.")
    plt.xlabel("Step")
    plt.ylabel("Energy density")
    plt.legend()
    print("Relative change in\n E: %.2e\n l: %.2e"%((E2[0] - E2[-1])/E2[0], (l2[0] - l2[-1])/l2[-1]))
    #plt.savefig('precessing.png')
    plt.show()




#radial
#Z0 = [0, 10, .1845, 0]

#Z1 = [10, -20, -0.0375, 0]
#Z2 = [0, 5, 0.4, 0]



#Z1 = [3, 10, -0.1375, -0.1]
#Z2 = [-3,10, 0.1375, -0.1]

Z1 = [0, 10, 0.2, 0]
Z2 = [0, 10, -0.2, 0]

#Z0 = [0, 10, .1849, 0]
n = 5000
tau_max = 520.0
Z11 = getOrbit(n, tau_max, Z1)
Z22 = getOrbit(n, tau_max, Z2)

plot2Orbit(Z11, Z22, 1.1)


'''
#plt.show()
#precessing
#Z0 = [0, 20, .1, 0]
#Z0 = [0, 10, 0.2, -0.2]
#Z0 = [0, 10, .25, 0]
#Z0 = [0, 10, 0.2, -.1]
Z0 = [0, 10, .2, 0]
n = 5000
tau_max = 5000
Z = getOrbit(n, tau_max, Z0)
plotOrbit(Z, 1.1)
'''
'''
#scaterring
#Z0 = [0, 100, 0.05, -0.5]
Z0 = [0, 10, 0.2, -.25]

n = 5000
tau_max = 1000
Z = getOrbit(n, tau_max, Z0)
plotOrbit(Z, 1.1)
'''
