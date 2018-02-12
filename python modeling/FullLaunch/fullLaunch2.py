''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Solves for the motion of a glider launched by a winch with a springy rope. 
The winch is connected to an engine through a torque converter. 

State varables:
x, horizontal glider CG position
xD, horizontal glider CG velocity
y, vertical glider CG position
yD, vertical glider CG velocity
theta, glider pitch angle above horizon
thetaD, glider pitch rotation rate
T, rope tension
vu, rope uptake speed (m/s)
ve, effective engine speed (m/s)
Sth, throttle setting
Me, pilot controlled moment (torque) from elevator
'''
import os, subprocess, sys, time
from numpy import pi, array, zeros,linspace
from matplotlib.pyplot import plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim
from scipy.integrate import odeint
g = 9.8

class glider():
    def __init__(self):
        # parameters
        self.vb = 30              #   speed of glider at best glide angle
        self.m = 400
        self.W = self.m*9.8          #   weight (N)
        self.Q = 20             #   L/D
        self.alphas = 6         #   stall angle vs glider zero
        Co = 0.75             #   Lift coefficient {} at zero glider AoA
        self.Lalpha = 2*pi*self.W/Co
        self.Ig = 600            #   Glider moment of inertia, kgm^2
        self.palpha = 3.0        #   change in air-glider pitch moment (/rad) with angle of attack 
        self.dv = 3.0            #   drag constant ()for speed varying away from vb
        self.dalpha = 1.8        #   drag constant (/rad) for glider angle of attack away from zero
        
        # state variables 
        self.x = None
        self.xD = None
        self.y = None    
        self.xD = None
        self.theta = None
        self.thetaD = None   
        self.vu = None
        return
    
class rope():
    def __init__(self):
        # Rope parameters  
        Rrope = 0.005/2     #  rope radius (m)
        self.A = pi*Rrope**2      #  rope area (m2)
        self.Y = 3e9             #  2400*9.8/pi*(0.005/2)^2/0.035  
                                 #  effective Young's modulus 10 GPa for rope from Dyneema
                                 #                 datasheet 3.5% average elongation at break,  
                                 #                 average breaking load of 2400 kg (5000 lbs)
        self.a = 0.2             #  horizontal distance (m) of rope attachment in front of CG
        self.b = 0.1             #  vertial distance (m) of rope attachment below CG
        self.lo = 1000           #  initial rope length (m)
        # state variables 
        self.T = None
        
class winch():
    def __init__(self):
        # Winch parameters  
        self.me = 40*9.8               #  Winch effective mass, kg
        #state variables
        self.vu = None            # Rope uptake speed

class engine():
    def __init__(self):
        # Engine parameters  
        self.hp = 350             # engine horsepower
        self.Pmax = 750*self.hp        # engine watts
        self.rpmpeak = 6000       # rpm for peak power
        self.omegapeak = self.rpmpeak*2*pi/60   #engine speed for peak power  
        self.me = 10            #  Engine effective mass (kg), effectively rotating at rdrum
        self.deltaEng = 1         #  time delay (sec) of engine power response to change in engine speed
        self.pe1 = 1.0; self.pe2 = 1.0; self.pe3 = 1. #  engine power curve parameters, gas engine
        # state variables 
        self.ve = None
        
class operator():
    def __init__(self):
        return

class pilot():
    def __init__(self):
        return
    
class plots():
    def __init__(self):
        return
    
##########################################################################
#                         Main script
##########################################################################                        
tStart = 0
tEnd = 10      # end time for simulation

gl = glider()
rp = rope()
wi = winch()
en = engine()
op = operator()



def stateSplit(S,gl,rp,wi,en):
    '''Splits the formal state vector S into the various state variables'''
    gl.x      = S[0]
    gl.xD     = S[1]
    gl.y      = S[2]
    gl.yD     = S[3]
    gl.theta  = S[4]
    gl.thetaD = S[5]
    rp.T      = S[6]
    wi.vu     = S[7]
    en.ve     = S[8]
    op.Sth    = S[9]
    pi.Me     = S[10]

def stateDer(S,t):
    '''Differential equations that give the first derivative of the state vector'''
    
# S0 = zeroes(11)
# S0[0] = 
# S0[]

# function that returns dy/dt
def model(y,t):
    k = 0.3
    dydt = -k * y
    return dydt

# initial condition
y0 = 5

# time points
t = linspace(0,20)

# solve ODE
y = odeint(model,y0,t)

# plot results
plot(t,y)
xlabel('time')
ylabel('y(t)')
show()







