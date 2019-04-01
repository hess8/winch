''' Bret Hess, bret.hess@gmail.com 
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
v, rope uptake speed (m/s)
ve, effective engine speed (m/s)
Sth, throttle setting
Me, pilot controlled moment (torque) from elevator
'''
import os, subprocess, sys, time
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tanh
from matplotlib.pyplot import plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim
from scipy.integrate import odeint
g = 9.8

class glider():
    def __init__(self):
        # parameters
        self.vb = 30              #   speed of glider at best glide angle
        self.vStall = 22          # stall speed fully loaded (m/s)
        self.stallLoss = 0.25     # loss of lift (fraction) post stall, see https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20140000500.pdf
        self.wStall = rad(5)      #transition width for stall loss
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
        # algebraic functions
        self.L = None
        self.D = None
        self.M = None
        # state variables 
        self.x = None
        self.xD = None
        self.y = None    
        self.yD = None
        self.theta = None
        self.thetaD = None   
        return
#    def L(self):
#        return (self.W + self.Lalpha*self.alpha) * (v/vb)^2
    
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
        self.me = 40              #  Winch effective mass, kg
        self.rdrum = 0.4          # Drum radius (m), needed only for conversion of rpm to engine peak effective speed
        #state variables
        self.v = None            # Rope uptake speed


class torqconv():
     def __init__(self):
         # TC parameters  
         self.Ko = 12
         self.dw = 0.13
     def invK(self,vrel):
         return 1/self.Ko * tanh((vrel-1/self.dw))

class engine():
    def __init__(self,rdrum):
        # Engine parameters  
        self.hp = 350             # engine horsepower
        self.Pmax = 750*self.hp        # engine watts
        self.rpmpeak = 6000       # rpm for peak power
        self.vpeak = self.rpmpeak*2*pi/60*rdrum   #engine effectivespeed for peak power  
        self.me = 10            #  Engine effective mass (kg), effectively rotating at rdrum
        self.deltaEng = 1         #  time delay (sec) of engine power response to change in engine speed
        self.pe1 = 1.0; self.pe2 = 1.0; self.pe3 = 1. #  engine power curve parameters, gas engine
        # state variables 
        self.v = None            #engine effective speed (m/s)
        self.Few = None           # effective force between engine and winch (could go in either engine or winch or in its own class)
    
    def Pavail(self,ve):            # power curve
        return self.Pmax*(self.pe1 * ve/self.vpeak + self.pe2 * (ve/self.vpeak)**2 - self.pe3 * (ve/self.vpeak)**3)
        
class operator():
    def __init__(self):
        def control(t,gl,rp,wi,en):
            dStr = 0
            return dStr
         
class pilot():
    def __init__(self):
        self.Me = None
        def control(t,gl):
            dMe = 0
            return dMe
    
class plots():
    def __init__(self):
        
        return 
        
def stateSplit(S,gl,rp,wi,tc,en,op,pl):
    '''Splits the formal state vector S into the various state variables'''
    gl.x      = S[0]
    gl.xD     = S[1]
    gl.y      = S[2]
    gl.yD     = S[3]
    gl.theta  = S[4]
    gl.thetaD = S[5]
    rp.T      = S[6]
    wi.v      = S[7]
    en.v      = S[8]
    op.Sth    = S[9]
    pi.Me     = S[10]
    return gl,rp,wi,tc,en,op,pl
    
def stateJoin(S,gl,rp,wi,tc,en,op,pl):
    '''Joins the various state variables into the formal state vector S '''
    S[0] = gl.x
    S[1] = gl.xD
    S[2] = gl.y        
    S[3] = gl.yD       
    S[4] = gl.theta    
    S[5] = gl.thetaD   
    S[6] = rp.T 
    S[7] = wi.v
    S[8] = en.v   
    S[9] = op.Sth  
    S[10] = pl.Me    
    return S

def stateDer(S,t,gl,rp,wi,tc,en,op,pl):
    '''Differential equations that give the first derivative of the state vector'''
    gl,rp,wi,en,tc,op,pl = stateSplit(S,gl,rp,wi,tc,en,op,pl)
    #----algebraic functions----#
    #rope
    thetarope = arctan(gl.y/(rp.lo-gl.x))
    lenrope = sqrt((rp.lo-gl.x)**2 + gl.y**2)
    #glider
    v = sqrt(gl.xD**2 + gl.yD**2) # speed
    gamma = arctan(gl.yD/gl.xD)  # climb angle
    alpha = gl.theta - gamma # angle of attack
    gl.L = (gl.W + gl.Lalpha*alpha) * (v/gl.vb)^2
    gl.D = gl.L/gl.Q*(1 + gl.dv*(v/gl.vb-1)^2 + gl.dalpha*alpha^2) #drag
    gl.M = (-gl.palpha*alpha - pl.Me) #torque on glider
    vgw = (gl.xD*(rp.lo - gl.x) - gl.yD*gl.y)/lenrope #velocity of glider toward winch
    #winch-engine
    vrel = wi.v/en.v
    Fee = tc.invK(vrel) * vrel**2 / wi.rdrum      
    Few = Fee* (2-vrel)  # effective force between engine and winch through torque converter
    #----derivatives of state variables----#
    dotx = gl.xD    
    dotxD = rp.T*cos(thetarope) - gl.D*cos(gamma) - gl.L*sin(gamma) #x acceleration
    doty = gl.yD    
    dotyD = gl.L*cos(gamma) - rp.T*sin(thetarope) - gl.D*sin(gamma) - gl.W #y acceleration
    dottheta = gl.thetaD    
    dotthetaD = 1/gl.Ig * (rp.T*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/rp.a)-gl.theta-thetarope) + gl.M)    
    dotT = rp.Y*rp.A*(wi.v - vgw)/lenrope
    dotvw =  1/wi.me * (Few - rp.T)
    dotve =  op.Sth * en.Pavail(en.v) / en.v - Few / (2 - wi.v - en.v)
    dotSth = op.control(t,gl,rp,wi,en)
    dotMe = pl.control(t,gl) 
    return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotT,dotvw,dotve,dotSth,dotMe]

##########################################################################
#                         Main script
##########################################################################                        
tStart = 0
tEnd = 10      # end time for simulation
t = linspace(tStart,tEnd)

gl = glider()
rp = rope()
wi = winch()
tc = torqconv()
en = engine(wi.rdrum)
op = operator()
pl = pilot()
  
S0 = zeros(11)
op.Sth = 0.5  #throttle setting
S0 = stateJoin(S0,gl,rp,wi,tc,en,op,pl)
S = odeint(stateDer,S0,t,gl,rp,wi,tc,en,op,pl)

# function that returns dy/dt
def aux(y,t):
    return 3*array([y[1],y[0]])*t
    
    

def model(y,t):
    k = 0.3
    dydt = -k * y * aux(y,t)
    return dydt

# initial condition
y0 = array([5,1])

# time points
t = linspace(0,20)

# solve ODE
y = odeint(model,y0,t)
#print(y[0])
# plot results
plot(t,y[:,0])
plot(t,y[:,1])
xlabel('time')
ylabel('y(t)')
show()






