from IPython import get_ipython
get_ipython().magic('reset -sf')

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
v, rope uptake speed (m/s)
ve, effective engine speed (m/s)
Sth, throttle setting
Me, pilot controlled moment (torque) from elevator
'''

import os, subprocess, sys, time
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tanh
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend
from scipy.integrate import odeint
g = 9.8
close("all") 


class glider:
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
        # algebraic functions
        self.L = 0
        self.D = 0
        self.M = 0
        # state variables 
        self.x = 0
        self.xD = 0
        self.y = 0    
        self.yD = 0
        self.theta = 0
        self.thetaD = 0   
        return
#    def L(self):
#        return (self.W + self.Lalpha*self.alpha) * (v/vb)^2
    
class rope:
    def __init__(self):
        # Rope parameters  
        Rrope = 0.005/2.0     #  rope radius (m)
        self.A = pi*Rrope**2      #  rope area (m2)
        self.Y = 3e9             #  2400*9.8/pi*(0.005/2)^2/0.035  
                                 #  effective Young's modulus 10 GPa for rope from Dyneema
                                 #                 datasheet 3.5% average elongation at break,  
                                 #                 average breaking load of 2400 kg (5000 lbs)
        self.a = 0.2             #  horizontal distance (m) of rope attachment in front of CG
        self.b = 0.1             #  vertial distance (m) of rope attachment below CG
        self.lo = 1000           #  initial rope length (m)
        # state variables 
        self.T = 0
        
class winch:
    def __init__(self):
        # Winch parameters  
        self.me = 40              #  Winch effective mass, kg
        self.rdrum = 0.4          # Drum radius (m), needed only for conversion of rpm to engine peak effective speed
        #state variables
        self.v = 0            # Rope uptake speed


class torqconv:
     def __init__(self):
         # TC parameters  
         self.Ko = 12             #TC capacity (rad/sec/(Nm)^1/2)
         self.dw = 0.13
     def invK(self,vrel):
         return 1/float(self.Ko) * tanh((1-vrel)/self.dw)

class engine:
    def __init__(self,rdrum):
        # Engine parameters  
        self.hp = 350             # engine horsepower
        self.Pmax = 750*self.hp        # engine watts
        self.rpmpeak = 6000       # rpm for peak power
#        self.vpeak = self.rpmpeak*2*pi/60*rdrum   #engine effectivespeed for peak power 
        self.vpeak = 30   # m/s engine effectivespeed for peak power.  This is determined by gearing, not the pure engine rpms
        self.me = 10.0            #  Engine effective mass (kg), effectively rotating at rdrum
        self.deltaEng = 1         #  time delay (sec) of engine power response to change in engine speed
        self.pe1 = 1.0; self.pe2 = 1.0; self.pe3 = 1.0 #  engine power curve parameters, gas engine
        # state variables 
        self.v = 0            #engine effective speed (m/s)
        self.Few = 0           # effective force between engine and winch (could go in either engine or winch or in its own class)
    
    def Pavail(self,ve):            # power curve
        vr = ve/float(self.vpeak)
        return self.Pmax*(self.pe1 * vr + self.pe2 * (vr)**2 - self.pe3 * (vr)**3)
        
class operator:
    def __init__(self):
        return     
    def control(self,t,gl,rp,wi,en):
        dStr = 0
        return dStr
         
class pilot:
    def __init__(self):
        self.Me = 0
        return
    def control(self,t,gl):
        dMe = 0
        return dMe
    
class plots:
    def __init__(self):
        return 
    def xy(self,x,ys,xlbl,ylbl,legendLabels,title):
        # plot results
        figure()
        for iy,y in enumerate(ys):
            plot(x,y,label = legendLabels[iy])
            xlabel(xlbl)
            ylabel(ylbl)
        legend()
        show()
        
def stateSplitMat(S,gl,rp,wi,tc,en,op,pl):
    '''Splits the formal state matrix S (each row a different time) into the various state variables'''
    gl.x      = S[:,0]
    gl.xD     = S[:,1]
    gl.y      = S[:,2]
    gl.yD     = S[:,3]
    gl.theta  = S[:,4]
    gl.thetaD = S[:,5]
    rp.T      = S[:,6]
    wi.v      = S[:,7]
    en.v      = S[:,8]
    op.Sth    = S[:,9]
    pl.Me     = S[:,10]
    return gl,rp,wi,tc,en,op,pl
    
    
def stateSplitVec(S,gl,rp,wi,tc,en,op,pl):
    '''Splits the formal state matrix S (each row a different time) into the various state variables'''
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
    pl.Me     = S[10]
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
    
    gl,rp,wi,tc,en,op,pl = stateSplitVec(S,gl,rp,wi,tc,en,op,pl)
    if 1.5<t<2.0:
        print 'pause'

    if gl.xD < 1e-6:
        gl.xD = 1e-6 #to handle v = 0 initial
    if en.v < 1e-6:
        en.v = 1e-6 #to handle v = 0 initial
    #----algebraic functions----#
    #rope
    thetarope = arctan(gl.y/float(rp.lo-gl.x))
    lenrope = sqrt((rp.lo-gl.x)**2 + gl.y**2)
    #glider
    v = sqrt(gl.xD**2 + gl.yD**2) # speed
    gamma = arctan(gl.yD/gl.xD)  # climb angle.  1d-6 handles the starting case of v=0
    alpha = gl.theta - gamma # angle of attack
    gl.L = (gl.W + gl.Lalpha*alpha) * (v/gl.vb)**2
    gl.D = gl.L/float(gl.Q)*(1 + gl.dv*(v/gl.vb-1)**2 + gl.dalpha*alpha**2) #drag
    gl.M = (-gl.palpha*alpha - pl.Me) #torque on glider
    vgw = (gl.xD*(rp.lo - gl.x) - gl.yD*gl.y)/float(lenrope) #velocity of glider toward winch
    #winch-engine
    vrel = wi.v/en.v
    Fee = tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3     
    Few = Fee* (2-vrel)  # effective force between engine and winch through torque converter
    #----derivatives of state variables----#
    dotx = gl.xD    
    dotxD = 1/float(gl.m) * (rp.T*cos(thetarope) - gl.D*cos(gamma) - gl.L*sin(gamma)) #x acceleration
    doty = gl.yD
    if gl.y < 0.1 and gl.L < gl.W: #on ground 
        dotyD = 1/float(gl.m) * (gl.L*cos(gamma) - rp.T*sin(thetarope) - gl.D*sin(gamma)) #y acceleration on ground
    else:
        dotyD = 1/float(gl.m) * (gl.L*cos(gamma) - rp.T*sin(thetarope) - gl.D*sin(gamma) - gl.W) #y acceleration
    dottheta = gl.thetaD    
    dotthetaD = 1/float(gl.Ig) * (rp.T*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/float(rp.a))-gl.theta-thetarope) + gl.M)    
    dotT = rp.Y*rp.A*(wi.v - vgw)/float(lenrope)
    dotvw =  1/float(wi.me) * (Few - rp.T)
    dotve =  1/float(en.me) *(op.Sth * en.Pavail(en.v) / float(en.v) - Few / (2 - vrel))
    dotSth = op.control(t,gl,rp,wi,en)
    dotMe = pl.control(t,gl) 
    return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotT,dotvw,dotve,dotSth,dotMe]

##########################################################################
#                         Main script
##########################################################################                        
tStart = 0
tEnd = 20      # end time for simulation
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
S = odeint(stateDer,S0,t,args=(gl,rp,wi,tc,en,op,pl))
#Split S (now a matrix with a state row for each time)

gl,rp,wi,tc,en,op,pl = stateSplitMat(S,gl,rp,wi,tc,en,op,pl)

# plot results
plts = plots()
plts.xy(t,[gl.x,gl.y],'time (sec)','position (m)',['x','y'],'glider position vs time')
plts.xy(t,[en.v,wi.v],'time (sec)','effective speed (m/s)',['engine','winch'],'Engine and winch speeds')


#figure()
#plot(t,gl.x)
#plot(t,gl.y)
#xlabel('time (sec)')
#ylabel('position')
#show()
#figure()
#plot(t,en.v)
#plot(t,wi.v+10)
#xlabel('time (sec)')
#ylabel('speed')
#show()

print 'Done'





