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

from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tanh,ceil,floor
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend
from scipy.integrate import odeint
g = 9.8
 

class timesteps:
    def __init__(self,tStart,tEnd,N):
        '''The way to do persistent variables in python is with class variables
        
        The solver does not keep a constant time step, and evaluates the stateDer 
        function usually two or more steps per time step'''
        
        self.oldt = -1.0
        self.i = -1   #counter for time step
        self.dt = (tEnd -tStart)/float((N-1))
    

class glider:
    def __init__(self,ntime):
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
        self.x = 0
        self.xD = 0
        self.y = 0    
        self.yD = 0
        self.theta = 0
        self.thetaD = 0 
        self.data = zeros(ntime,dtype = [('t', float),('x', float),('xD', float),('y', float),('yD', float),\
                                    ('v', float),('theta', float),('alpha', float)])
        return
    
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
    def __init__(self,ntime):
        self.Sth = 0
        self.data = zeros(ntime,dtype = [('t', float),('Sth', float)])
        return
    def control(self,t,gl,rp,wi,en):
        tramp = 4   #seconds
        thrmax = 0.5        
        if t < tramp:
            self.Sth =  thrmax/float(tramp) * t
        else:
            self.Sth =  thrmax
         
class pilot:
    def __init__(self,ntime):
        self.Me = 0
        self.ctrltype = ''
        self.data = zeros(ntime,dtype = [('t', float),('err', float),('Me', float)])
        return
        
    def control(self,t,ts,gl): 
        setpoint = 0.0
        tint = 0.5 #sec
        pp = 5; pd = 0; pint = 0
        Nint = ceil(tint/ts.dt)
        ctype = self.ctrltype
        if ctype != '': 
            #Find the error properties
            err = (gl.data[ts.i][ctype] - setpoint)
            #for the derivative, take the last two time steps
            if ts.i >= 2:  
                derr = (gl.data[ts.i][ctype] - gl.data[ts.i-2][ctype])/(gl.data[ts.i]['t'] - gl.data[ts.i-2]['t']) 
            else:
                derr = 0
            #for the integral, last Nint average
            if ts.i >= Nint:            
                interr = sum(gl.data[ts.i-Nint : ts.i][ctype])/(Nint + 1) - setpoint
            else:
                interr = 0
            self.Me = -(pp*err + pd*derr + pint*interr)
            pl.data[ts.i]['t'] = t #store error
            pl.data[ts.i]['err'] = err
            pl.data[ts.i]['Me'] = self.Me
        else:        
            self.Me = 0

    
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
    return S

def stateDer(S,t,gl,rp,wi,tc,en,op,pl):
    '''First derivative of the state vector'''
    gl,rp,wi,tc,en,op,pl = stateSplitVec(S,gl,rp,wi,tc,en,op,pl)
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
    gamma = arctan(gl.yD/gl.xD)  # climb angle.   
    alpha = gl.theta - gamma # angle of attack
    L = (W + gl.Lalpha*alpha) * (v/gl.vb)**2 #lift
    D = gl.L/float(gl.Q)*(1 + gl.dv*(v/gl.vb-1)**2 + gl.dalpha*alpha**2) #drag
    M = (-gl.palpha*alpha - pl.Me) #torque on glider
    vgw = (gl.xD*(rp.lo - gl.x) - gl.yD*gl.y)/float(lenrope) #velocity of glider toward winch
    #winch-engine
    vrel = wi.v/en.v
    Fee = tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3     
    Few = Fee * (2-vrel)  # effective force between engine and winch through torque converter
    # Update controls
    op.control(t,gl,rp,wi,en) 
    pl.control(t,ts,gl)
    #----derivatives of state variables----#
    dotx = gl.xD    
    dotxD = 1/float(gl.m) * (rp.T*cos(thetarope) - D*cos(gamma) - L*sin(gamma)) #x acceleration
    doty = gl.yD
    if gl.y < 0.1 and L < W: #on ground 
        dotyD = 1/float(gl.m) * (L*cos(gamma) - rp.T*sin(thetarope) - D*sin(gamma)) #y acceleration on ground
    else:
        dotyD = 1/float(gl.m) * (L*cos(gamma) - rp.T*sin(thetarope) - D*sin(gamma) - W) #y acceleration
    dottheta = gl.thetaD    
    dotthetaD = 1/float(gl.Ig) * (rp.T*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/float(rp.a))-gl.theta-thetarope) + pl.Me)    
    dotT = rp.Y*rp.A*(wi.v - vgw)/float(lenrope)
    dotvw =  1/float(wi.me) * (Few - rp.T)
    dotve =  1/float(en.me) *(op.Sth * en.Pavail(en.v) / float(en.v) - Few / (2 - vrel))

    # store data from this time step for use in controls.  The ode solver enters this routine
    # usually two or more times per time step.  We advance the counter only if the time has changed 
    # by close to a nominal time step    
    if t - ts.oldt > 0.9*ts.dt: 
        ts.i += 1 
        ts.oldt = t
        gl.data[ts.i]['t']  = t
        gl.data[ts.i]['x']  = gl.x
        gl.data[ts.i]['xD'] = gl.xD
        gl.data[ts.i]['y']  = gl.y
        gl.data[ts.i]['yD'] = gl.yD
        gl.data[ts.i]['v']  = v
        gl.data[ts.i]['theta']  = gl.theta
        gl.data[ts.i]['alpha']  = alpha
        op.data[ts.i]['t']   = t
        op.data[ts.i]['Sth'] = op.Sth
    return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotT,dotvw,dotve]

##########################################################################
#                         Main script
##########################################################################                        
tStart = 0
tEnd = 20      # end time for simulation
dt = 0.05       #nominal time step, sec
ntime = (tEnd - tStart)/dt + 1   # number of time steps expected
t = linspace(tStart,tEnd,num=ntime)
ts = timesteps(tStart,tEnd,ntime)

gl = glider(ntime)
rp = rope()
wi = winch()
tc = torqconv()
en = engine(wi.rdrum)
op = operator(ntime)
pl = pilot(ntime)
pl.ctrltype = 'alpha'
  
S0 = zeros(9)
gl.xD = 1e-6  #to avoid div/zero
gl.ve = 1e-6  #to avoid div/zero
S0 = stateJoin(S0,gl,rp,wi,tc,en,op,pl)
# S = odeint(stateDer,S0,t,args=(gl,rp,wi,tc,en,op,pl,),full_output = 1)
S = odeint(stateDer,S0,t,args=(gl,rp,wi,tc,en,op,pl,))
#Split S (now a matrix with a state row for each time)

gl,rp,wi,tc,en,op,pl = stateSplitMat(S,gl,rp,wi,tc,en,op,pl)

# plot results
close("all")
plts = plots()
plts.xy(t,[x,gl.y],'time (sec)','position (m)',['x','y'],'glider position vs time')
plts.xy(t,[en.v,wi.v],'time (sec)','effective speed (m/s)',['engine','winch'],'Engine and winch speeds')
plts.xy(t,[gl.xD,gl.yD],'time (sec)','velocity (m/s)',['vx','vy'],'glider velocity vs time')
gamma = arctan(gl.yD/gl.xD) 
alpha = gl.theta - gamma  
plts.xy(t,[180/pi*gl.theta,180/pi*gamma,180/pi*alpha],'time (sec)','angle (deg)',['pitch','climb','AoA'],'flight angles')  
vrel =wi.v/en.v 
Few = (2-vrel)*tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3
plts.xy(t,[rp.T,Few],'time (sec)','Force (N)',['rope','TC-winch'],'Forces between objects')
plts.xy(op.data[:ts.i]['t'],[op.data[:ts.i]['Sth']],'time (sec)','Throttle setting',['Throttle ',' '],'Throttle')
plts.xy(pl.data[:ts.i]['t'],[pl.data[:ts.i]['err'],pl.data[:ts.i]['Me']],'time (sec)','Error (rad/s)',['error','elevator moment'],'Pilot')


print 'Done'





