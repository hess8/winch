# from IPython import get_ipython  
# get_ipython().magic('reset -sf') #comment out this line if not running in Ipython; resets all variables to zero

''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Solves for the motion of a glider launched by a winch with a springy rope. 
The winch is connected to an engine through a torque converter. 

State varables:
x, horizontal glider CG position
xD, horizontal glider CG velocity
y, vertical glider CG position
#yD, vertical glider CG velocity
theta, glider pitch angle above horizon
thetaD, glider pitch rotation rate
T, rope tension
v, rope uptake speed (m/s)
ve, effective engine speed (m/s)
Sth, throttle setting
Me, pilot controlled moment (torque) from elevator
'''
import os
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,where,amin,amax,argmin,argmax,exp
from numpy import degrees as deg
from numpy import radians as rad
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend,title,grid
from scipy.integrate import odeint
g = 9.8
 
def stateSplitMat(S,gl,rp,wi,tc,en,op,pl):
    '''Splits the formal state matrix S (each row a different time) into the various state variables'''
    gl.x      = S[:,0]
    gl.xD     = S[:,1]
    gl.y      = S[:,2]
    gl.yD     = S[:,3]
    gl.theta  = S[:,4]
    gl.thetaD = S[:,5]
    pl.elev   = S[:,6]
    rp.T      = S[:,7]
    wi.v      = S[:,8]
    en.v      = S[:,9]
    return gl,rp,wi,tc,en,op,pl
    
    
def stateSplitVec(S,gl,rp,wi,tc,en,op,pl):
    '''Splits the formal state vector S into the various state variables'''
    gl.x      = S[0]
    gl.xD     = S[1]
    gl.y      = S[2]
    gl.yD     = S[3]
    gl.theta  = S[4]
    gl.thetaD = S[5]
    pl.elev   = S[6]
    rp.T      = S[7]
    wi.v      = S[8]
    en.v      = S[9]
    return gl,rp,wi,tc,en,op,pl    
    
def stateJoin(S,gl,rp,wi,tc,en,op,pl):
    '''Joins the various state variables into the formal state vector S '''
    S[0] = gl.x
    S[1] = gl.xD
    S[2] = gl.y        
    S[3] = gl.yD       
    S[4] = gl.theta    
    S[5] = gl.thetaD 
    S[6] = pl.elev
    S[7] = rp.T
    S[8] = wi.v
    S[9] = en.v     
    return S
    
def pid(var,time,setpoint,c,j,Nint):
    '''The array c contains the [p,i,d] coefficients'''
    # Find the error properties
    err = (var[j] - setpoint)
    # for the derivative, use the last two time steps
    if j >= 2:  
        derr = (var[j] - var[j-2])/(time[j]- time[j-2]) 
    else:
        derr = 0
    #for the integral, last Nint average
    if j >= Nint:            
        interr = sum(var[j-Nint : j])/(Nint + 1) - setpoint
    else:
        interr = 0        
#    interr = sum(var[:j+1])/(j+1) - setpoint
#    print 'interr',interr
    
#    print 't,varj,err,derr,interr',time[j],var[j],err,derr,interr

    return c[0]*err + c[1]*derr + c[2]*interr

class plots:
    def __init__(self,path):
        self.i = 0 #counter for plots, so each variable has a different color
        self.path = path
        self.linew = 3.0
        self.colorsList = ['palevioletred', u'#fc4f30', u'#6d904f','darkorange', u'#8b8b8b',
        u'#348ABD', u'#e5ae38', u'#A60628', u'#7A68A6', 'mediumaquamarine', u'#D55E00', 'darkviolet',
        u'#CC79A7',  u'#0072B2', u'#30a2da',u'#009E73','peru','slateblue'] # u'#F0E442',u'#467821','slateblue'      u'#56B4E9',
        return 
    
    def xy(self,xs,ys,xlbl,ylbl,legendLabels,titlestr):
        '''To allow different time (x) arrays, we require the xs to be a list'''
        if len(xs)<len(ys): #duplicate first x to every spot in list of correct length
            xs = [xs[0] for y in ys] 
        figure(figsize=(20, 10))
        for iy,y in enumerate(ys):
            plot(xs[iy],y,color=self.colorsList[self.i],linewidth=self.linew,label=legendLabels[iy])
            self.i += 1
            if self.i > len(self.colorsList)-1:
                self.i = 0
        xlabel(xlbl)
        ylabel(ylbl)
        grid(color='k', linestyle='--', linewidth=1)
        legend(loc ='lower right',framealpha=0.5)
#        legend(loc ='upper left')
        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
        ylim([ymin,ymax])
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        show(block = False)
#         show()
        
    def xyy(self,xs,ys,yscalesMap,xlbl,ylbls,legendLabels,titlestr):
        '''Allows plotting with two yscales, left (0) and right (1).  
        So yscalesMap is e.g. [0,0,1,0] for 4 datasets'''
        fig, ax0 = subplots(figsize=(20, 10))        
        ax1 = ax0.twinx()
        if len(xs)<len(ys): #duplicate first x to every spot in list of correct length
            xs = [xs[0] for y in ys] 
        ymaxs = [[],[]]; ymins = [[],[]]       
        for iy,y in enumerate(ys):
            if yscalesMap[iy] == 0:
                ax0.plot(xs[iy],y,color=self.colorsList[self.i],linewidth=self.linew,label=legendLabels[iy])  
            else:
                ax1.plot(xs[iy],y,color=self.colorsList[self.i],linewidth=self.linew,label=legendLabels[iy])
            ymaxs[yscalesMap[iy]].append(max(y))
            ymins[yscalesMap[iy]].append(min(y))
            self.i += 1 
            if self.i > len(self.colorsList)-1:
                self.i = 0
        ax0.set_xlabel(xlbl)
        ax0.set_ylabel(ylbls[0])
        ax0.grid(color='k', linestyle='--', linewidth=1)
        ax1.grid(color='k', linestyle='--', linewidth=1)
        ax1.set_ylabel(ylbls[1])
        #Set limits: force the zeros on both sides to be the same
        min0 = min(0,min(ymins[0]))
        max0 = 1.1 * max(ymaxs[0])
        ax0.set_ylim([min0,max0])
        max1 = 1.1 * max(ymaxs[1])
        if min0 < 0:
            shift = max1 * -min0/max0
        else:
            shift = 0
        min1 = min(ymins[1]) - shift
        ax1.set_ylim([min1,max1])
        ax0.legend(loc ='lower right',framealpha=0.5)
        ax1.legend(loc ='upper right',framealpha=0.5)
#        legend(loc ='upper left')
#        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
#        ylim([ymin,ymax])
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        show(block = False)   
#         show()          
    
class timeinfo:
    def __init__(self,tStart,tEnd,N):
        '''The way to do persistent variables in python is with class variables
        
        The solver does not keep a constant time step, and evaluates the stateDer 
        function usually two or more steps per time step'''
        
        self.oldt = -1.0
        self.i = -1   #counter for time step.  Must start at -1 because we want to increment it in stateDer, not pl.control.
        self.dt = (tEnd -tStart)/float((N-1))
        self.tprobed = zeros(N*100)  #stores a t for each time the integrator enters stateDer
        self.Nprobed = 0  #number of times integrator enters stateDer
        self.tEnd = tEnd
    

class glider:
    def __init__(self,theta0,ntime):
        # parameters
        self.vb = 32              #   speed of glider at best glide angle
        self.m = 600             # kg Grob, 2 pilots vs 400 for PIK20
        self.W = self.m*9.8          #   weight (N)
        self.Q = 20             #   L/D
        self.alphas = 6*3.14/180         #   stall angle vs glider zero
        self.Co = 0.75             #   Lift coefficient {} at zero glider AoA
        self.Lalpha = 2*pi*self.W/self.Co
        self.I = 600*6/4   #   Grob glider moment of inertia, kgm^2, scaled from PIK20E
        self.ls = 4           # distance(m) between cg and stabilizer center        
        self.palpha = 1.2     #   This is from estimation and xflr: (m/rad) coefficient for air-glider pitch moment from angle of attack (includes both stabilizer,wing, fuselage)
#        self.palpha = 10      #  increased 3x !!!!  to reflect no-stall conditions observed in rotation  (m/rad) coefficient for air-glider pitch moment from angle of attack
        self.pelev = 1.2     # (m/rad) coefficient for air-glider pitch moment from elevator deflection
#        self.pelev = 6     # (m/rad) coefficient for air-glider pitch moment from elevator deflection #increased to get 45 degree climb. 
        self.maxElev = rad(30)   # (rad) maximum elevator deflection
        self.dv = 3.0            #   drag constant ()for speed varying away from vb
#        self.dalpha = 40        #   drag constant (/rad) for glider angle of attack away from zero. 
        self.de = 0.025          #   drag constant (/m) for elevator moment
        self.SsSw = 0.11         # ratio of stabilizer area to wing area (Grob)
        #logic
        self.lastvy = 0
        self.vypeaked = False
        self.state = 'onGnd'
        # state variables 
        self.x = 0
        self.xD = 0
        self.y = 0    
        self.yD = 0
        self.theta = 0  #pitch vs horizontal
        self.thetaD = 0
        self.CDCL = [1.0,0.0, 160, 1800, 5800,130000] #y = 128142x5 - 5854.9x4 - 1836.7x3 + 162.92x2 - 0.9667x + 0.9905
#         self.deltar = 0.02  #ground contact force distance of action (m)
        self.deltar = 0.05  #ground contact force distance of action (m)
        self.d_m = 0.4  # distance main wheel to CG (m) 
        self.d_t = 3.1  # distance tail wheel to CG (m) 
        self.theta0 = rad(theta0) 
        #data
        self.data = zeros(ntime,dtype = [('t', float),('x', float),('xD', float),('y', float),('yD', float),('v', float),('theta', float),\
                                    ('gamma', float),('alpha', float),('vD', float),('vgw', float),('L', float),('D', float),('L/D',float),\
                                    ('gndTorq',float),('Fmain',float),('Ftail',float),('Malpha',float),\
                                    ('Pdeliv',float),('Edeliv',float),('Emech',float)])
    def findState(self,ti):
        '''Determine where in the launch the glider is'''
        gd = self.data[ti.i]
        #One-time switches:
        if not self.vypeaked and self.yD < self.lastvy:
            self.vypeaked = True #change state     self.vypeaked = True #change state 
        #state
        if self.y < 1.0 and gd['L'] < self.W:
            self.state = 'onGnd'
        elif not self.vypeaked and self.theta < rad(30):
            self.state = 'rotate'
        elif not self.vypeaked and self.theta >= rad(30):
            self.state = 'climb' 
        elif self.vypeaked and self.theta  < rad(40):
            self.state = 'roundout'  
        self.lastvy = gl.yD 
        
#         
#         elif not self.vypeaked and self.theta < 180/pi*30:
#             self.state = 'rotate'
#         elif not self.vypeaked and self.theta >= 180/pi*30:
#             self.state = 'climb' 
#         elif self.vypeaked and self.theta  < 180/pi*40:
    
    def gndForces(self,ti,gl):
#         ymain0 = self.d_m*sin(self.theta0)
#         ytail0 = -self.d_t*sin(self.theta0)
        del_ymain = self.y + self.d_m*(sin(self.theta) - sin(self.theta0))
        del_ytail = self.y - self.d_t*(sin(self.theta) - sin(self.theta0))
        damp = 0 #wheels damp, too.  If more than this it excites the short-period oscillation
        if del_ymain < self.deltar:
            Fmain =  self.W * self.d_t/(self.d_m + self.d_t) * (1-  del_ymain/self.deltar) - damp*gl.yD
        else:
            Fmain = 0
        if del_ytail < self.deltar:
            Ftail =  self.W * self.d_m/(self.d_m + self.d_t) * (1-  (del_ytail)/self.deltar) - damp*gl.yD
        else:
            Ftail = 0
        return [Fmain, Ftail]
        
#    def gndForces(self,ti): #averages, but makes program unstable
##         ymain0 = self.d_m*sin(self.theta0)
##         ytail0 = -self.d_t*sin(self.theta0)
#        del_ymain = self.y + self.d_m*(sin(self.theta) - sin(self.theta0))
#        del_ytail = self.y - self.d_t*(sin(self.theta) - sin(self.theta0))
#        Nint = 4 #average to damp out spring oscillations
#        if del_ymain < self.deltar:
#            Fmain =  self.W * self.d_t/(self.d_m + self.d_t) * (1 - del_ymain/self.deltar)
#        else:
#            Fmain = 0
#        if del_ytail < self.deltar:
#            Ftail =  self.W * self.d_m/(self.d_m + self.d_t) * (1 - (del_ytail)/self.deltar)
#        else:
#            Ftail = 0
#        self.data[ti.i]['Fmain']  = Fmain
#        self.data[ti.i]['Ftail']  = Ftail
#        FmainAvg = sum(self.data[ti.i-Nint:ti.i+1]['Fmain'])/Nint
#        FtailAvg = sum(self.data[ti.i-Nint:ti.i+1]['Ftail'])/Nint
#        return [FmainAvg, FtailAvg]
   
class rope:
    def __init__(self,tau):
        # Rope parameters  
        self.tau = tau #sec.  Artificial damping of oscillations.          
        self.d  = 0.005     #  rope diameter (m)
        self.A = pi*self.d**2/4      #  rope area (m2)
        self.Y = 30e9             #  2400*9.8/(pi*(0.005/2)^2)/0.035  
                                 #  effective Young's modulus 30 GPa for rope from Dyneema
                                 #                 datasheet 3.5% average elongation at break,  
                                 #                 average breaking load of 2400 kg (5000 lbs)
        self.a = 0.7             #  horizontal distance (m) of rope attachment in front of CG, for Grob
        self.b = 0.1            #  vertial distance (m) of rope attachment below CG (guess that CG is at wing root height with pilots
        self.lo = 6500 * 0.305         #  6500 ft to meters initial rope length (m)
#        self.lo = 1000 
        self.Cdr = 1.0           # rope drag coefficient
        self.mu = 0.015          # rope linear mass density (kg/meter)

#        self.lo = 1000         #  initial rope length (m)
        # state variables 
        self.T = 0
        # data
        self.data = zeros(ntime,dtype = [('T', float),('torq', float),('theta',float),('Pdeliv',float),('Edeliv',float)]) 
        
    def avgT(self,ti):          
        tint = 4*self.tau         
        Nint = min(0,ti.i,ceil(tint/ti.dt))
        if ti.i >= 0:
            print ti.i,sum(self.data[ti.i-Nint:ti.i+1]['T'])/(Nint + 1)
            return sum(self.data[ti.i-Nint:ti.i+1]['T'])/(Nint + 1)
            
        else:
            return 0
            
    def Tglider(self,thetarope):
        g = 9.8  
        return self.T + self.mu * g * sin(thetarope)
        
    def thetaRopeGlider(self,ti,thetarope,vtrans,lenrope):
        rho = 1.22          # air density (kg/m^3)
        g = 9.8            
        if self.T > 1 and thetarope > rad(1):  # 1 degree     
            dragCorr = self.Cdr * rho * vtrans**2 * self.d * lenrope / self.T / 8           
            weightCorr = 0.5 * self.mu * g *  lenrope * cos(thetarope) / (self.T + self.mu * g * sin(thetarope))
            return thetarope + dragCorr + weightCorr
        else:
            return thetarope
            
class winch:
    def __init__(self):
        # Winch parameters  
        self.me = 40              #  Winch effective mass, kg
        self.rdrum = 0.4          # Drum radius (m), needed only for conversion of rpm to engine peak effective speed
        #state variables
        self.v = 0            # Rope uptake speed
        #data
        self.data = zeros(ntime,dtype = [('Pdeliv',float),('Edeliv',float)])  #energy delivered to winch by TC

class torqconv:
     def __init__(self):
         # TC parameters  
#         self.Ko = 13             #TC capacity (rad/sec/(Nm)^1/2)  K = 142 rpm/sqrt(ft.lb) = 142 rpm/sqrt(ft.lb) * 2pi/60 rad/s/rpm * sqrt(0.74 ftlb/Nm) = 12.8 (vs 12 in current model!)
         self.Ko = 2.0          
         self.lowTorq = 2.0
         self.dw = 0.13
         #data
         self.data = zeros(ntime,dtype = [('Pdeliv',float),('Edeliv',float)])  #energy delivered to TC by engine (impeller) rotation
     def invK(self,vrel):
         return 1/float(self.Ko) * tanh((1-vrel)/self.dw)

class engine:
    def __init__(self,tcUsed,rdrum):
        # Engine parameters  
        self.tcUsed = tcUsed  #model TC, or bypass it (poor description of torque and energy loss)
        self.hp = 390             # engine rated horsepower
        self.Pmax = 0.95*750*self.hp        # engine watts.  0.95 is for other transmission losses besides TC
#        self.rpmpeak = 6000       # rpm for peak power
#        self.vpeak = self.rpmpeak*2*pi/60*rdrum   #engine effectivespeed for peak power 

#        self.vpeak = 20   #  Gear 2: m/s engine effectivespeed for peak power.  This is determined by gearing, not the pure engine rpms:  

        self.vpeak = 33   #  Gear 2: m/s engine effectivespeed for peak power.  This is determined by gearing, not the pure engine rpms:  
                          # 4500rpm /5.5 (gear and differential) = 820 rmp, x 1rad/sec/10rpm x 0.4m = 33m/s peak engine speed.
#        self.vpeak = 49   #  Gear 3:  m/s engine effectivespeed for peak power.  This is determined by gearing, not the pure engine rpms:  
                          # 4500rpm /3.7 (differential) = 1200 rmp, x 1rad/sec/10rpm x 0.4m = 49 m/s peak engine speed.

        self.me = 10.0            #  Engine effective mass (kg), effectively rotating at rdrum
        self.deltaEng = 1         #  time delay (sec) of engine power response to change in engine speed
        self.pe1 = 1.0; self.pe2 = 1.0; self.pe3 = 1.0 #  engine power curve parameters, gas engine
        self.vrelMax = (self.pe2 + sqrt(self.pe2**2 + 4*self.pe1*self.pe3))/2/self.pe3  #speed ratio where Pavail goes to zero.  for the pe's = 1, this is 1.62  
        # state variables 
        self.v = 0            #engine effective speed (m/s)
        self.Few = 0          # effective force between engine and winch (could go in either engine or winch or in its own class)
         #data
        self.data = zeros(ntime,dtype = [('v',float),('Pdeliv',float),('Edeliv',float)]) #energy delivered to engine rotating mass by pistons
    
    def Pavail(self,ve):            # power curve
        vr = ve/float(self.vpeak)
#        if vr > 1:
        pcurveEval = self.pe1 * vr + self.pe2 * (vr)**2 - self.pe3 * (vr)**3
#        else:
#            p = 1.5
#            pcurveEval = 1 - (1-vr)**p
        return self.Pmax*pcurveEval
        
class operator:
    def __init__(self,targetTmax,thrmax,tRampUp,ntime):
        self.targetTmax = targetTmax
        self.thrmax = thrmax        
        self.tRampUp = tRampUp
        self.Sth = 0
        self.data = zeros(ntime,dtype = [('t', float),('Sth', float)])
        self.angleMax = rad(80) #throttle goes to zero at this rope angle
        #logical
        self.storedState = 'onGnd'
        self.oldvTarget = 0
        self.currvTarget = None
        self.tSwitch = 0
    
    def preSet(self,t):
        ### Ramp up, hold, then decrease to steady value
        steadyThr = 0.5        
        tRampUp = self.tRampUp
        tHold = .5
        tDown = tRampUp + tHold
#        tRampDown = ti.tEnd - tRampUp - tHold
        tRampDown1 = 1 #sec...transition to steady
        tRampDown2 = 120 #longer ramp down
        if t <= tRampUp:
            self.Sth =  self.thrmax/float(tRampUp) * t
        elif tRampUp < t < tDown:
            self.Sth =  self.thrmax
        elif t >= tDown:
            self.Sth = max(steadyThr*(1-(t-tDown)/float(tRampDown2)),self.thrmax * (1-(t-tDown)/float(tRampDown1)))

    def linearDown(self,t):
        tRampUp = self.tRampUp
        tHold = 0
        area = self.thrmax * 5  #Fixed area in seconds = 1/2 * thrmax *(trampUp + trampDown)
        tDown = tRampUp + tHold
#        tRampDown = ti.tEnd - tRampUp - tHold
        tRampDown = 2*area - tRampUp
        if t <= tRampUp:
            self.Sth =  self.thrmax/float(tRampUp) * t
        elif tRampUp < t < tDown:
            self.Sth =  self.thrmax
        elif t >= tDown:
            self.Sth = max(0,self.thrmax * (1-(t-tDown)/float(tRampDown)))
            
    def control(self,t,ti,gl,rp,en):
        ''' The operator changes the throttle to give certain engine speeds 
        as a function of rpm history and rope angle (not climb angle).'''
        def targetEngSpeed(thetarope,en,gl,tauOp): 
            vengTargets = {'onGnd' : 1.0*en.vpeak,'rotate' : 0.9*en.vpeak, 'climb' : 0.8*en.vpeak,'roundout' : 0.7*en.vpeak}
            target = vengTargets[gl.state] 
            if gl.state != self.storedState: #have switched
                self.tSwitch = t
                self.oldvTarget = self.currvTarget
                self.storedState = gl.state
            else:
                self.currvTarget = target
            vsmooth = target +  (self.oldvTarget - target) *exp(-(t-self.tSwitch)/tauOp)
#            print 'vsmooth',vsmooth
            return vsmooth
#        angleMax = self.angleMax         
#        tRampUp = self.tRampUp 
        thetarope = rp.data[ti.i]['theta']
        tauOp = 1.0 #sec response time
        tint = tauOp 
        Nint = min(ti.i,ceil(tint/ti.dt))                  
        #throttle control
#         print 'thetarope,angleSwitch',thetarope,angleSwitch
#        if thetarope < angleSwitch:
#            pp = -0.1; pd = -.0; pint = -.2
#            c = array([pp,pd,pint]) 
#            time = self.data['t']
#            speedControl = min(self.thrmax,pid(en.data['v'],time,en.vpeak,c,ti.i,Nint))
#            if t <= tRampUp:
#                self.Sth = min(self.thrmax/float(tRampUp) * t, speedControl)
#            else:
#                self.Sth = speedControl
#        else: #angleSwitch < thetarope <= angleMax: 
#            self.Sth = max(0,self.thrmax * (1-(thetarope - angleSwitch)/(angleMax - angleSwitch)))
        pp = -0.1; pd = -.0; pint = -.2
        c = array([pp,pd,pint]) 
        time = self.data['t']
        veTarget = targetEngSpeed(thetarope,en,gl,tauOp)
        speedControl = min(self.thrmax,max(0,pid(en.data['v'],time,veTarget,c,ti.i,Nint)))
#        if t <= tRampUp:
#            self.Sth = min(self.thrmax/float(tRampUp) * t, speedControl)
#        else:
        self.Sth = speedControl
        
    def constT(self,t,ti,gl,rp,en): 
        tauOp = 1.0 #sec response time
        tint = tauOp 
        Nint = min(ti.i,ceil(tint/ti.dt)) 
        pp = -16; pd = -2; pint = -8
        c = array([pp,pd,pint]) 
        time = self.data['t']
        tRampUp = self.tRampUp
        if t <= tRampUp:
            targetT =  self.targetTmax/float(tRampUp) * t
        else:
            if gl.state == 'roundout':
                targetT = 1.0* self.targetTmax
            else: 
                targetT = self.targetTmax
        Tcontrol = min(self.thrmax,max(0,pid(rp.data['T']/gl.W,time,targetT,c,ti.i,Nint)))
        self.Sth = Tcontrol    

class pilot:
    def __init__(self,pilotType,ntime,ctrltype,setpoint):
        self.Me = 0
        self.MeTarget = 0
        self.err = 0
        self.type = pilotType
        self.ctrltype = ctrltype
        self.setpoint = setpoint
        self.currCntrl = 0
        self.tSwitch = None
        self.data = zeros(ntime,dtype = [('t', float),('err', float),('Me', float),('elev',float)])
        self.humanT = 0.5 #sec 
        #algebraic function
        self.elevTarget = 0   
        #state variable
        self.elev = 0  # elevator deflection, differs from elevTarget by about humanT
        
    def control(self,t,ti,gl):
        #all angles in routines are in radians (vs degrees on the main script input) 
        def alphaControl(time,setpoint,ti,Nint):
            if gl.state == 'onGnd' and gl.theta >= gl.theta0 - rad(1): #no elevator authority
                pp =  0; pd =   0; pint =  0
            elif gl.state == 'onGnd' and gl.theta < gl.theta0 - rad(1):
                pp = -1; pd = -16; pint = -0
            else:
                pp = -200; pd = -6; pint = -64.0
            c = array([pp,pd,pint]) * gl.I
            al = gl.data['alpha']
            return pid(al,time,setpoint,c,ti.i,Nint)   
        def vControl(time,setpoint,ti,Nint):
            if (gl.state == 'onGnd' and gl.theta >= gl.theta0 - rad(1)) or gl.data[ti.i]['v'] < 25: #no reason to control if going too slow
                pp =  0; pd =   0; pint =  0
            else: 
                pp = 8; pd = 4; pint = 8 #when speed is too high, pitch up
            c = array([pp,pd,pint])* gl.I/gl.vb
            varr = gl.data['v']
            return pid(varr,time,setpoint,c,ti.i,Nint)
        def vDDamp(time,setpoint,ti,Nint):
            pp = 0; pd = 4; pint = 0 #when speed is too high, pitch up
            c = array([pp,pd,pint])* gl.I/gl.vb
            v = gl.data['v']
            return pid(v,time,setpoint,c,ti.i,Nint)
        def limiter(x,xmax):
            if x > xmax:
                x = xmax
            elif x < -xmax:
                x = -xmax
            return x
        
        tint = 4.0 #sec
        Nint = ceil(tint/ti.dt)   
        time = gl.data['t']
        L = gl.data[ti.i]['L']
        alpha = gl.data[ti.i]['alpha']
        gamma = gl.data[ti.i]['gamma']
        crossAngle = rad(self.setpoint[2])  #climb angle to switch from one control to another
        Mdelev =  max(0.1,L * gl.pelev/(gl.Co + 2*pi*alpha)) #ratio between moment and elevator deflection 
        if '' in control:
            self.Me = 0
            self.elev = 0
        else:
            if self.currCntrl == 0 and gamma > crossAngle:  #switch to 2nd control
                self.currCntrl = 1 #switches only once
                self.tSwitch = t          
                print 'Turned on second control at {:3.1f} sec'.format(t)
            ctype = self.ctrltype[self.currCntrl]
            setpoint = self.setpoint[self.currCntrl]
            # determine the moment demanded by the control            
            if ctype == 'vDdamp': # speed derivative control only (damps phugoid) 
                self.MeTarget = vDDamp(time,setpoint,ti,Nint)
            elif ctype == 'v': #target v with setpoint'
                self.MeTarget = vControl(time,setpoint,ti,Nint)
            elif ctype == 'alpha': # control AoA  
                self.MeTarget =  alphaControl(time,rad(setpoint),ti,Nint)    
            # implement
            if self.type =='elevControl': 
                self.elevTarget = limiter(self.MeTarget/Mdelev,gl.maxElev) # determine the elevator setting 
                self.Me = Mdelev * self.elev #update the moment from the elevator
            elif self.type =='momentControl': # bypass pilot's control of elevator and simply set the moment required, and the elevator to the corresponding angle.
                self.Me = self.MeTarget
                self.elev = limiter(self.Me/Mdelev,gl.maxElev)            
        pl.data[ti.i]['t'] = t  
        pl.data[ti.i]['Me'] = self.Me  
        pl.data[ti.i]['elev'] = self.elev           

def stateDer(S,t,gl,rp,wi,tc,en,op,pl):
    '''First derivative of the state vector'''
    ti.Nprobed +=1
    ti.tprobed[ti.Nprobed-1] = t
    gl,rp,wi,tc,en,op,pl = stateSplitVec(S,gl,rp,wi,tc,en,op,pl)
#    print 't,e,eset,M ',t,pl.elev,pl.elevTarget,pl.Me
    if gl.xD < 1e-6:
        gl.xD = 1e-6 #to handle v = 0 initial
    if en.v < 1e-6:
        en.v = 1e-6 #to handle v = 0 initial

    #----algebraic functions----#
    #rope
    thetarope = arctan(gl.y/float(rp.lo-gl.x));
    if thetarope <-1e-6: thetarope += pi #to handle overflight of winch 
    lenrope = sqrt((rp.lo-gl.x)**2 + gl.y**2)
    #glider
    v = sqrt(gl.xD**2 + gl.yD**2) # speed
    vgw = (gl.xD*(rp.lo - gl.x) - gl.yD*gl.y)/float(lenrope) #velocity of glider toward winch
    vtrans = sqrt(v**2 - vgw**2 + 1e-6) # velocity of glider perpendicular to straight line rope
    thetaRG = rp.thetaRopeGlider(ti,thetarope,vtrans,lenrope) # rope angle at glider corrected for rope weight and drag
    Tg = rp.Tglider(thetarope) #tension at glider corrected for rope weight
    if gl.xD > 1: #avoid initial zeros problem
        gamma = arctan(gl.yD/gl.xD)  # climb angle.  
    else:
        gamma = 0
    alpha = gl.theta - gamma # angle of attack
    L = (gl.W + gl.Lalpha*alpha) * (v/gl.vb)**2 #lift       
    D = L/float(gl.Q)*(1 + gl.CDCL[2]*alpha**2+gl.CDCL[3]*alpha**3+gl.CDCL[4]*alpha**4+gl.CDCL[5]*alpha**5)# + gl.de*pl.Me #drag  
    if alpha > gl.alphas: #stall mimic
        L = 0.75*L
        D = 4*L/float(gl.Q)
    alphatorq = -L * gl.palpha * alpha/(gl.Co + 2*pi*alpha)
    [Fmain, Ftail] = gl.gndForces(ti,gl)
    gndTorq = Fmain*gl.d_m - Ftail*gl.d_t
    M = alphatorq + pl.Me + gndTorq  #torque of air and ground on glider
#     ropetorq = Tg*(rp.b - rp.a * tan(gl.theta + thetaRG)) #torque of rope on glider
#     ropetorq = Tg*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/float(rp.a))-gl.theta-thetaRG) #torque of rope on glider
    ropetorq = Tg*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/float(rp.a))-gl.data[ti.i-1]['theta']-thetaRG) #torque of rope on glider

#     if t>0.126362685837:
#         print 'ropetorq',ropetorq
    #winch-engine
    vrel = wi.v/en.v
    Fee = tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3     
    Few = Fee * (tc.lowTorq-vrel)  # effective force between engine and winch through torque converter
    #----derivatives of state variables----#
    dotx = gl.xD       
    dotxD = 1/float(gl.m) * (Tg*cos(thetaRG) - D*cos(gamma) - L*sin(gamma)) #x acceleration
    doty = gl.yD
    dotyD = 1/float(gl.m) * (L*cos(gamma) - Tg*sin(thetaRG) - D*sin(gamma) - gl.W  + Fmain + Ftail) #y acceleration
    dottheta = gl.thetaD    
    dotthetaD = 1/float(gl.I) * (ropetorq + M)
    if pl.type == 'elevControl':
        dotelev = 1/pl.humanT * (pl.elevTarget-pl.elev)
    else:
        dotelev = 0 

    dotT = rp.Y*rp.A*(wi.v - vgw)/float(lenrope) #- (rp.T-rp.avgT(ti))/rp.tau
    if en.tcUsed:
        dotvw =  1/float(wi.me) * (Few - rp.T)
        dotve =  1/float(en.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - Few / (2 - vrel))
    else: #no torque converter
        dotvw = 1/float(en.me + wi.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - rp.T)
        dotve = dotvw
    # The ode solver enters this routine
    # usually two or more times per time step.  We advance the time step counter only if the time has changed 
    # by close to a nominal time step    
    if t - ti.oldt > 2.0*ti.dt: 
        ti.i += 1 
#        print 't:{:8.3f} x:{:8.3f} xD:{:8.3f} y:{:8.3f} yD:{:8.3f} T:{:8.3f} L:{:8.3f} state {}'.format(t,gl.x,gl.xD,gl.y,gl.yD,rp.T,L,gl.state)
#        print t, 't:{:8.3f} x:{:8.3f} xD:{:8.3f} y:{:8.3f} yD:{:8.3f} D/L:{:8.3f}, L/D :{:8.3f}'.format(t,gl.x,gl.xD,gl.y,gl.yD,D/L,L/D)
#        print 't,Fmain+Ftail,gndTorq,atorq,rtorq,theta',t,Fmain+Ftail,gndTorq,alphatorq,ropetorq,gl.theta
#         if rp.T > 10:
#             print 'pause'
        # store data from this time step for use /in controls or plotting.  
        gl.data[ti.i]['t']  = t
        gl.data[ti.i]['x']  = gl.x
        gl.data[ti.i]['xD'] = gl.xD
        gl.data[ti.i]['y']  = gl.y
        gl.data[ti.i]['yD'] = gl.yD
        gl.data[ti.i]['v']  = v
        gl.data[ti.i]['theta']  = gl.theta
        gl.data[ti.i]['gamma']  = gamma
        gl.data[ti.i]['alpha']  = alpha
        gl.data[ti.i]['vD']  = sqrt(dotxD**2 + dotyD**2)
        gl.data[ti.i]['vgw']  = vgw
        gl.data[ti.i]['L']  = L
        gl.data[ti.i]['D']  = D
        if D > 10:
            gl.data[ti.i]['L/D']  = L/D
        else:
            gl.data[ti.i]['L/D']  = gl.Q
        gl.data[ti.i]['gndTorq']  = gndTorq
        gl.data[ti.i]['Malpha']  = alphatorq
        gl.data[ti.i]['Emech'] = 0.5*(gl.m * v**2 + gl.I * gl.thetaD**2) + gl.m  * g * gl.y  #only glider energy here
#        gl.data[ti.i]['Pdeliv'] = rp.T * vgw 
        gl.data[ti.i]['Pdeliv'] = Tg * v * cos(thetaRG + gl.theta) 
        gl.data[ti.i]['Edeliv'] = gl.data[ti.i - 1]['Edeliv'] + gl.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        rp.data[ti.i]['Pdeliv'] = rp.T * wi.v #+ 0.5*lenrope/rp.Y/rp.A * rp.T * dotT
        rp.data[ti.i]['Edeliv'] = rp.data[ti.i - 1]['Edeliv'] + rp.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        rp.data[ti.i]['T'] = rp.T
        rp.data[ti.i]['torq'] = ropetorq
#         if gl.y > 1:
#             print 'pause'
        rp.data[ti.i]['theta'] = thetarope
        wi.data[ti.i]['Pdeliv'] = Few * wi.v
        wi.data[ti.i]['Edeliv'] = wi.data[ti.i - 1]['Edeliv'] + wi.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        en.data[ti.i]['v'] = en.v  
        en.data[ti.i]['Pdeliv'] = op.Sth * en.Pavail(en.v)         
        en.data[ti.i]['Edeliv'] = en.data[ti.i - 1]['Edeliv'] + en.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        op.data[ti.i]['t']   = t
        op.data[ti.i]['Sth'] = op.Sth
        ti.oldt = t
    #---update things that we don't want updated everytime ODEint enters stateDer
        gl.findState(ti)
        # Update controls
    #    op.linearDown(t) 
#        op.control(t,ti,gl,rp,en)
#        op.preSet(t)
        pl.control(t,ti,gl)
        op.constT(t,ti,gl,rp,en)        
    return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotelev,dotT,dotvw,dotve]

##########################################################################
#                         Main script
##########################################################################                        
tRampUpList = [3] #If you only want to run one value
tStart = 0
tEnd = 40 # end time for simulation
dt = 0.05 # nominal time step, sec
targetTmax = 1.0
thrmax =  1.0
path = 'D:\\Winch launch physics\\results\\Feb28 2018 constant T, tramp loop'  #for saving plots
if not os.path.exists(path): os.mkdir(path)
#path = 'D:\\Winch launch physics\\results\\aoa control Grob USA winch'  #for saving plots
#control = ['alpha','alpha']  # Use '' for none
#setpoint = [3 ,3 , 20]  # deg,speed, deg last one is climb angle to transition to final control
#control = ['alpha','vDdamp']
#control = ['alpha','v']
#setpoint = [2,33, 20]  # deg,speed, deg last one is climb angle to transition to final control
#control = ['v','v']
control = ['v','v']
setpoint = [30,30, 90]  # deg,speed, deg last one is climb angle to transition to final control
# control = ['','']
#control = ['alpha','v']
#setpoint = [4 ,33, 20]  #deg,speed, deg last one is climb angle to transition to final control
#control = ['','']
#setpoint = [0 , 0, 30]  # deg,speed, deglast one is climb angle to transition to final control

ropetau = 0.0 #oscillation damping in rope, artificial
pilotType = 'momentControl'  # simpler model bypasses elevator...just creates the moments demanded
# pilotType = 'elevControl' # includes elevator and response time, and necessary ground roll evolution of elevator
tcUsed = True   # uses the torque controller
# tcUsed = False  #delivers a torque to the winch determined by Sthr*Pmax/omega

# control =/ 'v'  # Use '' for none
# setpoint = 1.0                    # for velocity, setpoint is in terms of vbest: vb
ntime = ((tEnd - tStart)/dt + 1 ) * 64.0   # number of time steps to allow for data points saved

# Loop over parameters for study, optimization
# tRampUpList = linspace(1,8,50)
data = zeros(len(tRampUpList),dtype = [('tRampUp', float),('xRoll', float),('tRoll', float),('yfinal', float),('vmax', float),('vDmax', float),('Sthmax',float),\
                                    ('alphaMax', float),('gammaMax', float),('thetaDmax', float),('Tmax', float),('yDfinal', float),('Lmax', float)])
yminLoop = 100 #if yfinal is less than this height, the run failed, so ignore this time point
for iloop,tRampUp in enumerate(tRampUpList):
    print '\nThrottle ramp up time', tRampUp    
    theta0 = 6   # deg resting angle of glider on ground    
    # create the objects we need from classes
    t = linspace(tStart,tEnd,num=ntime)    
    ti = timeinfo(tStart,tEnd,ntime) 
    gl = glider(theta0,ntime)
    rp = rope(ropetau) 
    wi = winch()
    tc = torqconv()
    en = engine(tcUsed,wi.rdrum)
    op = operator(targetTmax,thrmax,tRampUp,ntime)
    pl = pilot(pilotType,ntime,control,setpoint)
    # nonzero initial conditions
    gl.xD = 1e-6  #to avoid div/zero
    gl.yD = 1e-6  #to avoid div/zero
    gl.ve = 1e-6  #to avoid div/zero
    gl.theta  = rad(theta0)
    #initialize state vector to zero  
    S0 = zeros(10)
    #integrate the ODEs
    S0 = stateJoin(S0,gl,rp,wi,tc,en,op,pl)
    S = odeint(stateDer,S0,t,args=(gl,rp,wi,tc,en,op,pl,))
    #Split S (now a matrix with state variables in columns and times in rows)
    gl,rp,wi,tc,en,op,pl = stateSplitMat(S,gl,rp,wi,tc,en,op,pl)
    #If yD dropped below zero, the release occured.  Remove the time steps after that. 
    negvyTrigger = -0.1    
    if max(gl.yD)> 1 and min(gl.yD) < negvyTrigger :
        negyD =where(gl.yD < negvyTrigger)[0] #glider losing altitude
        itr = negyD[0]-1 #ode solver index for release time
    #    itr = argmin(gl.yD)
        t = t[:itr] #shorten
        if min(gl.data['yD']) < 0:
            negyD =where(gl.data['yD'] < 0)[0]
            ti.i = negyD[0]-1  #data index for release time  
        else:
            ti.i = argmax(gl.data['t'])
    else:
        itr = len(t)
    #Shortened labels and arrays for results
    tData = gl.data[:ti.i]['t']
    gData = gl.data[:ti.i]
    wData = wi.data[:ti.i]
    pData = pl.data[:ti.i]
    eData = en.data[:ti.i]
    oData = op.data[:ti.i]
    rData = rp.data[:ti.i]
    #misc results 
    vrel =wi.v/(en.v + 1e-6)
    Few = (2-vrel)*tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3
    
    #ground roll
    if max(gl.y) > 0.01:
        iEndRoll = where(gl.y > 0.01)[0][0]
        xRoll = gl.x[iEndRoll]
        tRoll = t[iEndRoll]
    else:
        xRoll = 0  #didn't get off ground
        tRoll = 0
    # final values     
    yfinal = gData['y'][-1]
    yDfinal  = gData['yD'][-1]
    if yDfinal < 0.5: yDfinal = 0      
    # max values
    thetaDmax = max(gl.thetaD[:itr])
    vmax = max(gData['v'])
    vDmax = max(gData['vD']) #max acceleration
    Sthmax = max(oData['Sth'])
    Tmax =  max(rp.T[:itr])/gl.W
    alphaMax = max(gData['alpha'])
    Lmax = max(gData['L'])/gl.W
    gammaMax = max(gData['gamma'])
    
    # Comments to user
    print 'Final height reached: {:5.0f} m, {:5.0f} ft.  Fraction of rope length: {:4.1f}%'.format(yfinal,yfinal/0.305,100*yfinal/float(rp.lo))
    print 'Maximum speed: {:3.0f} m/s, maximum rotation rate: {:3.1f} deg/s'.format(max(gData['v']),deg(max(gl.thetaD[:itr])))
    print 'Maximum Tension factor: {:3.1f}'.format(Tmax)
    print 'Ground roll: {:5.0f} m, {:5.1f} sec'.format(xRoll,tRoll)
    print 'Final vy: {:5.1f} m/s'.format(yDfinal)
    if abs(t[-1] - gl.data[ti.i]['t']) > 5*dt:
        print '\nWarning...the integrator struggled with this model.'
        print '\tIf some of the plots have a time axis that is too short vs others, '
        print '\t...try making smoother controls.'
    # check whether to keep this run as loop data
    if yfinal < yminLoop and iloop >0:
        print 'Bad run, not saving data'
        data[iloop] = data[iloop-1] #write over this data point with the last one. 
    else:
        # Store loop results
        data[iloop]['tRampUp'] = tRampUp
        data[iloop]['xRoll'] = xRoll
        data[iloop]['tRoll'] = tRoll
        data[iloop]['yfinal'] = yfinal
        data[iloop]['yDfinal'] = yDfinal
        data[iloop]['vmax'] = vmax
        data[iloop]['vDmax'] = vDmax
        data[iloop]['Lmax'] = Lmax    
        data[iloop]['Tmax'] = Tmax
        data[iloop]['Sthmax'] = Sthmax
        data[iloop]['alphaMax'] = deg(alphaMax)
        data[iloop]['gammaMax'] = deg(gammaMax)
        data[iloop]['thetaDmax'] = deg(thetaDmax)
# plot results for last one
close('all')
plts = plots(path)   
#plts.xy([tData],[gData['xD'],gData['yD'],gData['v']],'time (sec)','Velocity (m/s)',['vx','vy','v'],'Glider velocity vs time') 
#plts.xy([t],[deg(gl.theta[:itr],deg(gamma,deg(alpha],'time (sec)','angle (deg)',['pitch','climb','AoA'],'flight angles')  
#plts.xy([t],[en.v[:itr],wi.v[:itr]],'time (sec)','effective speed (m/s)',['engine','winch'],'Engine and winch speeds')
#plts.xy([tData],[gData['L']/gl.W,gData['D']/gl.W],'time (sec)','Force/W',['lift','drag'],'Aerodynamic forces')
#plts.xy([tData],],'time (sec)','Torque (Nm) ',['rope','alpha'],'Torques on glider')
#plts.xy([t],[rp.T[:itr]/gl.W,Few[:itr]/gl.W],'time (sec)','Force/weight',['tension','TC-winch force'],'Forces between objects')
#plts.xy([tData],[oData['Sth']],'time (sec)','Throttle setting',['Throttle ',' '],'Throttle')

#glider position vs time
plts.xy([t],[gl.x[:itr],gl.y[:itr]],'time (sec)','position (m)',['x','y'],'Glider position vs time')
plts.xy([gl.x[:itr]],[gl.y[:itr]],'x (m)','y (m)',['x','y'],'Glider y vs x')
#glider speed and angles
#plts.xy([tData,tData,tData,tData,t,t,t,t,tData],[gData['xD'],gData['yD'],gData['v'],deg(gData['alpha'],deg(gl.theta[:itr],deg(gamma,deg(pl.elev[:itr],deg(gl.thetaD[:itr],gData['L/D']],\
#        'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','angle of attack','pitch','climb','elevator','pitch rate (deg/sec)','L/D'],'Glider velocities and angles')
plts.xy([tData,tData,tData,tData,t,tData,tData,t],[gData['xD'],gData['yD'],gData['v'],deg(gData['alpha']),deg(gl.theta[:itr]),deg(gData['gamma']),deg(pData['elev']),deg(gl.thetaD[:itr])],\
        'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','angle of attack','pitch','climb','elevator','pitch rate (deg/sec)'],'Glider velocities and angles')

plts.i = 0 #restart color cycle
plts.xyy([tData,t,tData,t,tData,tData,tData,t],[gData['v'],wi.v[:itr],gData['y']/rp.lo,rp.T[:itr]/gl.W,gData['L']/gl.W,deg(gData['alpha']),deg(gData['gamma']),deg(gl.thetaD[:itr])],\
        [0,0,1,1,1,0,0,0],'time (sec)',['Velocity (m/s), Angles (deg)','Relative forces and height'],['v (glider)',r'$v_r$ (rope)','height/'+ r'$\l_o $','T/W', 'L/W', 'angle of attack','climb angle','rot. rate (deg/sec)'],'Glider and rope')
#lift,drag,forces
plts.xy([tData,tData,t,t],[gData['L']/gl.W,gData['D']/gl.W,rp.T[:itr]/gl.W,Few[:itr]/gl.W],\
        'time (sec)','Forces/Weight',['lift','drag','tension','TC-winch'],'Forces')
#torques
plts.xy([tData],[rData['torq'],gData['Malpha'],pData['Me'],gData['gndTorq']],'time (sec)','Torque (Nm)',['rope','stablizer','elevator','ground'],'Torques')
#Engine, rope and winch
plts.xy([t,t,tData,tData,tData],[en.v[:itr],wi.v[:itr],gData['vgw'],deg(rData['theta']),100*oData['Sth']],'time (sec)','Speeds (effective: m/s), Angle (deg), Throttle %',['engine speed','rope speed','glider radial speed','rope angle','throttle'],'Engine and rope')        
#Energy,Power
plts.xy([tData],[eData['Edeliv']/1e6,wData['Edeliv']/1e6,rData['Edeliv']/1e6,gData['Edeliv']/1e6,gData['Emech']/1e6],'time (sec)','Energy (MJ)',['to engine','to winch','to rope','to glider','in glider'],'Energy delivered and kept')        
plts.xy([tData],[eData['Pdeliv']/en.Pmax,wData['Pdeliv']/en.Pmax,rData['Pdeliv']/en.Pmax,gData['Pdeliv']/en.Pmax],'time (sec)','Power/Pmax',['to engine','to winch','to rope','to glider'],'Power delivered')        
# plot loop results
heightLoss = data['yfinal'] - max(data['yfinal'])#vs maximum in loop
plts.i = 0 #restart color cycle
plts.xyy([data['tRampUp']],[data['xRoll'],10*data['tRoll'],data['yfinal']/rp.lo*100,heightLoss,data['vmax'],data['vDmax']/g,data['Sthmax'],data['Tmax'],data['Lmax'],data['alphaMax'],data['gammaMax'],data['thetaDmax']],\
        [0,0,0,0,0,1,1,1,1,0,0,0],'throttle ramp-up time (sec)',['Velocity (m/s), Angles (deg), m, sec %',"Relative forces,g's"],\
        ['x gnd roll', 't gnd roll x 10','height/'+ r'$\l_o $%','Height diff',r'$v_{max}$',"max g's",'max throttle',r'$T_{max}/W$', r'$L_{max}/W$', r'$\alpha_{max}$',r'$\gamma_{max}$','rot. max (deg/sec)'],'Flight results vs throttle ramp-up time')
print 'Done'