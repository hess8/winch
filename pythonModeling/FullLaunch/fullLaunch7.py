# from IPython import get_ipython  
# get_ipython().magic('reset -sf') #comment out this line if not running in Ipython; resets all variables to zero

''' Bret Hess, bret.hess@gmail.com 
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
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,where,\
    amin,amax,argmin,argmax,exp,mean 
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
    return c[0]*err + c[1]*derr + c[2]*interr    
    
def smooth(data,time,N):
    '''Smooths data with a running average over tsmooth. data and time are arrays
    points are not evenly spaced in time. Smoothing is centered on t.
    Use a triangular weighting'''
    tsmooth = float(1.5) #sec (float is a precaution: it must not be an integer time)  
#    dt = 0.1 #sec
    smoothed = data
    tfinal = time[-1]
    for ism in range(N): #smooth any number of times.
        for i,t in enumerate(time):
             if tsmooth/2  < t < tfinal - tsmooth/2:
                 dsum = 0
                 totweight = 0
                 iearly = where( time < t - tsmooth/2)[0][-1] 
                 ilater = where( time > t + tsmooth/2)[0][0]                 
                 for it in range(iearly,ilater+1): 
                     weight = (time[it]-time[it-1])*(1-abs(time[it]-t)/(tsmooth/2))
                     dsum += smoothed[it] * weight
                     totweight += weight            
                 smoothed[i] = dsum/totweight
    return smoothed
        
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
        #data
        self.data = zeros(ntime,dtype = [('t', float)])
    

class glider:
    def __init__(self,theta0,ntime):
        # parameters
        self.vb = 32              #   speed of glider at best glide angle
        self.m = 600             # kg Grob, 2 pilots vs 400 for PIK20
        self.W = self.m*9.8          #   weight (N)
        self.Q = 30             #   L/D
        self.alphas = rad(6.8)         #  stall angle vs glider zero
        self.Co = 0.65             #   Lift coefficient {} at zero glider AoA
        self.CLalpha = 5.9        
        self.Lalpha = self.CLalpha*self.W/self.Co #from xflr.  A little less than 2 pi because part of the lift if from the elevator.  
        self.I = 600*6/4   #   Grob glider moment of inertia, kgm^2, scaled from PIK20E
        self.ls = 4           # distance(m) between cg and stabilizer center        
        self.palpha = 1.9     #   This is from xflr model: (m/rad) coefficient for air-glider pitch moment from angle of attack (includes both stabilizer,wing, fuselage)
        self.pelev =  1.2     # (m/rad) coefficient for air-glider pitch moment from elevator deflection
        self.maxElev = rad(20)   # (rad) maximum elevator deflection
        self.dv = 3.0            #   drag constant ()for speed varying away from vb
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
        self.deltar = 0.02  #ground contact force distance of action (m)
#        self.deltar = 0.05  #ground contact force distance of action (m)
        self.d_m = 0.25  # distance main wheel to CG (m) 
        self.d_t = 3.3  # distance tail wheel to CG (m) 
        self.theta0 = rad(theta0) 
        #data
        self.data = zeros(ntime,dtype = [('x', float),('xD', float),('y', float),('yD', float),('v', float),('theta', float),('thetaD', float),\
                                    ('gamma', float),('alpha', float),('vD', float),('vgw', float),('L', float),('D', float),('L/D',float),\
                                    ('gndTorq',float),('Fmain',float),('Ftail',float),('Malpha',float),\
                                    ('Pdeliv',float),('Edeliv',float),('Emech',float)])
    def findState(self,t,ti,rp):
        '''Determine where in the launch the glider is'''
        gd = self.data[ti.i]
        #One-time switches:
        if not self.vypeaked  and not self.state == 'onGnd' and self.yD < self.lastvy:
            self.vypeaked = True 
        #state
        if self.y < 1.0 and gd['L'] < self.W:
            self.state = 'onGnd'
        elif not self.vypeaked and self.theta < rad(10):
            self.state = 'rotate'
        elif not self.vypeaked and self.theta >= rad(10):
            self.state ='initClimb' 
        elif self.vypeaked and self.theta  < rad(40) and rp.data[ti.i]['theta'] < rp.thetaCutThr:
            self.state = 'mainClimb'  
        elif rp.data[ti.i]['theta'] >= rp.thetaCutThr :
            self.state = 'prepRelease'
        self.lastvy = gl.yD 
    
    def gndForces(self,ti,gl,rp):
        del_ymain = self.y + self.d_m*(sin(self.theta) - sin(self.theta0))
        del_ytail = self.y - self.d_t*(sin(self.theta) - sin(self.theta0))
        muS = 0.1  # static friction, including runway bumpiness
        muK = 0.05 # kinetic friction, including runway bumpiness
        damp = 0 #wheels damp, too.  If more than this it excites the short-period oscillation
        if del_ymain < self.deltar:
            Fmain =  self.W * self.d_t/(self.d_m + self.d_t) * (1-  del_ymain/self.deltar) - damp*gl.yD
        else:
            Fmain = 0
        if del_ytail < self.deltar:
            Ftail =  self.W * self.d_m/(self.d_m + self.d_t) * (1-  (del_ytail)/self.deltar) - damp*gl.yD
        else:
            Ftail = 0
        if gl.state == 'onGnd' :
            if gl.xD < 0.001 and rp.T < muS*gl.W: #not moving             
                Ffric = rp.T
            else: #static broken
                Ffric =  muK*gl.W
        else:
            Ffric = 0
        return [Fmain, Ftail,Ffric]
   
class rope:
    def __init__(self,thetaMax):
        # Rope parameters  
        self.thetaMax =rad(thetaMax) # release angle of rope         
        self.thetaCutThr =  0.8*self.thetaMax
        self.d  = 0.005     #  rope diameter (m)
        self.A = pi*self.d**2/4      #  rope area (m2)
        self.Ys = 30e9             #  2400*9.8/(pi*(0.005/2)^2)/0.035  
                                 #  effective static Young's modulus 30 GPa for rope from Dyneema
                                 #                 datasheet 3.5% average elongation at break,  
                                 #                 average breaking load of 2400 kg (5000 lbs)
        self.hystRate = 100     # hysteresis rate (1/sec) for dynamic stiffness. To turn this effect off, make this rate large, not small
        self.a = 0.7             #  horizontal distance (m) of rope attachment in front of CG, for Grob
        self.b = 0.3            #  vertial distance (m) of rope attachment below CG (guess that CG is at wing root height with pilots
        self.lo = 6500 * 0.305         #  6500 ft to meters initial rope length (m)
#        self.lo = 1000 
        self.Cdr = 1.0           # rope drag coefficient
        self.mu = 0.015          # rope linear mass density (kg/meter)
        # state variables 
        self.T = 0
        # data
        self.data = zeros(ntime,dtype = [('T', float),('torq', float),('theta',float),('Pdeliv',float),('Edeliv',float)]) 
        
        
    def avgT(self,ti):          
        tint = 4*self.tau         
        Nint = min(0,ti.i,ceil(tint/ti.dt))
        if ti.i >= 0:
            print(ti.i,sum(self.data[ti.i-Nint:ti.i+1]['T'])/(Nint + 1))
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
                    
    def chgT(self,vrad,vgw,lenrope):
        '''RHS of diff equation for tension change in terms of rope and glider speeds.  Allows hysteresis'''         
        p = 1    
#        print('hyst', 1/(1- (abs(vrad - vgw)/self.hystRate/lenrope)**p))
        return rp.Ys*rp.A*(vrad - vgw)/lenrope/(1- (abs(vrad - vgw)/self.hystRate/lenrope)**p)       
            
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
         self.Ko = 1.0          
         self.lowTorqR = 2.0
         self.dw = 0.13
         #data
         self.data = zeros(ntime,dtype = [('Pdeliv',float),('Edeliv',float)])  #energy delivered to TC by engine (impeller) rotation
     def invK(self,vrel):
         return 1/float(self.Ko) * tanh((1-vrel)/self.dw)

class engine:
    '''This version models the torque curve with parameters rather than the power curve, and allows setting the peak torque and peak power 
    rpms independently'''
    def __init__(self,tcUsed,rdrum):
        # Engine parameters  
        self.gear = 1.5    # 2nd gear ratio
        self.diff = 3.7 
        self.tcUsed = tcUsed  #model TC, or bypass it (poor description of torque and energy loss)
        self.hp = 300            # engine  horsepower
        self.Pmax = 0.95*750*self.hp        # engine watts.  0.95 is for other transmission losses besides TC
        self.torqMax = 550/0.74*self.Pmax/(0.95*750*390)  #ft lbs converted to Nm.        
        self.rpmPeak = 4500
        self.omegaPeak = self.rpmPeak*2*pi/60 #peak power ening speed
        self.vLimit = 5000*2*pi/60*wi.rdrum/self.gear/self.diff  #engine speed limiter
    #differential gear ratio        
        self.vpeakP = self.omegaPeak*wi.rdrum/self.gear/self.diff
        self.me = 10.0            #  Engine effective mass (kg), effectively rotating at rdrum
        self.deltaEng = 1         #  time delay (sec) of engine power response to change in engine speed
        self.vpeakTorq = 0.5 * self.vpeakP # engine speed for peak torque available.
        vr = self.vpeakTorq/self.vpeakP
        Pr = self.Pmax/(self.torqMax*self.omegaPeak)
        self.c0 = -((1 - 3* vr + 4* Pr* vr**2 - 2* Pr* vr**3)/(-1 + vr)**3) #coefficient for torque curve
        self.c1 = -((6* vr - 8* Pr* vr + Pr* vr**2 + Pr* vr**3)/(-1 + vr)**3)
        self.c2 = -((3 - 4* Pr + 3* vr - 4* Pr* vr + 2* Pr* vr**2)/(-1 + vr)**3)    
        self.c3 = -((2 - 3* Pr + Pr* vr)/(-1 + vr)**3)
        self.idle = self.vpeakP * 0.13

        # state variables 
        self.v = 0            #engine effective speed (m/s)
        self.Few = 0          # effective force between engine and winch (could go in either engine or winch or in its own class)
         #data
        self.data = zeros(ntime,dtype = [('v',float),('Pdeliv',float),('torq',float),('Edeliv',float)]) #energy delivered to engine rotating mass by pistons
     
    def Pavail(self,ve):            # power curve
        
        # First calculate available torque from pistons on engine rotating parts (torque curve).  
        vr = ve/self.vpeakP #relative engine speed vs peak power speed
        torqAvail = self.torqMax * (self.c0 + self.c1*vr - self.c2*vr**2 + self.c3*vr**3) #in Nm
        return torqAvail * (ve/self.vpeakP*self.omegaPeak)
                
class operator:
    def __init__(self,throttleType,targetT,thrmax,tRampUp,tHold,ntime):
        self.throttleType = throttleType        
        self.targetT = targetT
        self.thrmax = thrmax        
        self.tRampUp = tRampUp
        self.Sth = 0
        self.SthOld = 0
        self.angleMax = rad(80) #throttle goes to zero at this rope angle
        #logical
        self.storedState = 'onGnd'
        self.oldvTarget = 0
        self.currvTarget = None
        self.tSwitch = 0
        self.tSlackEnd = None
        self.tHold = tHold
        self.thrSlack = 0.1
        self.vSlackEnd = 1  #m/s
        #data
        self.data = zeros(ntime,dtype = [('Sth', float)])

    def control(self,t,ti,gl,rp,en):
        tRampUp = self.tRampUp        
        if self.tSlackEnd is None and gl.xD > self.vSlackEnd:  #one-time event
            self.tSlackEnd = t
        if self.throttleType == 'constT':
            if gl.xD < self.vSlackEnd: #take out slack
                self.Sth = self.thrSlack
            else:
                tauOp = 1.0 #sec response time
                tint = tauOp 
                Nint = min(ti.i,ceil(tint/ti.dt)) 
    #            pp = -16; pd = -3; pint = -8
    #            pp = -8; pd = -2; pint = -8
    #            pp = -2; pd = -4; pint = -8
    #            pp = -8; pd = -8; pint = -32 
                if gl.state != 'prepRelease':
                    pp = -8; pd = -8; pint = -32 
                else:
                    pp = -.0; pd = -0; pint = -0 
    #            pp = -8; pd = -8; pint = -32 
                c = array([pp,pd,pint]) 
                time = ti.data['t']
                if t <= tRampUp:
                    targetT =  self.targetT/float(tRampUp) * (t - self.tSlackEnd)
                else:
                    if gl.state == 'mainClimb':
                        targetT = 1.0* self.targetT
                    else: 
                        targetT = self.targetT
                Tcontrol = min(self.thrmax,max(0,pid(rp.data['T']/gl.W,time,targetT,c,ti.i,Nint)))
                #limit the throttle change to 40%/second
                maxrate = 0.4 #40%/sec
                rate = (Tcontrol - self.data[ti.i]['Sth'])/(t-ti.data[ti.i-1]['t'])
                if en.v > en.vLimit:
                    self.Sth = 0.9 * self.Sth
                elif rate > 0:
                    self.Sth = self.data[ti.i]['Sth'] + min(maxrate,rate) * (t-ti.data[ti.i-1]['t'])
                else:
                    self.Sth = self.data[ti.i]['Sth'] + max(-maxrate,rate) * (t-ti.data[ti.i-1]['t'])
        elif self.throttleType == 'preset':
            ### Ramp up, hold, then decrease to steady value
            tSlackEnd = self.tSlackEnd             
            steadyThr = 0.75        
            tDown = tRampUp + self.tHold
    #        tRampDown = ti.tEnd - tRampUp - self.tHold
            tRampDown1 = 1 #sec...transition to steady
            tRampDown2 = 120 #longer ramp down
            if gl.xD < self.vSlackEnd:
                self.Sth = self.thrSlack
            elif  tSlackEnd  <= t <  tSlackEnd + tRampUp:
                self.Sth = self.thrSlack + (self.thrmax - self.thrSlack) * (t - tSlackEnd)/float(tRampUp)
            elif tRampUp < t  < tDown + tSlackEnd:
                self.Sth =  self.thrmax
            elif t + tSlackEnd >= tDown:
                self.Sth = max(steadyThr*(1-((t + tSlackEnd)-tDown)/float(tRampDown2)),self.thrmax * (1-((t + tSlackEnd)-tDown)/float(tRampDown1)))
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
        self.data = zeros(ntime,dtype = [('err', float),('Me', float),('elev',float)])
        self.humanT = 0.5 #sec 
        #algebraic function
        self.elevTarget = 0   
        #state variable
        self.elev = 0  # elevator deflection, differs from elevTarget by about humanT
        
    def control(self,t,ti,gl):
        #all angles in routines are in radians (vs degrees on the main script input) 
        def alphaControl(t,time,setpoint,ti,Nint):
            if gl.state == 'onGnd' and gl.theta >= gl.theta0 - rad(1): #no elevator authority
                pp =  0; pd =   0; pint =  0
            elif gl.state == 'onGnd' and gl.theta < gl.theta0 - rad(1):
                pp = -1; pd = -16; pint = -0
            else:
                pp = -200; pd = -6; pint = -64.0
            c = array([pp,pd,pint]) * gl.I
            al = gl.data['alpha']
            return pid(al,time,setpoint,c,ti.i,Nint)   
        def vControl(t,time,setpoint,ti,Nint):
            if (gl.state == 'onGnd' and gl.theta >= gl.theta0 - rad(1)) or gl.data[ti.i]['v'] < 25: #no reason to control if going too slow
                pp =  0; pd =   0; pint =  0
            elif gl.state == 'onGnd':
                pp = 0; pd = 0; pint = 0 
            elif gl.state == 'rotate': 
                pp = 16; pd = 0; pint = 0   #when speed is too high, pitch up
            elif gl.state =='initClimb': 
                pp = 16; pd = 16; pint = 0  
            elif gl.state == 'mainClimb': 
                pp = 48; pd = 0; pint = 16 #These values when switching from thetaD control
            elif gl.state == 'prepRelease':
                pp =  0; pd =   0; pint =  0
            c = array([pp,pd,pint])* gl.I/gl.vb
            varr = gl.data['v']
            return pid(varr,time,setpoint,c,ti.i,Nint)

            return pid(varr,time,setpoint,c,ti.i,Nint)
        def vDDamp(t,time,setpoint,ti,Nint):
            pp = 0; pd = 4; pint = 0 #when speed is too high, pitch up
            c = array([pp,pd,pint])* gl.I/gl.vb
            v = gl.data['v']
            return pid(v,time,setpoint,c,ti.i,Nint)
            
        def thetaDContr(t,time,setpoint,ti,Nint):
            if gl.state in ['rotate','initClimb','mainClimb'] and gl.data[ti.i]['v'] > 25: #m/s:                         
#                 pp = 3.0; pd = 8; pint = 2
                pp = 3.0; pd = 8; pint = 0
            else:
                pp = 0; pd = 0; pint = 0
            c = array([pp,pd,pint])* gl.I/gl.vb
            varr = gl.data['thetaD']
            return pid(varr,time,setpoint,c,ti.i,Nint)
            
        def limiter(x,xmax):
            if x > xmax:
                x = xmax
            elif x < -xmax:
                x = -xmax
            return x
        tint = 4.0 #sec
        Nint = ceil(tint/ti.dt)   
        time = ti.data['t']
        L = gl.data[ti.i]['L']
        alpha = gl.data[ti.i]['alpha']
        gamma = gl.data[ti.i]['gamma']
        crossAngle = rad(self.setpoint[2])  #climb angle to switch from one control to another
        Mdelev =  max(0.1,L * gl.pelev/(gl.Co + gl.CLalpha*alpha)) #ratio between moment and elevator deflection.  Avoid zero velocity case with the 0.1.  
        maxMe = gl.maxElev * Mdelev
        if '' in control:
            self.Me = 0
            self.elev = 0
        else:
            if self.currCntrl == 0 and (gl.y>10 and (gamma > crossAngle or gl.thetaD<rad(10))):  #switch to 2nd control
                self.currCntrl = 1 #switches only once
                self.tSwitch = t          
                print('Turned on second control at {:3.1f} sec'.format(t))
            ctype = self.ctrltype[self.currCntrl]
            setpoint = self.setpoint[self.currCntrl]
            # determine the moment demanded by the control            
            if ctype == 'vDdamp': # speed derivative control only (damps phugoid) 
                self.MeTarget = vDDamp(t,time,setpoint,ti,Nint)
            elif ctype == 'v': #target v with setpoint'
                self.MeTarget = vControl(t,time,setpoint,ti,Nint)
            elif ctype == 'alpha': # control AoA  
                self.MeTarget =  alphaControl(t,time,rad(setpoint),ti,Nint)  
            elif ctype == 'thetaD': #control pitch rate
                self.MeTarget =  thetaDContr(t,time,rad(setpoint),ti,Nint)  
            thetaDContr
            # implement
            if self.type =='elevControl': 
                self.elevTarget = limiter(self.MeTarget/Mdelev,gl.maxElev) # determine the elevator setting 
                self.Me = Mdelev * self.elev #update the moment from the elevator
            elif self.type =='momentControl': # bypass pilot's control of elevator and simply set the moment required, and the elevator to the corresponding angle.
                self.Me = limiter(self.MeTarget,maxMe)
                self.elev = self.Me/Mdelev
        pl.data[ti.i]['Me'] = self.Me  
        pl.data[ti.i]['elev'] = self.elev           

def stateDer(S,t,gl,rp,wi,tc,en,op,pl):
    '''First derivative of the state vector'''  
    if rp.data[ti.i]['theta'] > rp.thetaMax: # glider released but the integrator must finish   
        return zeros(len(S))
    else: 
        gl,rp,wi,tc,en,op,pl = stateSplitVec(S,gl,rp,wi,tc,en,op,pl)
    #    print('t,e,eset,M ',t,pl.elev,pl.elevTarget,pl.Me)
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
        alphatorq = -L * gl.palpha * alpha/(gl.Co + gl.CLalpha*alpha)
        [Fmain, Ftail, Ffric] = gl.gndForces(ti,gl,rp)
        gndTorq = Fmain*gl.d_m - Ftail*gl.d_t
        M = alphatorq + pl.Me + gndTorq  #torque of air and ground on glider
        ropetorq = Tg*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/float(rp.a))-gl.theta-thetaRG) #torque of rope on glider
        #winch-engine
        vrel = wi.v/en.v
        Fee = tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3     
        Few = Fee * (tc.lowTorqR-vrel)  # effective force between engine and winch through torque converter
        #----derivatives of state variables----#
        dotx = gl.xD       
        dotxD = 1/float(gl.m) * (Tg*cos(thetaRG) - D*cos(gamma) - L*sin(gamma) - Ffric) #x acceleration
        doty = gl.yD
        dotyD = 1/float(gl.m) * (L*cos(gamma) - Tg*sin(thetaRG) - D*sin(gamma) - gl.W  + Fmain + Ftail) #y acceleration
        dottheta = gl.thetaD    
        dotthetaD = 1/float(gl.I) * (ropetorq + M)
        if pl.type == 'elevControl':
            dotelev = 1/pl.humanT * (pl.elevTarget-pl.elev)
        else:
            dotelev = 0 
        dotT = rp.chgT(wi.v,vgw,lenrope)                       
        if en.tcUsed:
            dotvw =  1/float(wi.me) * (Few - rp.T)
            dotve =  1/float(en.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - Fee)
        else: #no torque converter
            dotvw = 1/float(en.me + wi.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - rp.T)
            dotve = dotvw
        # The ode solver enters this routine
        # usually two or more times per time step.  We advance the time step counter only if the time has changed 
        # by close to a nominal time step    
        if t - ti.oldt > 1.0*ti.dt: 
            ti.i += 1 
    #         if t > 15 and gl.yD<0:
#             print('t:{:8.3f} x:{:8.3f} xD:{:8.3f} y:{:8.3f} yD:{:8.3f} T:{:8.3f} L:{:8.3f} state {}'.format(t,gl.x,gl.xD,gl.y,gl.yD,rp.T,L,gl.state))
#            print('t:{:8.3f} x:{:8.3f} xD:{:8.3f} y:{:8.3f} yD:{:8.3f} T:{:8.3f} L:{:8.3f} gndTorq: {:8.3f} state {:12s} '.format(t,gl.x,gl.xD,gl.y,gl.yD,rp.T,L,gndTorq,gl.state))
    #             print('pause')
    #        print(t, 't:{:8.3f} x:{:8.3f} xD:{:8.3f} y:{:8.3f} yD:{:8.3f} D/L:{:8.3f}, L/D :{:8.3f}'.format(t,gl.x,gl.xD,gl.y,gl.yD,D/L,L/D))
#            print('t,state,y',t,gl.state,gl.y)
    #         if rp.T > 10:
    #             print('pause')
            # store data from this time step for use /in controls or plotting.  
            ti.data[ti.i]['t']  = t
            gl.data[ti.i]['x']  = gl.x
            gl.data[ti.i]['xD'] = gl.xD
            gl.data[ti.i]['y']  = gl.y
            gl.data[ti.i]['yD'] = gl.yD
            gl.data[ti.i]['v']  = v
            gl.data[ti.i]['theta']  = gl.theta
            gl.data[ti.i]['thetaD']  = gl.thetaD
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
            gl.data[ti.i]['Pdeliv'] = Tg * v * cos(thetaRG + gl.theta) 
            gl.data[ti.i]['Edeliv'] = gl.data[ti.i - 1]['Edeliv'] + gl.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            rp.data[ti.i]['Pdeliv'] = rp.T * wi.v 
            rp.data[ti.i]['Edeliv'] = rp.data[ti.i - 1]['Edeliv'] + rp.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            rp.data[ti.i]['T'] = rp.T
            rp.data[ti.i]['torq'] = ropetorq
            rp.data[ti.i]['theta'] = thetarope
            wi.data[ti.i]['Pdeliv'] = Few * wi.v
            wi.data[ti.i]['Edeliv'] = wi.data[ti.i - 1]['Edeliv'] + wi.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            en.data[ti.i]['v'] = en.v  
            en.data[ti.i]['Pdeliv'] = op.Sth * en.Pavail(en.v)  
            en.data[ti.i]['torq'] =  op.Sth * en.Pavail(en.v)/(en.v/wi.rdrum*en.gear*en.diff)  #from pistons          
            en.data[ti.i]['Edeliv'] = en.data[ti.i - 1]['Edeliv'] + en.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            op.data[ti.i]['Sth'] = op.Sth
            ti.oldt = t
        #---update things that we don't need done ODEint enters stateDer
            gl.findState(t,ti,rp)
            # Update controls
            pl.control(t,ti,gl)
            op.control(t,ti,gl,rp,en) 
            op.SthOld = op.Sth
        return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotelev,dotT,dotvw,dotve]

##########################################################################
#                         Main script
##########################################################################                        
#tRampUpList = linspace(3,10,10)
tRampUpList = [3] #If you only want to run one value
tHold = 1.5
tStart = 0
tEnd = 15 # end time for simulation
dt = 0.05 # nominal time step, sec
targetT = 1.0
thrmax =  1.0
ropeThetaMax = 75 #degrees
smoothed = False
#smoothed = True
throttleType = 'constT'
#throttleType = 'preset'
#path = 'D:\\Winch launch physics\\results\\Mar5 2018 preset controlled v'  #for saving plots
#path = 'D:\\Winch launch physics\\results\\test'  #for saving plots
path = 'D:\\Winch launch physics\\results\\test2'
if not os.path.exists(path): os.mkdir(path)
#path = 'D:\\Winch launch physics\\results\\aoa control Grob USA winch'  #for saving plots
#control = ['alpha','alpha']  # Use '' for none
#setpoint = [3 ,3 , 90]  # deg,speed, deg last one is climb angle to transition to final control
control = ['thetaD','v']  # Use '' for none
setpoint = [10 ,30 , 45]  # deg,speed, deg last one is climb angle to transition to final control
#control = ['alpha','v']
#setpoint = [4,30, 20]  # deg,speed, deg last one is climb angle to transition to final control
#control = ['v','v']
#setpoint = [30,30, 90]  # deg,speed, deg last one is climb angle to transition to final control
#control = ['','']
#setpoint = [0 , 0, 30]  # deg,speed, deglast one is climb angle to transition to final control
ropetau = 0.0 #oscillation damping in rope, artificial
pilotType = 'momentControl'  # simpler model bypasses elevator...just creates the moments demanded
# pilotType = 'elevControl' # includes elevator and response time, and necessary ground roll evolution of elevator
tcUsed = True   # uses the torque controller
# tcUsed = False  #delivers a torque to the winch determined by Sthr*Pmax/omega

# control =/ 'v'  # Use '' for none
# setpoint = 1.0                    # for velocity, setpoint is in terms of vbest: vb
ntime = ((tEnd - tStart)/dt + 1 )   # number of time steps to allow for data points saved
# Loop over parameters for study, optimization

data = zeros(len(tRampUpList),dtype = [('tRampUp', float),('xRoll', float),('tRoll', float),('yfinal', float),('vmax', float),('vDmax', float),('Sthmax',float),\
                                    ('alphaMax', float),('gammaMax', float),('thetaDmax', float),('Tmax', float),('Tavg', float),('yDfinal', float),('Lmax', float)])
yminLoop = 100 #if yfinal is less than this height, the run failed, so ignore this time point
for iloop,tRampUp in enumerate(tRampUpList):
    theta0 = 6   # deg resting angle of glider on ground    
    # create the objects we need from classes
    t = linspace(tStart,tEnd,num=ntime)    
    ti = timeinfo(tStart,tEnd,ntime) 
    gl = glider(theta0,ntime)
    rp = rope(ropeThetaMax) 
    wi = winch()
    tc = torqconv()
    en = engine(tcUsed,wi.rdrum)
    op = operator(throttleType,targetT,thrmax,tRampUp,tHold,ntime)
    pl = pilot(pilotType,ntime,control,setpoint)
    # nonzero initial conditions
    gl.xD = 1e-6  #to avoid div/zero
    gl.yD = 1e-6  #to avoid div/zero
    gl.v = 1e-6  #to avoid div/zero
    en.v = en.idle
    gl.theta  = rad(theta0)
    #initialize state vector to zero  
    S0 = zeros(10)
    #integrate the ODEs
    S0 = stateJoin(S0,gl,rp,wi,tc,en,op,pl)
    S = odeint(stateDer,S0,t,args=(gl,rp,wi,tc,en,op,pl))
    #Split S (now a matrix with state variables in columns and times in rows)
    gl,rp,wi,tc,en,op,pl = stateSplitMat(S,gl,rp,wi,tc,en,op,pl)
    #Find where release occurred in data.  Remove the time steps after that. 
    if max(rp.data['theta'])  > rad(ropeThetaMax):
        ti.i = where(rp.data['theta'] > rad(ropeThetaMax))[0][0]
        tfinal = ti.data[ti.i]['t']
        itr = where(t >= tfinal)[0][0]
    else:
        ti.i = argmax(ti.data['t'])
        itr = len(t)
    #Find where release occurred in state data
    
    t = t[:itr] #shorten
    
# where(gl.yD < negvyTrigger/2)[0]    
#     if max(gl.yD)> 1 and min(gl.yD) < negvyTrigger/2 :
#         negyD =where(gl.yD < negvyTrigger/2)[0] #glider losing altitude
#         itr = negyD[0]-1 #ode solver index for release time
#     #    itr = argmin(gl.yD)
#         t = t[:itr] #shorten
#         if min(gl.data['yD']) < negvyTrigger/2 :
#             negyD =where(gl.data['yD'] < negvyTrigger/2 )[0]
#             ti.i = negyD[0]-1  #data index for release time  
#         else:
#             ti.i = argmax(ti.data['t'])
#     else:
        
    #Shortened labels and arrays for results
    tData = ti.data[:ti.i]['t']
    gData = gl.data[:ti.i]
    wData = wi.data[:ti.i]
    pData = pl.data[:ti.i]
    eData = en.data[:ti.i]
    oData = op.data[:ti.i]
    rData = rp.data[:ti.i]

    if smoothed:
    #define smoothed data arrays before plotting
        print('\nSmoothing data')
        x = smooth(gl.x[:itr],t,1)
        y = smooth(gl.y[:itr],t,1)
        xD = smooth(gl.xD[:itr],t,1)
        yD = smooth(gl.yD[:itr],t,1)
        v = smooth(gData['v'],tData,1)
        vD = smooth(gData['vD'],tData,1) 
        alpha = smooth(gData['alpha'],tData,3)
        theta = smooth(gl.theta[:itr],t,2)
        gamma = smooth(gData['gamma'],tData,1)
        elev = smooth(pData['elev'],tData,1)
        thetaD= smooth(gl.thetaD[:itr],t,3)
        wiv = smooth(wi.v[:itr],t,3)
        env = smooth(en.v[:itr],t,3)
        L = smooth(gData['L'],tData,3)
        D = smooth(gData['D'],tData,2)
        T = smooth(rp.T[:itr],t,3)
        vgw = smooth(gData['vgw'],tData,1)
        Malpha = smooth(gData['Malpha'],tData,1)
        Me = smooth(pData['Me'],tData,1)
        engP = smooth(eData['Pdeliv'],tData,1)
        engTorq = smooth(eData['torq'],tData,1)
        Sth = smooth(oData['Sth'],tData,3)
        winP = smooth(wData['Pdeliv'],tData,1)
        ropP = smooth(rData['Pdeliv'],tData,1)
        gliP = smooth(gData['Pdeliv'],tData,1)
        gndTorq = smooth(gData['gndTorq'],tData,1)
        ropeTheta = smooth(rData['theta'],tData,1)
        ropeTorq = smooth(rData['torq'],tData,1)
    else:
        #Shorten labels before plotting
        x = gl.x[:itr]
        y = gl.y[:itr]
        xD = gl.xD[:itr]
        yD = gl.yD[:itr]
        v = gData['v']
        vD = gData['vD'] 
        alpha = gData['alpha']
        theta = gl.theta[:itr]
        gamma = gData['gamma']
        elev = pData['elev']
        thetaD = gl.thetaD[:itr]
        wiv = wi.v[:itr]
        env = en.v[:itr]
        L = gData['L']
        D = gData['D']
        T = rp.T[:itr]
        vgw = gData['vgw']
        Malpha = gData['Malpha']
        Me = pData['Me']
        engP = eData['Pdeliv']
        engTorq = eData['torq']
        Sth = oData['Sth']
        winP = wData['Pdeliv']
        ropP = rData['Pdeliv']
        gliP = gData['Pdeliv']
        gndTorq = gData['gndTorq']
        ropeTheta = rData['theta']
        ropeTorq = rData['torq']
    
    #ground roll
    if max(gl.y) >gl.deltar:
        iEndRoll = where(gl.y > gl.deltar)[0][0]
        xRoll = gl.x[iEndRoll]
        tRoll = t[iEndRoll]
        iEndRollData =  where(ti.data['t']>tRoll)[0][0]
    else:
        xRoll = 0  #didn't get off ground
        tRoll = 0
    #misc results 
    vrel =wiv/(env + 1e-6)
    Few = (2-vrel)*tc.invK(vrel) * env**2 / float(wi.rdrum)**3
    Tavg = mean(rp.data[iEndRollData:ti.i]['T'])/gl.W
    # final values     
    yfinal = y[-1]
    yDfinal  = yD[-1]
    if yDfinal < 0.5: yDfinal = 0      
    # max values
    thetaDmax = max(thetaD)
    vmax = max(v)
    vDmax = max(vD) #max acceleration
    Sthmax = max(Sth)
    Tmax =  max(T)/gl.W
    alphaMax = max(alpha[iEndRollData:])
    Lmax = max(L)/gl.W
    gammaMax = max(gamma)
    
    # Comments to user
    print('Controls', control    )
    print('Throttle ramp up time (after slack is out)', tRampUp)
    print('Final height reached: {:5.0f} m, {:5.0f} ft.  Fraction of rope length: {:4.1f}%'.format(yfinal,yfinal/0.305,100*yfinal/float(rp.lo)))
    print('Maximum speed: {:3.0f} m/s, maximum rotation rate: {:3.1f} deg/s'.format(vmax,deg(thetaDmax)))
    print('Maximum Tension factor: {:3.1f}'.format(Tmax))
    print('Average Tension factor: {:3.1f}'.format(Tavg))
    print('Maximum angle of attack: {:3.1f} deg'.format(deg(alphaMax)))
    print('Ground roll: {:5.0f} m, {:5.1f} sec (includes about 1 sec of slack removal)'.format(xRoll,tRoll))
    print('Final vy: {:5.1f} m/s'.format(yDfinal))
    if abs(t[-1] - ti.data[ti.i]['t']) > 5*dt:
        print('\nWarning...the integrator struggled with this model.')
        print('\tIf some of the plots have a time axis that is too short vs others, ')
        print('\t...try making smoother controls.')
    # check whether to keep this run as loop data
    if yfinal < yminLoop and iloop >0:
        print('Bad run, not saving data')
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
        data[iloop]['Tavg'] = Tavg
        data[iloop]['Sthmax'] = Sthmax
        data[iloop]['alphaMax'] = deg(alphaMax)
        data[iloop]['gammaMax'] = deg(gammaMax)
        data[iloop]['thetaDmax'] = deg(thetaDmax)
# plot results for last one
close('all')
plts = plots(path)   

#plts.xy([tData],[xD,yD,v],'time (sec)','Velocity (m/s)',['vx','vy','v'],'Glider velocity vs time') 
#plts.xy([t],[deg(theta,deg(gamma,deg(alpha],'time (sec)','angle (deg)',['pitch','climb','AoA'],'flight angles')  
#plts.xy([t],[env,wiv],'time (sec)','effective speed (m/s)',['engine','winch'],'Engine and winch speeds')
#plts.xy([tData],[L/gl.W,gData['D']/gl.W],'time (sec)','Force/W',['lift','drag'],'Aerodynamic forces')
#plts.xy([tData],],'time (sec)','Torque (Nm) ',['rope','alpha'],'Torques on glider')
#plts.xy([t],[T/gl.W,Few[:itr]/gl.W],'time (sec)','Force/weight',['tension','TC-winch force'],'Forces between objects')
#plts.xy([tData],[Sth],'time (sec)','Throttle setting',['Throttle ',' '],'Throttle')

#plot engine power and torque curves from model parameters
engvel = linspace(0,1.1*en.vpeakP,100)
rpm = engvel/wi.rdrum*60/2/pi*en.diff*en.gear
omegaEng = [r_pm *2*pi/60 for r_pm in rpm]
powr = [en.Pavail(engv)/750 for engv in engvel] #in HP
torq = zeros(len(engvel),dtype = float)
for i in range(1,len(engvel)):
    torq[i] = en.Pavail(engvel[i])/omegaEng[i]*0.74 #leave zero speed at zero torque
torq[0] = torq[1] - (torq[2] - torq[1])*rpm[1]/(rpm[2] - rpm[1]) #Avoid zero speed torq calculation by extrapolating
plts.xy([rpm],[powr,torq],'Engine speed (rpm)','Power (HP), Torque(Ftlbs)',['Pistons power','Pistons torque'],'Engine curves')

#glider position vs time
plts.xy([t],[x,y],'time (sec)','position (m)',['x','y'],'Glider position vs time')
plts.xy([x],[y],'x (m)','y (m)',['x','y'],'Glider y vs x')
#glider speed and angles
#plts.xy([tData,tData,tData,tData,t,t,t,t,tData],[xD,yD,v,deg(alpha,deg(theta,deg(gamma,deg(pl.elev[:itr],deg(theta,gData['L/D']],\
#        'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','angle of attack','pitch','climb','elevator','pitch rate (deg/sec)','L/D'],'Glider velocities and angles')
plts.xy([t,t,tData,tData,t,tData,tData,t],[xD,yD,v,deg(alpha),deg(theta),deg(gamma),deg(elev),deg(thetaD)],\
        'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','angle of attack','pitch','climb','elevator','pitch rate (deg/sec)'],'Glider velocities and angles')

plts.i = 0 #restart color cycle
plts.xyy([tData,t,t,t,tData,tData,tData,t],[v,wiv,y/rp.lo,T/gl.W,L/gl.W,deg(alpha),deg(gamma),deg(thetaD)],\
        [0,0,1,1,1,0,0,0],'time (sec)',['Velocity (m/s), Angles (deg)','Relative forces and height'],['v (glider)',r'$v_r$ (rope)','height/'+ r'$\l_o $','T/W', 'L/W', 'angle of attack','climb angle','rot. rate (deg/sec)'],'Glider and rope')
#lift,drag,forces
plts.xy([tData,tData,t,t],[L/gl.W,D/gl.W,T/gl.W,Few/gl.W],\
        'time (sec)','Forces/Weight',['lift','drag','tension','TC-winch'],'Forces')
#torques
plts.xy([tData],[ropeTorq,Malpha,Me,gndTorq],'time (sec)','Torque (Nm)',['rope','stablizer','elevator','ground'],'Torques')
#Engine, rope and winch
plts.xy([t,t,tData,tData,tData],[env,wiv,vgw,deg(ropeTheta),100*Sth],'time (sec)','Speeds (effective: m/s), Angle (deg), Throttle %',['engine speed','rope speed','glider radial speed','rope angle','throttle'],'Engine and rope')        
#-British units-
plts.xy([t,tData,tData,t,tData,tData],[env/wi.rdrum*60/2/pi*en.diff*en.gear/10,engP/750,engTorq*0.74,wiv*1.94,vgw*1.94,100*Sth],\
    'time (sec)','Speeds (rpm,kts), Torque (ft-lbs), Throttle %',['eng rpm/10', 'pistons HP', 'pistons torque (ftlbs)','rope speed','glider radial speed','throttle'],'Engine British units')        
#Energy,Power
plts.xy([tData],[eData['Edeliv']/1e6,wData['Edeliv']/1e6,rData['Edeliv']/1e6,gData['Edeliv']/1e6,gData['Emech']/1e6],'time (sec)','Energy (MJ)',['to engine','to winch','to rope','to glider','in glider'],'Energy delivered and kept')        
plts.xy([tData],[engP/en.Pmax,winP/en.Pmax,ropP/en.Pmax,gliP/en.Pmax],'time (sec)','Power/Pmax',['to engine','to winch','to rope','to glider'],'Power delivered')        
#Specialty plot for presentation
plts.i = 0 #restart color cycle
plts.xyy([tData,t,t,t,tData,tData,tData,t],[1.94*v,1.94*wiv,y/0.305/10,T/gl.W,L/gl.W,deg(alpha),deg(gamma),deg(thetaD)],\
        [0,0,0,1,1,0,0,0],'time (sec)',['Velocity (kts), Height/10 (ft), Angle (deg)','Relative forces'],\
        ['v (glider)',r'$v_r$ (rope)','height/10','T/W', 'L/W', 'angle of attack','climb angle','rot. rate (deg/sec)'],'Normal rope torque')

if len(tRampUpList) > 1:
    # plot loop results
    heightLoss = data['yfinal'] - max(data['yfinal'])#vs maximum in loop
    plts.i = 0 #restart color cycle
    plts.xyy([data['tRampUp']],[data['xRoll'],10*data['tRoll'],data['yfinal']/rp.lo*100,heightLoss,data['vmax'],data['vDmax']/g,data['Sthmax'],data['Tmax'],data['Lmax'],data['alphaMax'],data['gammaMax'],data['thetaDmax']],\
            [0,0,0,0,0,1,1,1,1,0,0,0],'throttle ramp-up time (sec) ',['Velocity (m/s), Angles (deg), m, sec %',"Relative forces,g's"],\
            ['x gnd roll', 't gnd roll x 10','height/'+ r'$\l_o $%','Height diff',r'$v_{max}$',"max g's",'max throttle',r'$T_{max}/W$', r'$L_{max}/W$', r'$\alpha_{max}$',r'$\gamma_{max}$','rot. max (deg/sec)'],'Flight results vs throttle ramp-up time')
print('Done')