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
import os
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tanh,ceil,floor,where,amin,amax,argmin,argmax
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend,title
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

class plots:
    def __init__(self,path):
        self.i = 0 #counter for plots, so each variable has a different color
        self.path = path        
        return 
    
    def xy(self,xs,ys,xlbl,ylbl,legendLabels,titlestr):
        '''To allow different time (x) arrays, we require the xs to be a list'''
        if len(xs)<len(ys): #duplicate first x to every spot in new list
            xs = [xs[0] for y in ys] 
        colorsList = [u'#30a2da', u'#fc4f30', u'#6d904f','darkorange', u'#8b8b8b',
              u'#348ABD', u'#e5ae38', u'#A60628', u'#7A68A6', u'#467821', u'#D55E00', 'darkviolet',
              u'#CC79A7',  u'#0072B2',u'#56B4E9', u'#009E73','peru','slateblue'] # u'#F0E442'
        figure()
        for iy,y in enumerate(ys):
            plot(xs[iy],y,color=colorsList[self.i],linewidth=2.0,label=legendLabels[iy])
            self.i += 1
            if self.i > len(colorsList)-1:
                self.i = 0
            xlabel(xlbl)
            ylabel(ylbl)
        legend(loc='lower right')
        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
        ylim([ymin,ymax])
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        show()
    
class timeinfo:
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
        self.m = 600             # kg Grob, 2 pilots vs 400 for PIK20
        self.W = self.m*9.8          #   weight (N)
        self.Q = 20             #   L/D
        self.alphas = 6*3.14/180         #   stall angle vs glider zero
        Co = 0.75             #   Lift coefficient {} at zero glider AoA
        self.Lalpha = 2*pi*self.W/Co
        self.I = 600*6/4   #   Grob glider moment of inertia, kgm^2, scaled from PIK20E
        self.palpha = 1.8 * self.W        #   change in air-glider pitch moment (/rad) with angle of attack 
        self.dv = 3.0            #   drag constant ()for speed varying away from vb
        self.dalpha = 40        #   drag constant (/rad) for glider angle of attack away from zero. 
        self.de = 0.025          #   drag constant (/m) for elevator moment

        # state variables 
        self.x = 0
        self.xD = 0
        self.y = 0    
        self.yD = 0
        self.theta = 0
        self.thetaD = 0 
        
        #data
        self.data = zeros(ntime,dtype = [('t', float),('x', float),('xD', float),('y', float),('yD', float),\
                                    ('v', float),('theta', float),('alpha', float),('L', float),('Malpha',float),\
                                    ('D', float),('rptorq',float),('Pdeliv',float),('Edeliv',float),('Emech',float)])
        return
    
class rope:
    def __init__(self):
        # Rope parameters  
        Rrope = 0.005/2.0     #  rope radius (m)
        self.A = pi*Rrope**2      #  rope area (m2)
        self.Y = 30e9             #  2400*9.8/(pi*(0.005/2)^2)/0.035  
                                 #  effective Young's modulus 30 GPa for rope from Dyneema
                                 #                 datasheet 3.5% average elongation at break,  
                                 #                 average breaking load of 2400 kg (5000 lbs)
        self.a = 0.2             #  horizontal distance (m) of rope attachment in front of CG
        self.b = 0.1             #  vertial distance (m) of rope attachment below CG
        self.lo = 8000 * 0.305         #  initial rope length (m)
#        self.lo = 1000         #  initial rope length (m)
        # state variables 
        self.T = 0
        
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
         self.Ko = 13             #TC capacity (rad/sec/(Nm)^1/2)  K = 142 rpm/sqrt(ft.lb) = 142 rpm/sqrt(ft.lb) * 2pi/60 rad/s/rpm * sqrt(0.74 ftlb/Nm) = 12.8 (vs 12 in current model!)
         self.dw = 0.13
         #data
         self.data = zeros(ntime,dtype = [('Pdeliv',float),('Edeliv',float)])  #energy delivered to TC by engine (impeller) rotation
     def invK(self,vrel):
         return 1/float(self.Ko) * tanh((1-vrel)/self.dw)

class engine:
    def __init__(self,rdrum):
        # Engine parameters  
        self.hp = 350             # engine horsepower
        self.Pmax = 750*self.hp        # engine watts
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
        # state variables 
        self.v = 0            #engine effective speed (m/s)
        self.Few = 0           # effective force between engine and winch (could go in either engine or winch or in its own class)
         #data
        self.data = zeros(ntime,dtype = [('Pdeliv',float),('Edeliv',float)]) #energy delivered to engine rotating mass by pistons
    
    def Pavail(self,ve):            # power curve
        vr = ve/float(self.vpeak)
        return self.Pmax*(self.pe1 * vr + self.pe2 * (vr)**2 - self.pe3 * (vr)**3)
        
class operator:
    def __init__(self,ntime):
        self.Sth = 0
        self.data = zeros(ntime,dtype = [('t', float),('Sth', float)])
        return
    def control(self,t,gl,rp,wi,en):
        tRampUp = 1  #seconds
        tHold = 20
        tDown = tRampUp + tHold
        tRampDown = 60
        thrmax = 1.0 
        if t <= tRampUp:
            self.Sth =  thrmax/float(tRampUp) * t
        elif tRampUp < t < tDown:
            self.Sth =  thrmax
        elif t >= tDown:
            self.Sth = max(0,thrmax * (1-(t-tDown)/float(tRampDown)))
         
class pilot:
    def __init__(self,ntime,ctrltype,setpoint):
        self.Me = 0
        self.err = 0
        self.ctrltype = ctrltype
        self.setpoint = setpoint
        self.data = zeros(ntime,dtype = [('t', float),('err', float),('Me', float)])
        return
        
    def control(self,t,ti,gl): 
        setpoint = 0.0
        tint = 0.5 #sec
        Nint = ceil(tint/ti.dt)        
        cGo = True
        if '' in control:        
            cGo = False
        if len(self.ctrltype)>1:  # We have two types of control
            angleSwitch = self.setpoint[2]
            gamma = arctan(gl.yD/gl.xD)*180/pi
            if gamma < angleSwitch:
                ctype = self.ctrltype[0]
                setpoint = self.setpoint[0]
            else:
                ctype = self.ctrltype[1]
                setpoint = self.setpoint[1]                
        else:
            ctype = self.ctrltype
            setpoint = self.setpoint
        var = gl.data[ctype]
        if ctype == 'v':
            pp = 0; pd = 0; pint = 0
            var = var/gl.vb
            if gl.y < 1:
                cGo = False
        elif ctype =='alpha':
            pp = -100; pd = -20; pint = -40
        time = gl.data['t']
        if cGo: 
            #Find the error properties
            err = (var[ti.i] - setpoint)
            self.err = err
            #for the derivative, take the last two time steps
            if ti.i >= 2:  
                derr = (var[ti.i] - var[ti.i-2])/(time[ti.i]- time[ti.i-2]) 
            else:
                derr = 0
            #for the integral, last Nint average
            if ti.i >= Nint:            
                interr = sum(var[ti.i-Nint : ti.i])/(Nint + 1) - setpoint
            else:
                interr = 0
            self.Me = gl.I*(pp*err + pd*derr + pint*interr) #normalize error constants by Iglider. 
            pl.data[ti.i]['t'] = t #store error
            pl.data[ti.i]['err'] = err
            pl.data[ti.i]['Me'] = self.Me
        else:        
            self.Me = 0       

def stateDer(S,t,gl,rp,wi,tc,en,op,pl):
    '''First derivative of the state vector'''
    gl,rp,wi,tc,en,op,pl = stateSplitVec(S,gl,rp,wi,tc,en,op,pl)
    if gl.xD < 1e-6:
        gl.xD = 1e-6 #to handle v = 0 initial
    if en.v < 1e-6:
        en.v = 1e-6 #to handle v = 0 initial
    #----algebraic functions----#
    # Update controls
    op.control(t,gl,rp,wi,en) 
    pl.control(t,ti,gl)
    #rope
    thetarope = arctan(gl.y/float(rp.lo-gl.x))
    lenrope = sqrt((rp.lo-gl.x)**2 + gl.y**2)
    #glider
    v = sqrt(gl.xD**2 + gl.yD**2) # speed
    gamma = arctan(gl.yD/gl.xD)  # climb angle.   
    alpha = gl.theta - gamma # angle of attack
    L = (gl.W + gl.Lalpha*alpha) * (v/gl.vb)**2 #lift       
    D = L/float(gl.Q)*(1 + gl.Q*gl.dalpha*alpha**2)# + gl.de*pl.Me #drag
    if alpha > gl.alphas: #stall
        L = 0.75*L
        D = L/float(gl.Q)
    M = (-gl.palpha*alpha + pl.Me) #torque of air on glider
    ropetorq = rp.T*sqrt(rp.a**2 + rp.b**2)*sin(arctan(rp.b/float(rp.a))-gl.theta-thetarope) #torque of rope on glider
    vgw = (gl.xD*(rp.lo - gl.x) - gl.yD*gl.y)/float(lenrope) #velocity of glider toward winch
    
    #winch-engine
    vrel = wi.v/en.v
    Fee = tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3     
    Few = Fee * (2-vrel)  # effective force between engine and winch through torque converter
    #----derivatives of state variables----#
    dotx = gl.xD    
    dotxD = 1/float(gl.m) * (rp.T*cos(thetarope) - D*cos(gamma) - L*sin(gamma)) #x acceleration
    doty = gl.yD
    if gl.y < 1.0 and L < gl.W: #on ground 
        dotyD = 0 
        dottheta = 0
        dotthetaD = 0
    else:
        dotyD = 1/float(gl.m) * (L*cos(gamma) - rp.T*sin(thetarope) - D*sin(gamma) - gl.W) #y acceleration
        dottheta = gl.thetaD    
        dotthetaD = 1/float(gl.I) * (ropetorq + M)    
    dotT = rp.Y*rp.A*(wi.v - vgw)/float(lenrope)
    dotvw =  1/float(wi.me) * (Few - rp.T)
    dotve =  1/float(en.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - Few / (2 - vrel))

    # store data from this time step for use in controls or plotting.  The ode solver enters this routine
    # usually two or more times per time step.  We advance the time step counter only if the time has changed 
    # by close to a nominal time step    
    if t - ti.oldt > 0.9*ti.dt: 
        ti.i += 1 
        
        gl.data[ti.i]['t']  = t
        gl.data[ti.i]['x']  = gl.x
        gl.data[ti.i]['xD'] = gl.xD
        gl.data[ti.i]['y']  = gl.y
        gl.data[ti.i]['yD'] = gl.yD
        gl.data[ti.i]['v']  = v
        gl.data[ti.i]['theta']  = gl.theta
        gl.data[ti.i]['alpha']  = alpha
        gl.data[ti.i]['L']  = L
        gl.data[ti.i]['D']  = D
        gl.data[ti.i]['Malpha']  = -gl.palpha*alpha
        gl.data[ti.i]['rptorq'] = ropetorq
#        gl.data[ti.i]['Emech'] = 0.5*(gl.m * v**2 + gl.I * gl.thetaD**2 + en.me*en.v**2 + wi.me*wi.v**2) + gl.m  * g * gl.y 
        gl.data[ti.i]['Emech'] = 0.5*(gl.m * v**2 + gl.I * gl.thetaD**2) + gl.m  * g * gl.y  #only glider energy here
#        if rp.T > 1000:
#            print 'pause'
        gl.data[ti.i]['Pdeliv'] = rp.T * vgw 
        gl.data[ti.i]['Edeliv'] = gl.data[ti.i - 1]['Edeliv'] + gl.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        wi.data[ti.i]['Pdeliv'] = Few * wi.v
        wi.data[ti.i]['Edeliv'] = wi.data[ti.i - 1]['Edeliv'] + wi.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        en.data[ti.i]['Pdeliv'] = op.Sth * en.Pavail(en.v)         
        en.data[ti.i]['Edeliv'] = en.data[ti.i - 1]['Edeliv'] + en.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
        op.data[ti.i]['t']   = t
        op.data[ti.i]['Sth'] = op.Sth
        
        ti.oldt = t
    return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotT,dotvw,dotve]

##########################################################################
#                         Main script
##########################################################################                        
tStart = 0
tEnd = 90      # end time for simulation
dt = 0.05       #nominal time step, sec
path = 'D:\\Winch launch physics\\results\\test'  #for saving plots
#path = 'D:\\Winch launch physics\\results\\aoa control Grob USA winch'  #for saving plots
#control = 'alpha'  # Use '' for none
#setpoint = 2*pi/180   #alpha, 2 degrees

control = ['alpha','v']
#setpoint = [2*pi/180 , 1.0, 30]  #last one is climb angle to transition to final control
setpoint = [2*pi/180  , 35, 60]  #last one is climb angle to transition to final control


#control =/ 'v'  # Use '' for none
#setpoint = 1.0                    # for velocity, setpoint is in terms of vbest: vb


ntime = (tEnd - tStart)/dt + 1   # number of time steps expected
t = linspace(tStart,tEnd,num=ntime)

# create the objects we need from classes
ti = timeinfo(tStart,tEnd,ntime) 
gl = glider(ntime)
rp = rope()
wi = winch()
tc = torqconv()
en = engine(wi.rdrum)
op = operator(ntime)
pl = pilot(ntime,control,setpoint)
  
S0 = zeros(9)
# nonzero initial conditions
gl.xD = 1e-6  #to avoid div/zero
gl.ve = 1e-6  #to avoid div/zero
if control[0] == 'alpha':
    gl.theta = setpoint[0]
else:
    gl.theta = 2*pi/180   #initial AoA, 2 degrees

#integrate the ODEs
S0 = stateJoin(S0,gl,rp,wi,tc,en,op,pl)
S = odeint(stateDer,S0,t,args=(gl,rp,wi,tc,en,op,pl,))
#Split S (now a matrix with state variables in columns and times in rows)
gl,rp,wi,tc,en,op,pl = stateSplitMat(S,gl,rp,wi,tc,en,op,pl)

#If yD drops below zero, the sim went too long.  Remove time points after that. 
if min(gl.yD) < 0 :
    negyD = where(gl.yD < 0)[0]
    itr = negyD[0]-1 #ode solver index for release time
    t = t[:itr] #shorten
    negyDData = where(gl.data['yD']<0)[0]
    ti.i = negyDData[0] - 1  #data index for release time
else:
    itr = len(t)

# plot results
tData = gl.data[:ti.i]['t']
gData = gl.data[:ti.i]
wData = wi.data[:ti.i]
pData = pl.data[:ti.i]
eData = en.data[:ti.i]
oData = op.data[:ti.i]

close("all")
plts = plots(path)   
plts.xy([t],[gl.x[:itr],gl.y[:itr]],'time (sec)','position (m)',['x','y'],'Glider position vs time')
plts.xy([tData],[gData['xD'],gData['yD'],gData['v']],'time (sec)','Velocity (m/s)',['vx','vy','v'],'Glider velocity vs time')
gamma = arctan(gl.yD[:itr]/gl.xD[:itr]) 
alpha = gl.theta[:itr] - gamma  
plts.xy([t],[180/pi*gl.theta[:itr],180/pi*gamma,180/pi*alpha],'time (sec)','angle (deg)',['pitch','climb','AoA'],'flight angles')  
plts.xy([t],[en.v[:itr],wi.v[:itr]],'time (sec)','effective speed (m/s)',['engine','winch'],'Engine and winch speeds')
plts.xy([tData],[gData['L']/gl.W,gData['D']/gl.W],'time (sec)','Force/W',['lift','drag'],'Aerodynamic forces')
plts.xy([tData],[gData['rptorq'],gData['Malpha']],'time (sec)','Torque (Nm) ',['rope','alpha'],'Torques on glider')
plts.xy([tData],[100*pData['err'],pData['Me']],'time (sec)',' ',['errorx100 (rad/s)','elevator moment (Nm)'],'Pilot')
vrel =wi.v/en.v 
Few = (2-vrel)*tc.invK(vrel) * en.v**2 / float(wi.rdrum)**3
plts.xy([t],[rp.T[:itr]/gl.W,Few[:itr]/gl.W],'time (sec)','Force/weight',['tension','TC-winch force'],'Forces between objects')
plts.xy([tData],[oData['Sth']],'time (sec)','Throttle setting',['Throttle ',' '],'Throttle')
   #glider speed and angles
plts.xy([tData,tData,tData,t,t,t,t],[gData['xD'],gData['yD'],gData['v'],180/pi*gl.theta[:itr],180/pi*gamma,180/pi*alpha,180/pi*gl.thetaD[:itr]],\
        'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','pitch','climb','AoA','pitch rate (deg/sec)'],'Glider velocities and angles')
   #lift,drag,forces
plts.xy([tData,tData,t,t],[gData['L']/gl.W,gData['D']/gl.W,rp.T[:itr]/gl.W,Few[:itr]/gl.W],\
        'time (sec)','Force/W ',['lift','drag','tension','TC-winch'],'Forces')
    #Engine and winch
plts.xy([t,t,tData],[en.v[:itr],wi.v[:itr],100*oData['Sth']],'time (sec)','Speeds (effective: m/s), Throttle',['engine speed','winch speed','throttle %'],'Engine and winch')        
    #Energy,Power
plts.xy([tData],[eData['Edeliv']/1e6,wData['Edeliv']/1e6,gData['Edeliv']/1e6,gData['Emech']/1e6],'time (sec)','Energy (MJ)',['to engine','to winch','to glider','in glider'],'Energy delivered and kept')        
plts.xy([tData],[eData['Pdeliv']/en.Pmax,wData['Pdeliv']/en.Pmax,gData['Pdeliv']/en.Pmax],'time (sec)','Power/Pmax',['to engine','to winch','to glider'],'Power delivered')        
# Comments to user
yfinal = gData['y'][-1]  
vyfinal  = gData['yD'][-1]
if vyfinal < 0.5: vyfinal = 0
print 'Final height reached: {:5.0f} m, {:5.0f} ft.  Fraction of rope length: {:4.1f}%'.format(yfinal,yfinal/0.305,100*yfinal/float(rp.lo))
print 'Maximum speed: {:3.0f} m/s, maximum rotation rate: {:3.1f} deg/s'.format(max(gData['v']),180/pi*max(gl.thetaD[:itr]))
print 'Maximum Tension factor: {:3.2f}'.format(max(rp.T[:itr]/gl.W))
print 'Final vy: {:5.1f} m/s'.format(vyfinal)

if abs(t[-1] - gl.data[ti.i]['t']) > 5*dt:
    print '\nWarning...the integrator had a harder time with this model.'
    print '\tIf some of the plots have a time axis that is too short vs others, '
    print '\t...try making smoother controls.'

print 'Done'





