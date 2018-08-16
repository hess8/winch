import IPython

# from IPython import magic
# magic('reset -sf') #comment out this line if not running in Ipython; resets all variables to zero

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
elev, pilot controlled elevator deflection
'''
import matplotlib
# matplotlib.use("GTKAgg")
import os,sys
from numpy import pi,array,zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,where,\
    amin,amax,argmin,argmax,exp,mean,divide 
from numpy.linalg import norm
from numpy import degrees as deg
from numpy import radians as rad

from matplotlib.pyplot import ion,figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend,title,grid
from scipy.integrate import odeint
import logging
logging.basicConfig()
# ion() #turn on interactive mode 
print 'Backend {}, interactive {}'.format(matplotlib.get_backend(),matplotlib.is_interactive())
g = 9.8

def readfile(filepath):
    file1 = open(filepath,'r')
    lines = file1.readlines()
    file1.close()
    return lines

def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

def downsample(arr,maxN):
    if len(arr)/float(maxN) < 2:
        return arr
    else:
        intDown = int(floor(len(arr)/float(maxN)))
        indx = 0
        keep = []
        while indx < len(arr) - 1:
            keep.append(arr[indx])
            indx += intDown      
        return array(keep)
 
def stateSplitMat(S,gl,rp,wi,tc,en,op,pl):
    '''Splits the formal state matrix S (each row a different time) into the various state variables'''
    gl.x      = S[:,0]
    gl.xD     = S[:,1]
    gl.y      = S[:,2]
    gl.yD     = S[:,3]
    gl.theta  = S[:,4]
    gl.thetaD = S[:,5]
    pl.elev   = S[:,6]
    op.Sth    = S[:,7]
    rp.T      = S[:,8]
    wi.v      = S[:,9]
    en.v      = S[:,10]
    return gl,rp,wi,tc,en,op,pl

def writeState(S,gl,rp,wi,tc,en,op,pl,path):
    '''Writes state to a file with one float number per line'''
#     gl.x      = S[:,0]
#     gl.xD     = S[:,1]
#     gl.y      = S[:,2]
#     gl.yD     = S[:,3]
#     gl.theta  = S[:,4]
#     gl.thetaD = S[:,5]
#     pl.elev   = S[:,6]
#     rp.T      = S[:,7]
#     wi.v      = S[:,8]
#     en.v      = S[:,9]
    writefile(['{:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n {:12.8f}\n'\
                .format(gl.x, gl.xD, gl.y, gl.yD, gl.theta, gl.thetaD, pl.elev, op.Sth, rp.T, wi.v, en.v)], path)
    
def readState(path,gl,rp,wi,tc,en,op,pl):
    '''Reads each line into a state variable'''
    lines = readfile(path)
    gl.x      = float(lines[0])
    gl.xD     = float(lines[1])
    gl.y      = float(lines[2])
    gl.yD     = float(lines[3])
    gl.theta  = float(lines[4])
    gl.thetaD = float(lines[5])
    pl.elev   = float(lines[6])
    op.Sth    = float(lines[7])
    rp.T      = float(lines[8])
    wi.v      = float(lines[9])
    en.v      = float(lines[10]) 
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
    op.Sth    = S[7]
    rp.T      = S[8]
    wi.v      = S[9]
    en.v      = S[10]
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
    S[7] = op.Sth
    S[8] = rp.T
    S[9] = wi.v  
    S[10] = en.v        
    return S
    
def pid(var,time,setpoint,c,j,Nint):
    '''The array c contains the [p,i,d] coefficients'''
    # Find the error properties
    err = (var[j] - setpoint)
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
#     if 10<time[j] <13:
#         print 't',time[j],c[0]*err , c[1]*derr , c[2]*interr,c[0]*err + c[1]*derr + c[2]*interr   
    return c[0]*err + c[1]*derr + c[2]*interr    
    
def smooth(data,time,N):
    '''Smooths data with a running average over tsmooth. data and time are arrays
    points are not evenly spaced in time. Smoothing is centered on t.
    Use a triangular weighting'''
    if N == 0: return data
    tsmooth = float(1.0) #sec (float is a precaution: it must not be an integer time)  
#    dt = 0.1 #sec
    smoothed = data
    tfinal = time[-1]
    toExtrapolate = []
    smoothedList = []
    for ism in range(N): #smooth any number of times.
        for i,t in enumerate(time):
            if t > tsmooth/2 and t < tfinal - tsmooth/2:
                iearly = where( time < t - tsmooth/2)[0][-1]
                ilater = where( time > t + tsmooth/2)[0][0] 
            elif t < tsmooth/2: # don't include as many points in average
                iearly = 1 #start here because we need time[iearly-1]
                deli = where( time > t + tsmooth/2)[0][0] - iearly + 1 #+1 to avoid div/zero
                ilater = iearly + int(deli * (t/(tsmooth/2)))
            elif t > tfinal - tsmooth/2:
                ilater = len(data) - 2
                deli = ilater - where( time < t - tsmooth/2)[0][-1] + 1
                iearly = ilater - int(deli * ((tfinal - t)/(tsmooth/2)))
            dsum = 0; totweight = 0 
            for it in range(iearly,ilater+1): 
                weight = (time[it]-time[it-1])*(1-abs(time[it]-t)/(tsmooth/2))
                dsum += smoothed[it] * weight
                totweight += weight  
#             if totweight == 0.0:
#                 'pause'                    
            smoothed[i] = dsum/totweight 
            
    return smoothed
        
class plots:
    def __init__(self,path):
        self.iColor = 0 #counter for plots, so each variable has a different color
        self.path = path
        self.linew = 3.0
        self.colorsList = ['palevioletred', u'#fc4f30', u'#6d904f','darkorange', 'darkviolet', 
        u'#348ABD', u'#e5ae38', u'#A60628', u'#7A68A6', 'mediumaquamarine', u'#D55E00', 'violet',
        u'#CC79A7',  u'#0072B2', u'#30a2da',u'#009E73','peru','slateblue',u'#8b8b8b',] # u'#F0E442',u'#467821','slateblue'      u'#56B4E9',
        return 
    
    def xy(self,holdOpen,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None):
        '''To allow different time (x) arrays, we require the xs to be a list'''
        matplotlib.rcParams.update({'font.size': 14})
        fig, ax0 = subplots(figsize=(14, 7))
        if len(xs)<len(ys): #duplicate first x to every spot in list of correct length
            xs = [xs[0] for y in ys] 
        if len(legendLabels) != len(ys):
            sys.exit('Stop: number of legend labels and curves do not match for curve {}'.format(titlestr))
        for iy,y in enumerate(ys):
            if len(xs[iy]) != len(y):
                sys.exit('Stop: curve {} on plot "{}" has x with dimensions different from y.'.format(iy+1,titlestr))            
            ax0.plot(xs[iy],y,color=self.colorsList[self.iColor],linewidth=self.linew,label=legendLabels[iy])
            self.iColor += 1
            if self.iColor > len(self.colorsList)-1:
                self.iColor = 0
        ax0.set_xlabel(xlbl)
        ax0.set_ylabel(ylbl)
        ax0.grid(color='k', linestyle='--', linewidth=1)
        ax0.legend(loc ='lower right',framealpha=0.5)
#        legend(loc ='upper left')
        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
        ax0.set_ylim([ymin,ymax])
        #Set x limits
        if not xmin is None or not xmax is None:
            ax0.set_xlim([xmin,xmax])
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        if holdOpen: print 'Graphs ready'
        show(block = holdOpen)
#         show()
        
    def xyy(self,holdOpen,xs,ys,yscalesMap,xlbl,ylbls,legendLabels,titlestr,xmin=None,xmax=None):
        '''Allows plotting with two yscales, left (0) and right (1).  
        So yscalesMap is e.g. [0,0,1,0] for 4 datasets'''
        matplotlib.rcParams.update({'font.size': 14})
        linestyles = ['-', '--', '-.', ':']
        fig, ax0 = subplots(figsize=(14, 7))        
        ax1 = ax0.twinx()
        if len(xs)<len(ys): #duplicate first x to every spot in list of correct length
            sys.exit('Stop: number of horizontal (x) arrays does not equal the number of vertical(y) arrays for plot "{}"'.format(titlestr))
        if len(yscalesMap) != len(legendLabels) != len(ys):
            sys.exit('Stop: number of entries in the 0/1 map, legend labels and curves do not match for plot "{}"'.format(titlestr))
        ymaxs = [[],[]]; ymins = [[],[]]       
        for iy,y in enumerate(ys):
            lstyle = '-'
            colorCurve = self.colorsList[self.iColor]
            if 'margin' in legendLabels[iy]:
                lstyle = '--'
            if 'recovery' in legendLabels[iy].lower():
                lstyle = ':'
            if 'vel gust' in legendLabels[iy].lower():
                lstyle = '-'
                colorCurve = 'k'
            if len(xs[iy]) != len(y):
                sys.exit('Stop: curve {} on plot "{}" has x with dimensions different from y.'.format(iy+1,titlestr))
            if yscalesMap[iy] == 0:
                ax0.plot(xs[iy],y,color=colorCurve,linewidth=self.linew,linestyle=lstyle,label=legendLabels[iy])  
            else:
                ax1.plot(xs[iy],y,color=colorCurve,linewidth=self.linew,linestyle=lstyle,label=legendLabels[iy])
            ymaxs[yscalesMap[iy]].append(max(y))
            ymins[yscalesMap[iy]].append(min(y))
            self.iColor += 1 
            if self.iColor > len(self.colorsList)-1:
                self.iColor = 0
        ax0.set_xlabel(xlbl)
        ax0.set_ylabel(ylbls[0])
        ax0.grid(color='k', linestyle='--', linewidth=1)
        ax1.grid(color='k', linestyle='--', linewidth=1)
        ax1.set_ylabel(ylbls[1])
        #Set y limits: force the zeros on both sides to be the same
        
#         xxx
#         min0 = min(0,min(ymins[0]))
#         max0 = 1.1 * max(ymaxs[0])
#         ax0.set_ylim([min0,max0])
#         max1 = 1.1 * max(ymaxs[1])
#         if min0 < 0:
#             shift = max1 * -min0/max0
#         else:
#             shift = 0
#         min1 = min(ymins[1]) - shift
#         
#         xxx
        
        max0 = 1.1 * max(ymaxs[0]); min0 = min(ymins[0])    
        max1 = 1.1 * max(ymaxs[1]); min1 = min(ymins[1])
#         scale0 = max0 - min0; scale1 =  max1 - min1
        if min0 >= 0 and min1 >= 0:
            ymin0 = 0; ymin1 = 0
        elif min0 >= 0 and min1 < 0:
            ymin0 = min1*max0/max1; ymin1 = min1
        elif min0 < 0 and min1 >= 0:
            ymin0 = min0; ymin1 = min0*max1/max0
        else: #both negative
            if abs(min0/max0) > abs(min1/max1):
                ymin0 = min0;
                ymin1 = min0*max1/max0
            else:
                ymin0 = min1*max0/max1
                ymin1 = min1
        ax0.set_ylim([ymin0,max0])
        ax1.set_ylim([ymin1,max1])
        #Set x limits
        if not xmin is None or not xmax is None:
            ax0.set_xlim([xmin,xmax])
        ax0.legend(loc ='upper left',framealpha=0.5)
        ax1.legend(loc ='upper right',framealpha=0.5)
#        legend(loc ='upper left')
#        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
#        ylim([ymin,ymax])
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        if holdOpen: print 'Graphs ready'
        show(block = holdOpen)  
#         show()  

class logger:
    def __init__(self,path,loopWhich):
        self.terminal = sys.stdout
        if loopWhich is None:
            self.log = open('{}{}log.dat'.format(path,os.sep), 'w')
        else:
            self.log = open('{}{}log{}.dat'.format(path,os.sep,loopWhich), 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)          
    
class timeinfo:
    def __init__(self,tStart,tEnd,ntime):
        '''The way to do persistent variables in python is with class variables
        
        The solver does not keep a constant time step, and evaluates the stateDer 
        function usually two or more steps per time step'''
        self.oldt = -1.0
        self.i = -1   #counter for time step.  Must start at -1 because we want to increment it in stateDer, not pl.control.
        self.dt = (tEnd -tStart)/float((ntime-1))
        self.tEnd = tEnd
        #data
        self.data = zeros(ntime,dtype = [('t', float)])
        self.tRelease = None
    

class glider:
    def __init__(self,name,theta0,ntime):
        # parameters
        if name == 'GrobAstir':
            self.vb = 32              # (m/s)  speed of glider at best glide angle
            self.xCG = 0.32  # m back from leading edge
            #wing
            self.xCP = [0.5023,-0.0656, 0.0131,-0.0015,0.000066] #x center of pressure from the wing leading edge vs alpha IN DEGREES for polynomial about zero.0.000066*H3^4 - 0.0015*H3^3 + 0.0131*H3^2 - 0.0656*H3 + 0.5023
            alphaVbestVsWing = rad(1.25) # deg
            alphaStallVsWing = rad(9.0)
            self.alphaStall = alphaStallVsWing - alphaVbestVsWing       #  stall angle vs glider zero
            self.stallLoss = 0.25     # loss of lift (fraction) post stall, see https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20140000500.pdf
            alphaNegStallvsWing = rad(-10) #stall angle for negative lift
            self.alphaNegStall = alphaNegStallvsWing - alphaVbestVsWing
            self.m = 650             #650 kg: max load # 600 kg: Grob, 2 pilots vs 400 for PIK20
            self.W = self.m*9.8          #   weight (N)
            self.n1 = 5.3           # max wing lift factor before danger of structural failure   
            self.n4 = 2.65           # max negative wing lift factor before danger of structural failure         
            self.Co = 0.521             #   Lift coefficient {} of wing at zero glider AoA
            
            self.CLalpha = 5.2    # CL slope /rad: A little less than the airfoil's 2 pi because part of the lift is from the elevator.      
    #         self.Lalpha = self.CLalpha*self.W/self.Co #from xflr5.  
            self.CLmax = self.Co + self.CLalpha * self.alphaStall
            self.vStall =  self.vb*sqrt(self.Co/self.CLmax)   # stall speed fully loaded (m/s)
            self.ralpha =  self.CLalpha/self.Co    # (1/rad)  * wing lift slope/Co 
            #stabilizer
            SsSw = 0.11         # ratio of stabilizer area to wing area (Grob)
            self.Cos = -0.052  # stabilizer lift at best glide alpha, tilt already included
            self.CLSalphas = 3.9  # CL of stab slope /rad from xflr directly
            self.CLSdelelev = 1.8  # CL of stab slope /rad 
            self.I = 1780   #  kg m^2 Grob glider moment of inertia Iyy from xflr5
            self.ls = 4.7           # distance(m) between cg and stabilizer center  
            self.ralphas = SsSw * self.CLSalphas/self.Co   #   (1/rad) ratio of areas * stab lift slope/Co for air-glider pitch moment from angle of attack. 
            self.rsw0 = SsSw * self.Cos/self.Co
            self.rdelev =  SsSw * self.CLSdelelev/self.Co    # (1/rad) ratio of areas * stab lift slope/Co for air-glider pitch moment from elevator deflection           
            self.maxElev = rad(20)   # (rad) maximum elevator deflection
            #other
            self.Agear = 0.02        # drag area (m^2) of main gear
            self.CD = [0.0147,0.0021,0.0003,0.000015] #Drag parameters vs alpha IN DEGREES for polynomial about zero. y = 1E-05x3 + 0.0003x2 + 0.0021x + 0.0147
            self.deltar = 0.02  #ground contact force distance of action (m)
    #        self.deltar = 0.05  #ground contact force distance of action (m)
            self.d_m = 0.25  # distance main wheel to CG (m) 
            self.d_t = 3.3  # distance tail wheel to CG (m) 
            self.theta0 = rad(theta0)
            self.mw = 0.3*self.m
            self.yG = 0.35
            self.yL = 0.45
            
        elif name == 'ASW27':
            sys.exit('stop; havent loaded parameters yet')
          
        
        #variables
#        self.lastvy = 0
        self.oldState = None
        self.vypeaked = False
        self.state = 'onGnd'
        # state variables 
        self.x = 0
        self.xD = 0
        self.y = 0    
        self.yD = 0
        self.theta = 0  #pitch vs horizontal
        self.thetaD = 0
        #data
        self.data = zeros(ntime,dtype = [('x', float),('xD', float),('y', float),('yD', float),('v', float),('vD', float),('yDD', float),('theta', float),('thetaD', float),\
                                    ('gamma', float),('alpha', float),('alphaStall', float),('vgw', float),('L', float),('D', float),('L/D',float),\
                                    ('gndTorq',float),('Fmain',float),('Ftail',float),('Malpha',float),\
                                    ('Pdeliv',float),('Edeliv',float),('Emech',float),('smRecov', float),('smStruct', float),('smStall', float)])
                    
    def findState(self,t,ti,rp):
        '''Determine where in the launch the glider is'''
        gd = self.data[ti.i]
        #One-time switches:
        minHeightRotate = 2
        if not self.vypeaked and self.state != 'onGnd' and self.y > minHeightRotate and self.yD > 10 and self.thetaD < 0:
            print 'vy first peaked at t {:6.2f}'.format(t)            
            self.vypeaked = True #one-time switch
        #state
        if self.y < minHeightRotate and gd['L'] < self.W:
            self.state = 'onGnd'
        elif not self.vypeaked and self.theta < rad(10):
            self.state = 'preClimb'
        elif not self.vypeaked and self.theta >= rad(10):
            self.state ='initClimb' 
        elif self.vypeaked and self.thetaD < 0 :#and t > 7.5: #rad(30) < self.theta  < rad(40) and rp.data[ti.i]['theta'] < rp.thetaCutThr:
            self.state = 'mainClimb'  
        elif rp.data[ti.i]['theta'] >= rp.thetaCutThr :
            self.state = 'prepRelease'
        if self.state != self.oldState:
            print 'Glider state : {} at {:3.1f} s'.format(self.state,t)
            self.oldState = self.state
    
    def gndForces(self,ti,rp):
        del_ymain = self.y + self.d_m*(sin(self.theta) - sin(self.theta0))
        del_ytail = self.y - self.d_t*(sin(self.theta) - sin(self.theta0))
        muS = 0.1  # static friction, including runway bumpiness
        muK = 0.05 # kinetic friction, including r unway bumpiness
        damp = 0 #wheels damp, too.  If more than this it excites the short-period oscillation
        if del_ymain < self.deltar:
            Fmain =  self.W * self.d_t/(self.d_m + self.d_t) * (1-  del_ymain/self.deltar) - damp*gl.yD
        else:
            Fmain = 0
        if del_ytail < self.deltar:
            Ftail =  self.W * self.d_m/(self.d_m + self.d_t) * (1-  (del_ytail)/self.deltar) - damp*gl.yD
        else:
            Ftail = 0
        if self.state == 'onGnd' :
            if self.xD < 0.001 and rp.T < muS*self.W: #not moving             
                Ffric = rp.T
            else: #static broken
                Ffric = muK*self.W
        else:
            Ffric = 0
        return [Fmain, Ftail,Ffric]
                
    def Lnl(self,v,alpha,alphaStab):
        '''Nonlinear lift including post stall'''      
        Lpeak = self.W * (self.vb/self.vStall)**2
        Lsteady = Lpeak * .67
        A = Lpeak-Lsteady        
        w = 0.36 * self.alphaStall
        p = 2
        alphaTrans = rad(6.0-1.25)
        #after stall
        w2 = 1.0 * self.alphaStall
        p2 = 2.0 
        #final drop coefficient
        c = 0.6
        alphaTrans3 = rad(25)
        #negative lift stall
        cn2 = 100.0
        B = 0.5
        shift = rad(2)
        if alpha < self.alphaNegStall + shift:
            x = abs(alpha - (self.alphaNegStall+shift))
#             *Lsteady +
            L = self.W * (-B*(1-exp(-x/w)) + exp(-x/w)*(1 + self.ralpha * alpha + cn2*x**2) )* (v/self.vb)**2 #negative stall 
#             sys.exit('Stop.  Negative lift stall.  Lift not calculated at angle {} deg'.format(deg(alpha)))
        elif self.alphaNegStall <alpha < alphaTrans:
            L = self.W * (1 + self.ralpha * alpha ) * (v/self.vb)**2
        elif alpha < self.alphaStall:
            L = (Lsteady + A * exp( -(abs((alpha- self.alphaStall))/w)**p2  ) ) * (v/self.vb)**2
        elif alpha < alphaTrans3:
            L = (Lsteady + A * exp( -(abs((alpha- self.alphaStall))/w2)**p2  ) )  * (v/self.vb)**2
        else: #a decline toward 90 degrees
            L = (Lsteady * (1 - c * (alpha - alphaTrans3)**1.5) + A * exp( -(abs((alpha- self.alphaStall))/w2)**p2  ) ) * (v/self.vb)**2         
        return L
    
    def xcp(self,alpha):
        '''wing center of pressure measured vs leading edge'''
        alphaMin = -2
        alphaMax = 9.5
        if alpha < rad(alphaMin): alpha = rad(alphaMin)
        elif alpha > rad(alphaMax): alpha = rad(alphaMax)
        return gl.xCP[0] + gl.xCP[1]*deg(alpha) +  gl.xCP[2]*deg(alpha)**2 + gl.xCP[3]*deg(alpha)**3

    def smRecov(self,Vbr,Lbr,alphabr,gammabr,pl): 
        '''Height change after recovery from rope break or power failure'''
        g = 9.8
        Vr = 1.3*self.vStall   #recovery speed
        np0 = 1.5              #recovery pullout g's (constant during pullout, so lift changes during pullout, but not needed to model here, as it's in the diffEqns) 
        p = 0.5
        #Vball is velocity after delay, entering ballistic pushover 
        Vball,gammaBall,ygainDelay = self.delayOutcomes(Vbr,Lbr,alphabr,gammabr,pl) # quantities after delay, beginning of ballistic recovery  
        ygainBallistic = Vr**2/9.8/2 * ((Vball/Vr)**2 -1)
        Vx = cos(gammaBall)*Vball; Vy = sin(gammaBall)*Vball
        if Vx > Vr:  #recovery without dive
            type = 'noDive'
            ygainSlow = sqrt(Vx**2 - Vr**2)/2/g  #slowing from Vx to Vr.
            sm = self.y + ygainBallistic + ygainSlow
        else:
            np = np0
            stall = True
            while stall:
                type = 'Dive'
                gammadive = arctan(sqrt(Vr**2 - Vx**2)/Vx)
                gp = gammadive**p
                Gamma = gp * sin(gp/2)/(np - 0.5 * gp * sin(gp/2)) #big Gamma, see paper
                sm = self.y + ygainDelay + ygainBallistic - Vr**2/9.8 * (Gamma + Gamma**2)
                Vf = Vr + g*gammadive*Vr*sin(gammadive/2)/(np-gammadive/2*sin(gammadive/2))
                if (Vf/self.vStall)**2 >= np + 1:
                    stall = False
                else:
                    np -= 0.1 
        return sm,Vball,gammaBall,ygainDelay,type
    
    def delayOutcomes(self,Vbr,Lbr,alphabr,gammabr,pl):
        '''Assumes a quadratic gamma(t) during delay.  I trust this more than the 12 version below '''
        g = 9.8
        m = self.m
        tdelay = pl.recovDelay
        c1 = (Lbr - m*g*cos(gammabr))/m/Vbr
        c2 = 0.5/m/Vbr * (self.W*self.ralpha * (self.thetaD - c1) - 2*Lbr * g*sin(gammabr)/Vbr)
        gammaAvg = gammabr + 1/2.0*c1*tdelay + 1/3.0*c2*tdelay**2
        Vball = Vbr - g*sin(gammaAvg) * tdelay
        gammaBall = min(pi/2,gammabr + c1*tdelay + c2*tdelay**2)
        ygainDelay = (Vbr - 0.5*g*sin(gammaAvg) * tdelay) * sin(gammaAvg)*tdelay
        return Vball,gammaBall,ygainDelay      
        
class air:
    def __init__(self,vhead,vupdr,hupdr,vgustSM,startGust,widthGust,vGustPeak,gl):
        # wind parameters 
        self.vhead = vhead #headwind
        ## - old updraft quantities
        self.vupdr = vupdr #updraft
        self.hupdr = hupdr #height to turn on updraft
        ## - safety margin quantities
        mu = gl.vStall**2 / 9.8 * gl.CLmax/gl.CLalpha #for safety margin calcs
        k = 0.88 * mu/(5.3 + mu) #for safety margin calcs
        self.vFactors = k * vgustSM  / gl.vStall**2  * gl.CLalpha/gl.CLmax #for safety margin calcs
        ## - dynamic gust simulation quantities 
        self.rStart = None
        self.startGust = startGust
        if not startGust is None:
            self.gustStartType = startGust.split()[1] #either 'm' or 's' for height or time
            self.gustStartValue = float(startGust.split()[0]) # height or seconds float
        self.widthGust = widthGust
        self.vGustPeak = vGustPeak #peak gust velocity, can be positive or negative
        #data
        self.data = zeros(ntime,dtype = [('vGust', float)]) 
        
    def gustLoad(self,gl,v): #for safety margin calcs
        '''Includes gust mediation formula from JAR22'''
        ng = self.vFactors * v
        return ng
    
    def gustDynamic(self,t,vglider,gamma,x,y):
        '''Apply gust perpendicular to glider velocity, in a 1-cos shape.  However
         it takes some distance above the ground for a vertical component of wind to develop.
         If the glider y < widthGust, rather than applying the gust perpendicular to glider velocity it is more of a headwind
         We change this linearly from 0 height to widthGust'''
        if self.startGust is None:
            return 0,0,0,0,0 
        else: 
            if self.rStart is None and ((self.gustStartType == 'm' and y >= self.gustStartValue) or (self.gustStartType == 's' and t >= self.gustStartValue)):
                 #one time switch on gust
                self.rStart = array([x,y])
                return 0,0,0,0,0
            elif self.rStart is None:
                return 0,0,0,0,0
            else:
                d = norm(array([x,y])-self.rStart)
                if d < 2*self.widthGust:
                    vgust = 0.5*vGustPeak*(1-cos(pi*d/self.widthGust)) #can be positive or negative varying with sign of vGustPeak
                else: 
                    vgust = 0
                dStab = d - gl.ls
                if 0 < dStab < 2*self.widthGust:
                    vgustStab = 0.5*vGustPeak*(1-cos(pi*dStab/self.widthGust)) 
                else: 
                    vgustStab = 0  
                thetaGust = gamma + pi * (1- 0.5*min(y/self.widthGust,1))  # 180 degrees when glider is on the ground          
                return -vgust*cos(thetaGust),vgust*sin(thetaGust),d,-vgustStab*cos(thetaGust),vgustStab*sin(thetaGust)  

class rope:
    def __init__(self,lo,thetaMax,breakAngle=None,breakTime=None):
        # Rope parameters  
        self.thetaMax =rad(thetaMax) # release angle of rope  
        self.breakAngle = breakAngle
        self.breakTime = breakTime
        if not breakAngle is None:
            self.breakAngle = rad(breakAngle) # break angle of rope
        if not breakTime is None:
            self.breakTime = breakTime# break angle of rope
        self.broken = False
        self.thetaCutThr =  0.8*self.thetaMax
        self.d  = 0.005     #  rope diameter (m)
        self.A = pi*self.d**2/4      #  rope area (m2)
        self.Apara = 0.02        # drag area (m^2) of parachute along rope
        self.Ys = 30e9             #  2400*9.8/(pi*(0.005/2)^2)/0.035  
                                 #  effective static Young's modulus 30 GPa for rope from Dyneema
                                 #                 datasheet 3.5% average elongation at break,  
                                 #                 average breaking load of 2400 kg (5000 lbs)
        self.hystRate = 100     # hysteresis rate (1/sec) for dynamic stiffness. To turn this effect off, make this rate large, not small
        self.a = 0.7  
        self.b = 0.3            #  vertical distance (m) of rope attachment below CG                     #  horizontal distance (m) of rope attachment in front of CG, for Grob
#        print'b set equal to zero!'
#        self.b = 0.0            #  vertical distance (m) of rope attachment below CG
#        print 'b set equal to 0.4!'
#        self.b = 0.5 

        
        self.lo = lo
        print 'Rope length', lo 
        self.Twl = 8500           #Newtons strength of weak link (Brown)
        self.Cdr = 1.0           # rope drag coefficient
        self.mu = 0.015          # rope linear mass density (kg/meter)
        # state variables 
        self.T = 0
        # data
        self.data = zeros(ntime,dtype = [('T', float),('Tglider', float),('torq', float),\
                       ('theta',float),('Pdeliv',float),('Edeliv',float),('sm', float)]) 
        
        
#     def avgT(self,ti):          
#         tint = 4*self.tau         
#         Nint = min(0,ti.i,ceil(tint/ti.dt))
#         if ti.i >= 0:
#             print ti.i,sum(self.data[ti.i-Nint:ti.i+1]['T'])/(Nint + 1)
#             return sum(self.data[ti.i-Nint:ti.i+1]['T'])/(Nint + 1)       
#         else:
#             return 0
            
    def Tglider(self,thetarope,lenrope):
        g = 9.8 
        return self.T + self.mu * g * lenrope * sin(thetarope)
        
    def thetaRopeGlider(self,t,ti,thetarope,vtrans,lenrope):
        rho = 1.22          # air density (kg/m^3)
        g = 9.8  
        dragCorr = self.Cdr * rho * vtrans**2 * self.d * lenrope / (self.T + 1e-6) / 8           
        weightCorr = 0.5 * self.mu * g *  lenrope * cos(thetarope) / (self.T + self.mu * g * sin(thetarope)+ 1e-6)
        return thetarope + dragCorr + weightCorr
                    
    def chgT(self,vrad,vgw,lenrope):
        '''RHS of diff equation for tension change in terms of rope and glider speeds.  Allows hysteresis'''         
        p = 1    
#        print 'hyst', 1/(1- (abs(vrad - vgw)/self.hystRate/lenrope)**p)
        return rp.Ys*rp.A*(vrad - vgw)/lenrope/(1- (abs(vrad - vgw)/self.hystRate/lenrope)**p)       
            
class winch:
    def __init__(self):
        # Winch parameters  
        self.me = 40              #  Winch effective mass, kg
        self.rdrum = 0.4          # Drum radius (m), needed only for conversion of rpm to engine peak effective speed
        #state variables
        self.v = 0            # Rope uptake speed
        #data
        self.data = zeros(ntime,dtype = [('v',float),('Pdeliv',float),('Edeliv',float)])  #energy delivered to winch by TC

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
    '''Models the torque curve with parameters rather than the power curve, and allows setting the peak torque and peak power 
    rpms independently'''
    def __init__(self,tcUsed,rdrum):
        # Engine parameters  
        self.gear = 1.5    # 2nd gear ratio
#        self.gear = 1.0    # 3rd gear ratio
        self.diff = 3.7 
        self.tcUsed = tcUsed  #model TC, or bypass it (poor description of torque and energy loss)
        self.hp = 300            # engine  horsepower
        self.Pmax = 0.95*750*self.hp        # engine watts.  0.95 is for other transmission losses besides TC
        self.torqMax = 550/0.74*self.Pmax/(0.95*750*390)  # Scaled torque from 390 HP graphs; ft lbs converted to Nm.        
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
    def __init__(self,throttleType,targetT,dipT,thrmax,tRampUp,tHold,ntime):
        #parameters
        self.humanT = 0.5 
        self.targetT = targetT
        self.dipT = dipT
        self.thrmax = thrmax        
        self.tRampUp = tRampUp
        self.angleMax = rad(80) #throttle goes to zero at this rope angle
        self.thrSlack = 0.1
        self.vRopeMax = 30 #m/s 60 kts #Dip period starts when a certain rope speed is reached
        #output
        self.SthTarget = 0
       
        #logical and flow
        self.throttleType = throttleType  
        self.storedState = 'onGnd'
        self.oldvTarget = 0
        self.currvTarget = None
        self.tSwitch = 0
        self.tSlackEnd = None
        self.tHold = tHold
        self.vSlackEnd = 0  #m/s
        self.vRopeReachedTrig = False
        #data
#         self.data = zeros(ntime,dtype = [('Sth', float)])
        #state variables
        self.Sth = 0

    def control(self,t,ti,gl,rp,wi,en):
        tRampUp = self.tRampUp
        tint = 1.0  
        Nint = int(min(ti.i,ceil(tint/ti.dt)))
#        maxrate = 1000 #basically no limit. 
        if self.tSlackEnd is None and gl.xD > self.vSlackEnd:  #one-time event
            self.tSlackEnd = t
        if self.throttleType in ['constT','constTdip']:
            tSlackEnd = self.tSlackEnd              
            tEndRamp = tSlackEnd + tRampUp
            if gl.xD < self.vSlackEnd: #take out slack
                self.SthTarget = self.thrSlack
            else:
                if tSlackEnd  <= t <  tEndRamp and gl.state == 'onGnd': #onGnd so that we can test for sim that is airstart
                    targetT =  self.targetT * (t - self.tSlackEnd)/float(tRampUp)  
                    pp = -1; pd = -1; pint = -2
                elif gl.state == 'onGnd':
                    targetT = self.targetT
                    pp = -8; pd = -32; pint = -16
                elif gl.state == 'prepRelease':
                    targetT = 0
                    pp = -.0; pd = -0; pint = -0                     
                else: #all other states and conditions            
                    targetT = self.targetT
                    pp = -16; pd = -32; pint = -32
                    #overwrite during dip period:
#                   #Dip period starts when a certain rope speed is reached
                    if self.throttleType == 'constTdip':
                        vDipTrigger = 1.2 *gl.vStall #m/s
                        if not self.vRopeReachedTrig:
                            if wi.v > vDipTrigger:
                                print 'vRope reached dip trigger of {} m/s at time {} s'.format(vDipTrigger,t)
                                self.vRopeReachedTrig = True #one-time switch 
                        if self.vRopeReachedTrig and gl.state in ['preClimb','initClimb']:
                            targetT = self.dipT                   
                c = array([pp,pd,pint]) 
                time = ti.data['t']
                Tcontrol = min(self.thrmax,max(0,pid(rp.data['T']/gl.W,time,targetT,c,ti.i,Nint)))
                self.SthTarget = Tcontrol #if not controlling rate of change

        elif self.throttleType == 'preset':
            ### Ramp up, hold, then decrease to steady value
            tSlackEnd = self.tSlackEnd             
            steadyThr = 0.65            
            tRampDown1 = 1 #sec...transition to steady
            tRampDown2 = 120 #longer ramp down
            if tSlackEnd is None:
                self.SthTarget = self.thrSlack
            else:
                tFull = tSlackEnd + tRampUp #when accomplished
                tDown = tFull + self.tHold
                tSlowDown = tDown + tRampDown1
                if  tSlackEnd  <= t <  tFull:
                    self.SthTarget = self.thrSlack + (self.thrmax - self.thrSlack) * (t - tSlackEnd)/float(tRampUp)
                elif tFull < t  < tDown:
                    self.SthTarget =  self.thrmax
                elif tDown <= t < tSlowDown:
                    self.SthTarget = self.thrmax - (self.thrmax -steadyThr) * (t - tDown)/float(tRampDown1)
                else:
                    self.SthTarget = max(0, steadyThr*(1-(t - tSlowDown)/float(tRampDown2)))
        elif self.throttleType == 'constThr':
            tSlackEnd = self.tSlackEnd                        
            if tSlackEnd is None:
                self.SthTarget = self.thrSlack
            else:
                tSlackEnd = self.tSlackEnd                  
                tEndRamp = tSlackEnd + tRampUp
                if  tSlackEnd  <= t <  tEndRamp:
                    self.SthTarget = self.thrSlack + (self.thrmax - self.thrSlack) * (t - tSlackEnd)/float(tRampUp)
                else:
                    self.SthTarget =  self.thrmax
class pilot:
    def __init__(self,pilotType,ntime,ctrltype,setpoint,recovDelay):
        self.Me = 0
        self.MeTarget = 0
        self.err = 0
        self.type = pilotType
        self.ctrltype = ctrltype
        self.setpoint = setpoint
        self.currCntrl = 0
        self.tSwitch = None
        self.tStartClimb = None
        self.pilotStart = None
        self.data = zeros(ntime,dtype = [('err', float),('Me', float),('elev',float)])
        self.humanT = 1.0 #sec
        self.cOld = [0,0,0] #previous PID coefficiens
        self.recovDelay = recovDelay
        #algebraic function
        self.elevTarget = 0   
        #state variable
        self.elev = 0  # elevator deflection, differs from elevTarget by about humanT
        
    def control(self,t,ti,gl,alphaStab):
        #all angles in routines are in radians (vs degrees on the main script input) 
        def alphaControl(t,time,setpoint,ti,Nint):
            if gl.state == 'onGnd' and gl.theta >= gl.theta0 - rad(1): #no elevator authority
                pp =  0; pd =   0; pint =  0
            elif gl.state == 'onGnd' and gl.theta < gl.theta0 - rad(1):
#                pp = -1; pd = -16; pint = -0
                pp = -0; pd = -0; pint = -0
            else:
                if self.tStartClimb is None:
                    self.tStartClimb = t  #one-time switch to climb
                if gl.state == 'preClimb':
                    pp = -512; pd = -1024; pint = -512
                elif gl.state == 'initClimb':
                    pp = -512; pd = -1024; pint = -512
                elif gl.state == 'mainClimb': 
                    pp = -1024; pd = -512; pint = -1024
                elif gl.state == 'prepRelease':
                    pp = 0; pd = 0; pint = 0 
            c = array([pp,pd,pint]) * gl.I
            if not self.tStartClimb is None:
                self.cCurr = c
                cSmooth = self.cCurr - (self.cCurr - self.cOld) * exp(-(t-self.tStartClimb)/self.humanT)
            else:
                cSmooth = c
                self.cOld = c
            al = gl.data['alpha']
            return pid(al,time,setpoint,cSmooth,ti.i,Nint)   

        def vDControl(t,time,setpoint,ti,Nint):
            '''Only derivative control, for use with another control'''
            pp =  0; pd = 256; pint =  0 
            c = array([pp,pd,pint])* gl.I/gl.vb
            varr = gl.data['v']
            return pid(varr,time,setpoint,c,ti.i,Nint)
                
        def vControl(t,time,setpoint,ti,Nint):
            if (gl.state == 'onGnd' and gl.theta >= gl.theta0 - rad(1)) or gl.data[ti.i]['v'] < 25: #no reason to control if going too slow
                pp =  0; pd =   0; pint =  0
            elif gl.state == 'onGnd':
                pp = 0; pd = 0; pint = 0 
            elif gl.state == 'preClimb':
                if gl.thetaD>setpoint:
                    pp = 8; pd = 0; pint = 0
                else:
                    pp = 0; pd = 0; pint = 0
            elif gl.state =='initClimb':
                pp = 64; pd = 0; pint = 0  
            elif gl.state == 'mainClimb': 
                pp = 32; pd = 16; pint = 32
            elif gl.state == 'prepRelease':
                pp =  0; pd =   0; pint =  0   
            c = array([pp,pd,pint])* gl.I/gl.vb
            varr = gl.data['v']
            return pid(varr,time,setpoint,c,ti.i,Nint)
                              
        def limiter(x,xmax):
            if x > xmax:
                x = xmax
            elif x < -xmax:
                x = -xmax
            return x
        tint = self.humanT #sec
        Nint = int(ceil(tint/ti.dt))  
        time = ti.data['t']
        v = gl.data[ti.i]['v']
        alpha = gl.data[ti.i]['alpha']
        gamma = gl.data[ti.i]['gamma']
        crossAngle = rad(self.setpoint[2])  #climb angle to switch from one control to another
        Mdelev =  max(0.001,  gl.ls * gl.W * gl.rdelev * (v/gl.vb)**2) #ratio between moment and elevator deflection.  Avoid zero velocity case with the 0.1.  
                                                                         # alphaStab can be different from wing alpha because of gusts.
        maxMe = gl.maxElev * Mdelev
#         if '' in control:
#             self.Me = 0
#             self.elev = 0
#         else:
#        if self.currCntrl == 0 and (gl.y>10 and (gamma > crossAngle or gl.yD<0)):  #switch to 2nd control # or gl.thetaD < rad(10)  
        if self.currCntrl == 0 and (gl.y>1 and (gamma > crossAngle or gl.data[ti.i]['yDD']<0)):  #switch to 2nd control # or gl.thetaD < rad(10)  
            self.currCntrl = 1 #switches only once
            self.tSwitch = t
            print 'Second control ({},{:3.1f}) turned on  at {:3.1f} s.'.format(self.ctrltype[1],float(self.setpoint[1]),t)
        ctype = self.ctrltype[self.currCntrl]
        setpoint = self.setpoint[self.currCntrl]
        # determine the moment demanded by the control            
        if ctype == 'vD': # speed derivative control only (damps phugoid) 
            self.MeTarget = vDControl(t,time,setpoint,ti,Nint)
        elif ctype == 'v': #target v with setpoint'
            self.MeTarget = vControl(t,time,setpoint,ti,Nint)
        elif ctype == 'alpha': # control AoA  
            self.MeTarget =  alphaControl(t,time,rad(setpoint),ti,Nint)  
        elif ctype == 'alphaVd': # control AoA  
            self.MeTarget =  0.5*alphaControl(t,time,rad(setpoint),ti,Nint)
            self.MeTarget += 3.0*vDControl(t,time,0,ti,Nint)
#             self.MeTarget = vDControl(t,time,0,ti,Nint)
        # implement
        if self.type =='elevControl': 
            self.elevTarget = limiter(self.MeTarget/Mdelev,gl.maxElev) # determine the elevator planned setting 
            self.Me = Mdelev * self.elev #update the moment from the elevator
        elif self.type =='momentControl': # bypass pilot's control of elevator and simply fix the moment required and the elevator to the corresponding angle.
            self.Me = limiter(self.MeTarget,maxMe)
            self.elev = self.Me/Mdelev
        pl.data[ti.i]['Me'] = self.Me  
        pl.data[ti.i]['elev'] = self.elev  

def stateDer(S,t,gl,ai,rp,wi,tc,en,op,pl,save):
    '''First derivative of the state vector'''          
    if rp.data[ti.i]['theta'] > rp.thetaMax: # glider released but the integrator must finish   
        return zeros(len(S))
    else: 
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
        if not rp.broken: #break rope if chosen for simulation:
            if not rp.breakAngle is None and (thetarope > rp.breakAngle or t>rp.breakTime):
                rp.broken = True
                rp.tRelease = t
                print 'Rope broke at rope angle {:4.1f} deg, {:4.1f} s.'.format(deg(thetarope),t)
                print 'Horizontal velocity {:4.1f} m/s.'.format(gl.xD)
        else:
            rp.T = 0.0
        lenrope = sqrt((rp.lo-gl.x)**2 + gl.y**2)
        #steady wind
        vwx = ai.vhead
        if gl.y > ai.hupdr:
            vwy = ai.vupdr
        else:
            vwy = 0.0    
        #glider
        v = sqrt(gl.xD**2 + gl.yD**2) # speed
        vecAir = array([gl.xD+vwx,gl.yD-vwy]) #note minus sign for y, similar to above for vwy.  Vy reduces glider speed vs air
        vAir = norm(vecAir) # speed vs air
        vgw = (gl.xD*(rp.lo - gl.x) - gl.yD*gl.y)/float(lenrope) #velocity of glider toward winch
        vtrans = sqrt(v**2 - vgw**2 + 1e-6) # velocity of glider perpendicular to straight line rope
        thetaRG = rp.thetaRopeGlider(t,ti,thetarope,vtrans,lenrope) # rope angle at glider corrected for rope weight and drag
        Tg = rp.Tglider(thetarope,lenrope) #tension at glider corrected for rope weight
        if gl.xD > 1: #avoid initial zeros problem
            gamma = arctan(gl.yD/gl.xD)  # climb angle.  
            gammaAir = arctan((gl.yD-vwy)/(gl.xD+vwx))  # climb angle vs air
        else:
            gamma = 0
            gammaAir = 0
        #gusts -- only one 1-cos gust per simulation
        vgx,vgy,d,vgxStab,vgyStab = ai.gustDynamic(t,vAir,gammaAir,gl.x,gl.y)       
        vecGustWing = array([vgx,-vgy])  #note minus sign for y, similar to above for vwy
        gammaAirWing = gammaAir; vAirWing = vAir
        if norm(vecGustWing) > 0:
            vecAirWing = vecAir + vecGustWing
            vAirWing = norm(vecAirWing) # speed vs air
            gammaAirWing = arctan(vecAirWing[1]/vecAirWing[0])  # climb angle vs air with gust
        alpha = gl.theta - gammaAirWing # angle of attack for wing
        vecGustStab = array([vgxStab,-vgyStab])  #note minus sign for y, similar to above for vwy
        gammaAirStab = gammaAir; vAirStab = vAir
        if norm(vecGustStab) > 0:
            vecAirStab = vecAir + vecGustStab
            vAirStab = norm(vecAirStab)  
            gammaAirStab = arctan(vecAirStab[1]/vecAirStab[0])  
        alphaStab = gl.theta - gammaAirStab # angle of attack for elevator
        #forces on glider   
        Lwing =  gl.Lnl(vAirWing,alpha,alphaStab) #lift 
        Lstab = gl.W * (gl.rsw0 + gl.ralphas *  alphaStab) * (vAirStab/gl.vb)**2 
        L = Lwing + Lstab 
        CD = gl.CD[0] + gl.CD[1]*deg(alpha) +  gl.CD[2]*deg(alpha)**2 +  gl.CD[3]*deg(alpha)**3
        D = gl.W * CD/gl.Co*(vAirWing/gl.vb)**2 \
           + 0.5 * 1.22 * (gl.Agear * vAirWing**2 + rp.Apara * vgw**2)  # + gl.de*pl.Me #drag  
        wingtorq = Lwing * (gl.xCG - gl.xcp(alpha)) 
        stabtorq = - gl.ls * Lstab
        [Fmain, Ftail, Ffric] = gl.gndForces(ti,rp)
        gndTorq = Fmain*gl.d_m - Ftail*gl.d_t
        M = wingtorq + stabtorq + pl.Me + gndTorq  #torque of air and ground on glider
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
        dotSth = 1/op.humanT * (op.SthTarget-op.Sth) 
        dotT = rp.chgT(wi.v,vgw,lenrope)                       
        if en.tcUsed:
            dotvw =  1/float(wi.me) * (Few - rp.T)
            dotve =  1/float(en.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - Fee)
        else: #no torque converter
            dotvw = 1/float(en.me + wi.me) * (op.Sth * en.Pavail(en.v) / float(en.v) - rp.T)
            dotve = dotvw
        # The ode solver enters this routine usually two or more times per time step.  
        # We advance the time step counter only if the time has changed by close to a nominal time step

        if t - ti.oldt > 1.0*ti.dt: 
            ti.i += 1 
#             print 't,d,vgx,vgy', t,d,vgx,vgy
#             if t > 10: 
#             print 't:{:8.3f}| x:{:8.3f}| xD:{:8.3f}| y:{:8.3f}| yD:{:8.3f}| T:{:8.3f}| Lw:{:8.3f}| Ls:{:8.3f}| vgw:{:8.3f}|state {}|'.format(t,gl.x,gl.xD,gl.y,gl.yD,rp.T,Lwing,Lstab,vgw,gl.state)
#             if t>4.43:
#                 'pause'
#             print 't:{:8.3f}| x:{:8.3f}| xD:{:8.3f}| y:{:8.3f}| yD:{:8.3f}| v:{:8.3f} vAirWing:{:8.3f}| T:{:8.3f}| L:{:8.3f}| alpha: {:8.3f}| gammaW: {:8.3f}| theta: {:8.3f}| thetaD: {:8.3f}| dotthetaD: {:8.3f}|Me: {:8.3f}|alphatorq {:8.3f}| CD {:8.3f}|'.format(t,gl.x,gl.xD,gl.y,gl.yD,v,vAirWing,rp.T,L,deg(alpha),deg(gammaAirWing),deg(gl.theta),deg(gl.thetaD),deg(dotthetaD),pl.Me,alphatorq,CD)
    #             print 'pause'          
    #        print t, 't:{:8.3f}| x:{:8.3f}| xD:{:8.3f}| y:{:8.3f}| yD:{:8.3f}| D/L:{:8.3f}|, L/D :{:8.3f}|'.format(t,gl.x,gl.xD,gl.y,gl.yD,D/L,L/D)
#            print 't,elev',t,deg(pl.elev)           
            if t>4000:
                pl.elev = rad(20)
                print 'FORCING ELEVATOR TO FULL!!'
                
            if not save is None and t > save[1] and not os.path.exists(statefile):
                writeState(S,gl,rp,wi,tc,en,op,pl,save[0])
            ti.data[ti.i]['t']  = t
            gl.data[ti.i]['x']  = gl.x
            gl.data[ti.i]['xD'] = gl.xD
            gl.data[ti.i]['y']  = gl.y
            gl.data[ti.i]['yD'] = gl.yD
            gl.data[ti.i]['v']  = v
            gl.data[ti.i]['vD']  = sqrt(dotxD**2 + dotyD**2)
            gl.data[ti.i]['yDD']  = dotyD
            gl.data[ti.i]['theta']  = gl.theta
            gl.data[ti.i]['thetaD']  = gl.thetaD
            gl.data[ti.i]['gamma']  = gamma
            gl.data[ti.i]['alpha']  = alpha
            gl.data[ti.i]['alphaStall']  = gl.alphaStall #constant stall angle
            if gl.y>0.01: 
                sm,Vball,gammaBall,ygainDelay,type = gl.smRecov(v,L,alpha,gamma,pl)
#                 print 't:{:8.3f} type:{:8s} x:{:8.3f} y:{:8.3f} ygnDelay:{:8.3f} sm:{:8.3f} v:{:8.3f} vball:{:8.3f} gam:{:8.3f} gamball:{:8.3f}  '.format(t,type,gl.x,gl.y,ygainDelay,sm,v,Vball,deg(gamma),deg(gammaBall))
                if gl.y>0.3 or sm > 0: gl.data[ti.i]['smRecov']  = sm #gl.smRecov(v,L,alpha,gamma,pl)
            ngust = ai.gustLoad(gl,v) 
            if Lwing > 0:
                gl.data[ti.i]['smStall']  = (vAirWing/gl.vStall)**2 - Lwing/gl.W - ngust #safety margin vs stall (g's)
                gl.data[ti.i]['smStruct']  = gl.n1*sqrt((1-gl.mw*gl.yG/gl.m/gl.yL*(1-1/gl.n1))) - Lwing/gl.W - ngust #safety margin vs structural damage (g's)   
            else:
                gl.data[ti.i]['smStall']  = (vAirWing/gl.vStall)**2 - ngust #safety margin vs stall (g's).  No negative lift is considered.
                gl.data[ti.i]['smStruct']  = gl.n4*sqrt((1-gl.mw*gl.yG/gl.m/gl.yL*(1-1/gl.n1))) - abs(Lwing)/gl.W - ngust #safety margin vs structural damage (g's)          
#             gl.data[ti.i]['smStruct']  = gl.n1*(1-0) - Lglider/gl.W - ngust #safety margin vs structural damage (g's)            

            gl.data[ti.i]['vgw']  = vgw
            gl.data[ti.i]['L']  = L
            gl.data[ti.i]['D']  = D
            gl.data[ti.i]['gndTorq']  = gndTorq
            gl.data[ti.i]['Malpha']  = wingtorq + stabtorq
            gl.data[ti.i]['Emech'] = 0.5*(gl.m * v**2 + gl.I * gl.thetaD**2) + gl.m  * g * gl.y  #only glider energy here 
            gl.data[ti.i]['Pdeliv'] = Tg * v * cos(thetaRG + gl.theta) 
            gl.data[ti.i]['Edeliv'] = gl.data[ti.i - 1]['Edeliv'] + gl.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            ai.data[ti.i]['vGust'] = norm(vecGustWing)
            rp.data[ti.i]['Pdeliv'] = rp.T * wi.v 
            rp.data[ti.i]['Edeliv'] = rp.data[ti.i - 1]['Edeliv'] + rp.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            rp.data[ti.i]['T'] = rp.T
            rp.data[ti.i]['Tglider'] = Tg            
            rp.data[ti.i]['torq'] = ropetorq
            rp.data[ti.i]['theta'] = thetarope
            rp.data[ti.i]['sm'] = (rp.Twl - Tg)/gl.W
            wi.data[ti.i]['v'] = wi.v
            wi.data[ti.i]['Pdeliv'] = Few * wi.v
            wi.data[ti.i]['Edeliv'] = wi.data[ti.i - 1]['Edeliv'] + wi.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
            en.data[ti.i]['v'] = en.v  
            en.data[ti.i]['Pdeliv'] = op.Sth * en.Pavail(en.v)  
            en.data[ti.i]['torq'] =  op.Sth * en.Pavail(en.v)/(en.v/wi.rdrum*en.gear*en.diff)  #from pistons          
            en.data[ti.i]['Edeliv'] = en.data[ti.i - 1]['Edeliv'] + en.data[ti.i]['Pdeliv'] * (t-ti.oldt) #integrate
#             op.data[ti.i]['Sth'] = op.Sth
            ti.oldt = t
        #---update things that we don't need done ODEint enters stateDer
            gl.findState(t,ti,rp)
            # Update controls
            pl.control(t,ti,gl,alphaStab)
            op.control(t,ti,gl,rp,wi,en) 
#         if t>4:
#             print 'test', dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotelev,dotSth,dotT,dotvw,dotve
        return [dotx,dotxD,doty,dotyD,dottheta,dotthetaD,dotelev,dotSth,dotT,dotvw,dotve]

##########################################################################
#                         Main script
##########################################################################                        
#--- rope lenght
# lo = 6500 * 0.305         #  6500 ft to meters initial rope length (m)
lo = 2000                 # m initial rope length

#--- plotting
smoothed = True
path = 'D:\\Winch launch physics\\results\\test{}m'.format(lo)
if not os.path.exists(path): os.mkdir(path)
#--- loop
loopWhich = 0 #which setpoint element to change in loop
# loopWhich = 'both'
#--- logging
sys.stdout = logger(path,loopWhich) #log screen output to file log.dat

#--- integration
print
tfactor = 16; print 'Time step reduction by factor', tfactor
# tfactor = 160; print 'Time step reduction by factor', tfactor
dt = 0.05/float(tfactor) # nominal time step, sec

#glider
name = 'GrobAstir'
if name == 'GrobAstir':
    theta0 = 6   # deg resting angle of glider on ground 
elif name == 'ASW27':
    theta0 = 4   # deg resting angle of glider on ground 

#--- time
tStart = 0
tEnd = 65 # end time for simulation
ntime = int((tEnd - tStart)/dt) + 1  # number of time steps to allow for data points saved

#--- air
#headwind
vhead = 0
#standard 1-cosine dynamic gust perpendicular to glider path
startGust = None
# startGust = '7 s'
# startGust = '20 m'
widthGust = 12  #m, halfwidth
vGustPeak = 15 * (widthGust/float(110))**(1/float(6))  #m/s
#updraft step function at a single time  
vupdr = 0 
hupdr = 1e6 #m At what height to turn the updraft on for testing
vgustSM = 0  #m/s Keep this zero if simulating dynamic gusts.  Always perpendicular to flight path, for safety margin calcs. 
if abs(vhead) > 0: print 'Headwind', vhead   # m/s
if abs(vupdr) > 0: print 'Updraft of {} m/s, starting at {} m'.format(vupdr,hupdr), vupdr   # m/s

#--- throttle and engine
tRampUp = 2 # time to ramp throttle up
tHold = 0.5
targetT = 1.0
dipT = 0.8
thrmax =  1.0
tcUsed = True   # uses the torque controller
# tcUsed = False  #delivers a torque to the winch determined by Sthr*Pmax/omega
throttleType = 'constT'
# throttleType = 'constTdip'
# throttleType = 'constThr'
#throttleType = 'preset'
if throttleType == 'constThr': print 'Constant throttle',thrmax
elif 'constT' in throttleType: print 'targetT',targetT
if 'dip' in throttleType: print 'dipT',dipT

#--- rope other
ropeThetaMax = 75 #release angle degrees
ropeBreakAngle = 100 #rope angle for break
ropeBreakTime = 100 #sec
#ropeBreakAngle = None
if ropeBreakAngle < ropeThetaMax: print 'Rope break simulation at angle {} deg'.format(ropeBreakAngle)  #
if ropeBreakTime < tEnd: print 'Rope break simulation at time {} sec'.format(ropeBreakTime)  #

#--- pilot
# pilotType = 'momentControl'  # simpler model bypasses elevator...just creates the moments demanded
pilotType = 'elevControl' # includes elevator and response time, and necessary ground roll evolution of elevator
recovDelay = 0.5

#loopParams = [''] #Alpha
# control = ['alpha','alpha']  # Use '' for none
# setpoint = [3 ,3 , 90]  # deg,speed, deg last one is climb angle to transition to final control
control = ['alpha','alphaVd']  # Use '' for none
setpoint = [3,2,30]  # deg,speed, deg last one is climb angle to transition to final control

# control = ['alpha','alpha']  # Use '' for none
# setpoint = [3,3,30]  # deg,speed, deg last one is climb angle to transition to final control
# control = ['alpha','Vd']  # Use '' for none
# setpoint = [3 ,0 , 30]  # deg,speed, deg last one is climb angle to transition to final control

# control = ['alpha','v']  # Use '' for none
# setpoint = [3 ,43 , 30]  # deg,speed, deg last one is climb angle to transition to final control

#control = ['thetaD','v']  # Use '' for none
# setpoint = [5 ,40 , 18]  # deg,speed, deg last one is climb angle to transition to final control
# control = ['alpha','v']
# setpoint = [3,43, 20]  # deg,speed, deg last one is climb angle to transition to final control
#control = ['','vgrad']
#setpoint = [0,35,15]  # deg,speed, deg last one is trigger climb angle to gradually raise the target velocity to setpoint
# control = ['v','v']
# setpoint = [30,30,10]  # deg,speed, deg last one is trigger climb angle to gradually raise the target velocity to setpoint
# control = ['','']
# setpoint = [0,30,10]  # deg,speed, deg last one is trigger climb angle to gradually raise the target velocity to setpoint
#control = ['v','v']
#setpoint = [30,30, 90]  # deg,speed, deg last one is climb angle to transition to final control
# control = ['','']
# setpoint = [0 , 0, 30]  # deg,speed, deglast one is climb angle to transition to final control

# Loop over parameters for study, optimization
loop = True
if loop:
    loopParams = linspace(0,8,45) #Alpha\
else:
    loopParams = [setpoint[0]]
#
data = zeros(len(loopParams),dtype = [('alphaLoop', float),('xRoll', float),('tRoll', float),('yfinal', float),('vmax', float),('vDmax', float),('Sthmax',float),\
                                    ('alphaMax', float),('gammaMax', float),('thetaDmax', float),('Tmax', float),('Tavg', float),('yDfinal', float),('Lmax', float),\
                                    ('smRopeMin', float),('smStallMin', float),('smStructMin', float),('smRecovMin', float)])
yminLoop = 20 #if yfinal is less than this height, the run failed, so ignore this time point
for iloop,param in enumerate(loopParams): 
    alphaLoop = param
    if len(loopParams) > 1:
        if loopWhich == 0:
            resultsDir = path + '\\loop0--{:3.1f},3'.format(param)
            setpoint = [param,3,30]  # deg,speed, deg last one is climb angle to transition to final control
        elif loopWhich == 1:
            resultsDir = path + '\\loop1--3,{:3.1f}'.format(param)
            setpoint = [3,param,30]  # deg,speed, deg last one is climb angle to transition to final control
        elif loopWhich == 'both':
            resultsDir = path + '\\loop01--{:3.1f},{:3.1f}'.format(param,param)
            setpoint = [param,param,30]  # deg,speed, deg last one is climb angle to transition to final control
        if not os.path.exists(resultsDir): os.mkdir(resultsDir)
    else:
        resultsDir = path       
    print '\nInitial pilot control: ({},{:4.1f})'.format(control[0],setpoint[0])    
   
    # create the objects we need from classes
    t = linspace(tStart,tEnd,num=ntime)    
    ti = timeinfo(tStart,tEnd,ntime) 
    gl = glider(name,theta0,ntime)
    ai = air(vhead,vupdr,hupdr,vgustSM,startGust,widthGust,vGustPeak,gl)
    rp = rope(lo,ropeThetaMax,ropeBreakAngle,ropeBreakTime) 
    wi = winch()
    tc = torqconv()
    en = engine(tcUsed,wi.rdrum)
    op = operator(throttleType,targetT,dipT,thrmax,tRampUp,tHold,ntime)
    pl = pilot(pilotType,ntime,control,setpoint,recovDelay)
    # standard nonzero initial conditions, some to avoid div/zero
    gl.xD = 1e-6; gl.yD = 1e-6; gl.v = 1e-6; en.v = en.idle; gl.theta  = rad(theta0)
    # save state
    statefile = resultsDir + '\\state'
    
    ##### Advanced block for air start
    save = None
#     save = [statefile, 6.2 ] # 2nd entry is time to save state
    if not save is None:     
        if os.path.exists(statefile): os.remove(statefile)
    # read initial state
#     loadState = True
    loadState = False
    if loadState: #
        print 'loading state from {}'.format(statefile)     
        gl,rp,wi,tc,en,op,pl = readState(statefile,gl,rp,wi,tc,en,op,pl)
        #Change some of the initial conditions
        rp.T = targetT * gl.W
#         v = 43
#         gl.xD = v * cos(gl.theta); gl.yD = v* sin(gl.theta)
#         gl.y = 1
   
    ##### End advanced block for air start
    
    #initialize state vector to zero  
    S0 = zeros(11)
    #integrate the ODEs
    S0 = stateJoin(S0,gl,rp,wi,tc,en,op,pl)
    S = odeint(stateDer,S0,t,mxstep=500,args=(gl,ai,rp,wi,tc,en,op,pl,save))
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
    #Shorten state data
    t = t[:itr] #shorten
    tmax = t[-1]
    print 'Length of integrator time array:', len(t)
    print 'Length of data arrays:', ti.i    
#     sys.exit('stop')
    
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
    maxSmoothSize = 1000 #downsample because we don't need more than this for plotting, so plotting and smoothing will go faster
    tData = downsample(ti.data[:ti.i]['t'],maxSmoothSize)
    gData = downsample(gl.data[:ti.i],maxSmoothSize)
    wData = downsample(wi.data[:ti.i],maxSmoothSize)
    pData = downsample(pl.data[:ti.i],maxSmoothSize)
    eData = downsample(en.data[:ti.i],maxSmoothSize)
#     oData = op.data[:ti.i]
    rData = downsample(rp.data[:ti.i],maxSmoothSize)
    aData = downsample(ai.data[:ti.i],maxSmoothSize)

    if smoothed:
    #define smoothed data arrays before plotting
        print '\nSmoothing data'
        t = downsample(t,maxSmoothSize)
        x = smooth(downsample(gl.x[:itr],maxSmoothSize),t,1)
        y = smooth(downsample(gl.y[:itr],maxSmoothSize),t,1)
        xD = smooth(downsample(gl.xD[:itr],maxSmoothSize),t,1)
        yD = smooth(downsample(gl.yD[:itr],maxSmoothSize),t,1)
        v = smooth(gData['v'],tData,1)
        vD = smooth(gData['vD'],tData,1) 
        alpha = smooth(gData['alpha'],tData,1)
        theta = smooth(downsample(gl.theta[:itr],maxSmoothSize),t,1)
        thetaD= smooth(downsample(gl.thetaD[:itr],maxSmoothSize),t,1)
        gamma = smooth(gData['gamma'],tData,1)
        elev = smooth(pData['elev'],tData,3)
        Sth = smooth(downsample(op.Sth[:itr],maxSmoothSize),t,2)
        wiv = smooth(downsample(wi.v[:itr],maxSmoothSize),t,2)
        env = smooth(downsample(en.v[:itr],maxSmoothSize),t,2)
        L = smooth(gData['L'],tData,1)
        D = smooth(gData['D'],tData,1)
        T = smooth(downsample(rp.T[:itr],maxSmoothSize),t,1)
        Tg = smooth(rData['Tglider'],tData,1)
        vgw = smooth(gData['vgw'],tData,1)
        Malpha = smooth(gData['Malpha'],tData,1)
        Me = smooth(pData['Me'],tData,1)
        engP = smooth(eData['Pdeliv'],tData,2)
        engTorq = smooth(eData['torq'],tData,2)
        winP = smooth(wData['Pdeliv'],tData,2)
        ropP = smooth(rData['Pdeliv'],tData,2)
        gliP = smooth(gData['Pdeliv'],tData,1)
        gndTorq = smooth(gData['gndTorq'],tData,1)
        ropeTheta = smooth(rData['theta'],tData,1)
        ropeTorq = smooth(rData['torq'],tData,1)
        #don't smooth safety margins more than once!
        smRope = smooth(rData['sm'],tData,1) 
        smStall =smooth(gData['smStall'],tData,1) 
        smStruct = smooth(gData['smStruct'],tData,1) 
        smRecov = smooth(gData['smRecov'],tData,1)
        vGust = smooth(aData['vGust'],tData,1)
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
        thetaD = gl.thetaD[:itr]
        gamma = gData['gamma']
        elev = pData['elev']
        Sth = op.Sth[:itr]
        wiv = wi.v[:itr]
        env = en.v[:itr]
        L = gData['L']
        D = gData['D']
        T = rp.T[:itr]
        Tg = rData['Tglider']
        vgw = gData['vgw']
        Malpha = gData['Malpha']
        Me = pData['Me']
        engP = eData['Pdeliv']
        engTorq = eData['torq']
        winP = wData['Pdeliv']
        ropP = rData['Pdeliv']
        gliP = gData['Pdeliv']
        gndTorq = gData['gndTorq']
        ropeTheta = gData['theta']
        ropeTorq = rData['torq']
        smRope = rData['sm'] 
        smStall = gData['smStall'] 
        smStruct = gData['smStruct'] 
        smRecov = gData['smRecov']
        vGust = aData['vGust']
    #ground roll
    if max(y) >gl.deltar:
        iEndRoll = where(y > gl.deltar)[0][0]
        xRoll = x[iEndRoll]
        tRoll = t[iEndRoll]
        iEndRollData =  where(tData>tRoll)[0][0]
    else:
        xRoll = 0  #didn't get off ground
        tRoll = 0
        iEndRollData = 1
    #misc results 
    vrel =wiv/(env + 1e-6)
    Few = (2-vrel)*tc.invK(vrel) * env**2 / float(wi.rdrum)**3
    Tavg = mean(T[iEndRollData:])/gl.W
    # final values     
    yfinal = y[-1]
    yDfinal  = yD[-1]
#     if yDfinal < 0.5: yDfinal = 0      
    # max values
    ymax = max(y)
    thetaDmax = max(thetaD)
    vmax = max(v)
    vDmax = max(vD) #max acceleration
    Sthmax = max(Sth)
    Tmax =  max(T)/gl.W
    alphaMax = max(alpha[iEndRollData:])
    Lmax = max(L)/gl.W
    gammaMax = max(gamma)
    # sm minimum values taken when glider is 1m or more above the ground 
    minHeight = 1 #m to be airborne
    iminHeight = where(y > minHeight)[0][0]
    tminHeight = t[iminHeight]
    iminHeightData = where(tData>tminHeight)[0][0]
    smRopeMin = min(smRope[iminHeightData:])
    smStallMin = min(smStall[iminHeightData:])
    smStructMin = min(smStruct[iminHeightData:])
    smRecovMin = min(smRecov[iminHeightData:])
    # Comments to user
#    print 'Controls', control    
    print 'Throttle ramp up time (after slack is out)', tRampUp
    print 'Final vy: {:5.1f} m/s, energy added as height'.format(yDfinal)
    if yDfinal > 0:
        ymax += 0.5* yDfinal**2/g
    print 'Final height reached: {:5.0f} m, {:5.0f} ft.  Fraction of rope length: {:4.1f}%'.format(yfinal,yfinal/0.305,100*yfinal/float(rp.lo))
    print 'Maximum speed: {:3.0f} m/s, maximum rotation rate: {:3.1f} deg/s'.format(vmax,deg(thetaDmax))
    print 'Maximum tension factor: {:3.1f}'.format(Tmax)
    print 'Average tension factor: {:3.1f}'.format(Tavg)
    print 'Maximum angle of attack: {:3.1f} deg'.format(deg(alphaMax))
    print 'Ground roll: {:5.0f} m, {:5.1f} sec'.format(xRoll,tRoll)
    if abs(t[-1] - tData[-1]) > 0.05: #sec
        print '\nWarning...the integrator struggled with this model.'
        print '\tIf some of the plots have a time axis that is too short vs others, '
        print '\t...try making smoother controls.'
    # check whether to keep this run as loop data
    if yfinal < yminLoop and iloop >0:
        print 'Bad run, not saving data'
        data[iloop] = data[iloop-1] #write over this data point with the last one. 
    else:
        # Store loop results
        data[iloop]['alphaLoop'] = alphaLoop
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
        data[iloop]['smRopeMin'] = smRopeMin
        data[iloop]['smStallMin'] = smStallMin
        data[iloop]['smStructMin'] = smStructMin
        data[iloop]['smRecovMin'] = smRecovMin
        
    # plot results for each one
    close('all')
    plts = plots(resultsDir)   
    
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
    plts.xy(False,[rpm],[powr,torq],'Engine speed (rpm)','Power (HP), Torque(Ftlbs)',['Pistons power','Pistons torque'],'Engine curves')
    
    ##plot lift curve vs alpha at v_best
    alphaList = linspace(-20,80,100)
    lift = array([gl.Lnl(gl.vb,rad(alph),rad(alph)) for alph in alphaList])
    plts.xy(False,[alphaList],[lift/gl.W],\
            'Angle of attack (deg)','Lift/Weight',[''],'Lift vs angle of attack')
    
    #glider position vs time
    plts.xy(False,[t],[x,y],'time (sec)','position (m)',['x','y'],'Glider position vs time')
    plts.xy(False,[x],[y],'x (m)','y (m)',[''],'Glider y vs x')
    #glider speed and angles
    #plts.xy([tData,tData,tData,tData,t,t,t,t,tData],[xD,yD,v,deg(alpha,deg(theta,deg(gamma,deg(pl.elev[:itr],deg(theta,gData['L/D']],\
    #        'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','angle of attack','pitch','climb','elevator','pitch rate (deg/sec)','L/D'],'Glider velocities and angles')
    plts.xy(False,[t,t,tData,tData,t,tData,tData,t],[xD,yD,v,deg(alpha),deg(theta),deg(gamma),deg(elev),deg(thetaD)],\
            'time (sec)','Velocity (m/s), Angles (deg)',['vx','vy','v','angle of attack','pitch','climb','elevator','pitch rate (deg/sec)'],'Glider velocities and angles')
    
    # plts.iColor = 0 #restart color cycle
    # plts.xyy(False,[tData,t,t,tData,tData,tData,tData,tData,t],[v,wiv,y/rp.lo,Tg/gl.W,L/gl.W,divide(L,D),deg(alpha),deg(gamma),deg(thetaD)],\
    #         [0,0,1,1,1,0,0,0,0],'time (sec)',['Velocity (m/s), Angles (deg), L/D','Relative forces and height'],['v (glider)',r'$v_r$ (rope)','height/'+ r'$\l_o $','T/W at glider', 'L/W', 'L/D','angle of attack','climb angle','rot. rate (deg/sec)'],'Glider and rope')
    #Forces
    plts.xy(False,[tData,tData,t,tData,t],[L/gl.W,D/gl.W,T/gl.W,Tg/gl.W,Few/gl.W],\
            'time (sec)','Forces/Weight',['lift','drag','tension at winch','tension at glider','TC-winch'],'Forces')
    #torques
    plts.xy(False,[tData],[ropeTorq,Malpha,Me,gndTorq],'time (sec)','Torque (Nm)',['rope','stablizer+wing','elevator','ground'],'Torques')
    #Engine, rope and winch
    plts.xy(False,[t,t,tData,tData,t],[env,wiv,vgw,deg(ropeTheta),100*Sth],'time (sec)','Speeds (effective: m/s), Angle (deg), Throttle %',['engine speed','rope speed','glider radial speed','rope angle','throttle'],'Engine and rope')        
    #-British units-
    plts.xy(False,[t,tData,tData,t,tData,t],[env/wi.rdrum*60/2/pi*en.diff*en.gear/10,engP/750,engTorq*0.74,wiv*1.94,vgw*1.94,100*Sth],\
        'time (sec)','Speeds (rpm,kts), Torque (ft-lbs), Throttle %',['eng rpm/10', 'pistons HP', 'pistons torque (ftlbs)','rope speed','glider radial speed','throttle'],'Engine British units')        
    #Energy,Power
    plts.xy(False,[tData],[eData['Edeliv']/1e6,wData['Edeliv']/1e6,rData['Edeliv']/1e6,gData['Edeliv']/1e6,gData['Emech']/1e6],'time (sec)','Energy (MJ)',['to engine','to winch','to rope','to glider','in glider'],'Energy delivered and kept')        
    plts.xy(False,[tData],[engP/en.Pmax,winP/en.Pmax,ropP/en.Pmax,gliP/en.Pmax],'time (sec)','Power/Pmax',['to engine','to winch','to rope','to glider'],'Power delivered')        
    #zoom in on a time range same plot as above
    #Specialty plots for presentations
    # zoom = True
    # if zoom:
    #     t1 = 6 ; t2 = 17
    #     plts.xyy(False,[tData,t,t,tData,tData,tData,tData,tData,tData,t,t],[1.94*v,1.94*wiv,y/0.305/10,Tg/gl.W,L/gl.W,10*deg(alpha),10*deg(gData['alphaStall']),deg(gamma),10*deg(elev),Sth,env/wi.rdrum*60/2/pi*en.diff*en.gear/100],\
    #             [0,0,0,1,1,0,0,0,0,1,0],'time (sec)',['Velocity (kts), Height/10 (ft), Angle (deg)','Relative forces'],\
    #             ['v (glider)',r'$v_r$ (rope)','height/10','T/W', 'L/W', 'angle of attack x10','stall angle x10','climb angle','elev deflection x10','throttle','rpm/100'],'Glider and engine expanded',t1,t2)
    #     plts.iColor = 0 #restart color cycle
    #     plts.xyy(False,[tData,t,tData,tData,tData,tData,tData,tData],[1.94*v,y/0.305,deg(gamma),L/gl.W,smStruct,smStall,smRope,smRecov/0.305],\
    #             [0,0,0,1,1,1,1,0],'time (sec)',['Velocity (kts), Height (ft), Angle (deg)',"Relative forces"],\
    #             ['v','height','Climb angle','L/W','Struct margin','Stall margin','Rope margin','Recovery margin'],'Glider and safety margins',t1,t2)
    # plts.iColor = 0 #restart color cycle
    # plts.xyy(False,[tData,t,t,tData,tData,tData,tData,tData,tData,t,t],[1.94*v,1.94*wiv,y/0.305/10,Tg/gl.W,L/gl.W,10*deg(alpha),10*deg(gData['alphaStall']),deg(gamma),10*deg(elev),Sth,env/wi.rdrum*60/2/pi*en.diff*en.gear/100],\
    #         [0,0,0,1,1,0,0,0,0,1,0],'time (sec)',['Velocity (kts), Height/10 (ft), Angle (deg)','Relative forces'],\
    #         ['v (glider)',r'$v_r$ (rope)','height/10','T/W', 'L/W', 'angle of attack x10','stall angle x10','climb angle','elev deflection x10','throttle','rpm/100'],'Glider and engine')
    # plts.iColor = 0 #restart color cycle
    # plts.xyy(True,[tData,t,tData,tData,tData,tData,tData,tData],[1.94*v,y/0.305/10,deg(gamma),L/gl.W,smStruct,smStall,smRope,smRecov/0.305/10],\
    #         [0,0,0,1,1,1,1,0],'time (sec)',['Velocity (kts), Height (ft), Angle (deg)',"Relative forces"],\
    #         ['v','height/10','climb angle','L/W','Struct margin','Stall margin','Rope margin','Recovery margin'],'Glider and safety margins')
    # Imperial units:
    # if vGustPeak > 0 and ymax >  hGust:
    #     plts.xyy(False,[tData,t,tData,tData,tData,tData,tData,tData,tData],[1.94*v,y/0.305,deg(gamma),L/gl.W,smStruct,smStall,smRope,smRecov/0.305, 1.94*vGust*10],\
    #         [0,0,0,1,1,1,1,0,0],'time (sec)',['Velocity (kts), Height (ft), Angle (deg)',"Relative forces"],\
    #         ['velocity','height','climb angle','L/W','struct margin','stall margin','rope margin','recovery margin','vel gust x10'],'Glider and safety margins')
    # else:
    #     plts.xyy(False,[tData,t,tData,tData,tData,tData,tData,tData],[1.94*v,y/0.305,deg(gamma),L/gl.W,smStruct,smStall,smRope,smRecov/0.305],\
    #         [0,0,0,1,1,1,1,0],'time (sec)',['Velocity (kts), Height (ft), Angle (deg)',"Relative forces"],\
    #         ['velocity','height','climb angle','L/W','struct margin','stall margin','rope margin','recovery margin'],'Glider and safety margins')
    
    #metric units
    zoom = True
    
    if zoom and not startGust is None and 's' in startGust:
        t1 = float(startGust.split()[0])-0.5
        t2 = t1 + 2*widthGust/25.0 + 2.0 -1.5
    else:
        t1 = 5
        t2 = 7
        
    if vGustPeak > 0 and not startGust is None:
        #without safety margins  
        plts.iColor = 0 #restart color cycle
        plts.xyy(False,[tData,t,tData,tData,t,tData],[v*10,y,deg(gamma)*10,L/gl.W,T/gl.W,vGust*10],\
            [0,0,0,1,1,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
            ['velocity x10','height','climb angle x10','L/W','T/W','vel gust x10'],'Winch launch')
        if zoom:
            plts.iColor = 0 #restart color cycle
            plts.xyy(False,[tData,t,tData,tData,t,tData],[v*10,y,deg(gamma)*10,L/gl.W,T/gl.W,vGust*10],\
                [0,0,0,1,1,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
                ['velocity x10','height','climb angle x10','L/W','T/W','vel gust x10'],'Winch launch expanded',t1,t2)    
        #with safety margins (T/W is redundant)
        plts.iColor = 0 #restart color cycle
        plts.xyy(not loop and not zoom,[tData,t,tData,tData,tData,tData,tData,tData,tData,tData],[v*10,y,deg(gamma)*10,L/gl.W,smStruct,smStall,smRope,smRecov, vGust*10],\
            [0,0,0,1,1,1,1,0,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
            ['velocity x10','height','climb angle x10','L/W','struct margin','stall margin','rope margin','recovery margin','vel gust x10'],'Winch launch and safety margins')
        if zoom:
            plts.iColor = 0 #restart color cycle
            plts.xyy(not loop,[tData,t,tData,tData,tData,tData,tData,tData,tData,tData],[v*10,y,deg(gamma)*10,L/gl.W,smStruct,smStall,smRope,smRecov, vGust*10],\
                [0,0,0,1,1,1,1,0,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
                ['velocity x10','height','climb angle x10','L/W','struct margin','stall margin','rope margin','recovery margin','vel gust x10'],'Winch launch and safety margins expanded',t1,t2) 
    else:
        #without safety margins
        if zoom:
            plts.iColor = 0 #restart color cycle
            plts.xyy(False,[tData,t,tData,tData,t],[v*10,y,deg(gamma)*10,L/gl.W,T/gl.W],\
            [0,0,0,1,1,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
            ['velocity x10','height','climb angle x10','L/W','T/W'],'Winch launch expanded',t1,t2)
    
        plts.iColor = 0 #restart color cycle
        plts.xyy(False,[tData,t,tData,tData,t],[v*10,y,deg(gamma)*10,L/gl.W,T/gl.W],\
            [0,0,0,1,1,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
            ['velocity x10','height','climb angle x10','L/W','T/W'],'Winch launch')
    
        #with safety margins
        if zoom:
            plts.iColor = 0 #restart color cycle
            plts.xyy(False,[tData,t,tData,tData,tData,tData,tData,tData],[v*10,y,deg(gamma)*10,L/gl.W,smStruct,smStall,smRope,smRecov],\
                [0,0,0,1,1,1,1,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
                ['velocity x10','height','climb angle x10','L/W','struct margin','stall margin','rope margin','recovery margin'],'Winch launch and safety margins expanded',t1,t2)
        plts.iColor = 0 #restart color cycle
        plts.xyy(not loop,[tData,t,tData,tData,tData,tData,tData,tData],[v*10,y,deg(gamma)*10,L/gl.W,smStruct,smStall,smRope,smRecov],\
            [0,0,0,1,1,1,1,0],'time (sec)',['Velocity (m/s), Height (m), Angle (deg)',"Relative forces"],\
            ['velocity x10','height','climb angle x10','L/W','struct margin','stall margin','rope margin','recovery margin'],'Winch launch and safety margins')
    
# plot loop results (saved in last parameter's folder
if loop:
    heightDiff = data['yfinal'] - min(data['yfinal'])#vs minimum in loop
    plts.iColor = 0 #restart color cycle
    plts.xyy(False,[data['alphaLoop']]*12,[data['xRoll'],10*data['tRoll'],data['yfinal']/rp.lo*100,heightDiff,data['vmax'],data['vDmax']/g,data['Sthmax'],data['Tmax'],data['Lmax'],data['alphaMax'],data['gammaMax'],data['thetaDmax']],\
            [0,0,0,0,0,1,1,1,1,0,0,0],'Angle of attack setting (deg)',['Velocity (m/s), Angles (deg), m, sec, %',"Relative forces,g's"],\
            ['x gnd roll', 't gnd roll x 10','height/rope %','height diff',r'$v_{max}$',"max g's",'max throttle',r'$T_{max}/W$', r'$L_{max}/W$', r'$\alpha_{max}$',r'$\gamma_{max}$','rot. max (deg/sec)'],'Flight (all) vs target angle of attack')
    #fewer results, for presentation
    plts.iColor = 0 #restart color cycle
    plts.xyy(False,[data['alphaLoop']]*4,[data['vmax'],heightDiff,data['gammaMax'],data['Lmax']],\
            [0,0,0,1],'Angle of attack target (deg)',['Velocity (m/s), Angle (deg), Height (m), %',"Relative forces"],\
            [r'$v_{max}$','height diff','climb angle max', r'$L_{max}/W$' ],'Flight vs target angle of attack')
    #With safety margins for presentation
    plts.iColor = 4 #match colors in previous plot
    plts.xyy(True,[data['alphaLoop']]*4,[data['smStructMin'],data['smStallMin'],data['smRopeMin'],data['smRecovMin']],\
            [1,1,1,0],'Target angle of attack (deg)',['Recovery height margin (m)',"Safety margin (g's)"],\
            ['struct margin','stall margin','rope margin','recovery margin'],'Safety margins vs target angle of attack')
    
print 'Done'