# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com 
Statistical analysis of gust and wind speeds and extreme gust probabilities
'''
import os,sys
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,rint,where,\
    amin,amax,argmin,argmax,exp,mean,mod,int32,sum,log,log10,log1p,float32,transpose,savetxt,loadtxt
import matplotlib
from matplotlib.pyplot import ion,figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,\
    legend,title,grid,matshow,colorbar
   
from datetime import datetime
import time

def readfile(filepath):
    with open(filepath) as f:
        lines = f.read().splitlines() #strips the lines of \n
    return lines

def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

class readData:
    '''Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/download.phtml?network=AWOS.
    '''
    def __init__(self):
        return 
    
    def readStatePaths(self,inPath,state):
        statePaths = []
        dirList = os.listdir(inPath)
        for path in dirList:
            if path.split('_')[0] == state:
                statePaths.append(path)
        return statePaths

    def readStationData(self,stationPath,verbose):
        '''Reads data from one station file at stationPath.  Puts date, time and wind and gust speed into an array
        Format station,valid,lon,lat,tmpf, dwpf, relh, drct, sknt, p01i, alti, mslp, vsby, gust, skyc1, skyc2, skyc3, skyc4, skyl1, skyl2, skyl3, skyl4, wxcodes, metar'''
        dtFormat = "%Y-%m-%d %H:%M"
        lines = readfile(stationPath)
        data = zeros(len(lines),dtype = [('year', int),('month', int),('day', int),('hour', int),('min', int),('totmin', int32),('wind',int),('direc', int),('gust', int)])        
        idata = 0
        nheaders = 6
        for iline,line in enumerate(lines[nheaders:]):
            lineOK = True
            if mod(idata,10000)==0 and verbose:
                print('{} '.format(idata),)
            if mod(idata,100000)==0 and verbose:
                print()
            try:
                info = line.split(',')
            except:
                print('line {} could not be split'.format(iline),line)
                lineOK = False
            try:
                [date,hrmin] = info[1].split()       
                d = datetime.strptime(date+' '+hrmin, dtFormat)
            except:
                if verbose:
                    print('Date {} could not be parsed'.format(iline),info)
                lineOK = False
            try:
                totmin = floor(int(time.mktime(d.timetuple()))/60.0) #since beginning of epoch
            except:
                print('totmin {} could not be calculated'.format(iline),d)
                lineOK = False                
            try:
                wind = int(rint(float(info[8]))) #knots
            except:
                if verbose:
                    print('wind {}, position 8 was not valid'.format(iline),line)
                lineOK = False
            try:
                if info[13] =='M':
                    gust = 0
                else:
                    gust = int(rint(float(info[13]))) #knots
            except:
                if verbose:
                    print('gust {}, position 13 was not valid'.format(iline),line)
                lineOK = False    
            if lineOK:
                #don't use the following that we aren't using in analysis
#                 dateInfo = date.split('-')
#                 data[idata]['year'] = dateInfo[0]
#                 data[idata]['month'] = dateInfo[1]
#                 data[idata]['day'] = dateInfo[2]        
#                 timeInfo = hrmin.split(':')
#                 data[idata]['hour'] = timeInfo[0]
#                 data[idata]['min'] = timeInfo[1]
                data[idata]['totmin'] = totmin
                data[idata]['wind'] = wind    
                data[idata]['gust'] = gust
                idata += 1  
#             except:
#                 print('Data line {} with length {} is incomplete; skipped'.format(iline,len(info)),info)
        if verbose:print()
        return idata,data[:idata] #don't pass zeros from skipped lines

class correlate:   
    def __init__(self,probFloor):
        self.probFloor = probFloor
     
    def nEventsCount(self,t1,t2,dt,nData,data,vmax,verbose):
        '''For each data point with time greater than t1, find vmax during the past t1 (label _l).  
        Then during the time t2 past this point, find vmax (label _k)'''
        if verbose: print('Counting past and future wind and gust events')
        nWindEvents = zeros((vmax+1,vmax+1),dtype = int32)
        nGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
        n1 = int(rint(t1/float(dt)))
        n2 = int(rint(t2/float(dt)))
        for i in range(nData):
            if mod(i,10000)==0 and verbose:
                print('{} '.format(i),)
            if mod(i,100000)==0 and verbose:
                print()
            # check that there are no time skips in the region of t-t1 to t+t2:
            index1 = i - n1
            index2 = i + n2
            if index1 >= 0 and index2 <= nData - 1:
                if data[index1]['totmin'] != data[i]['totmin'] - t1:
                    if verbose: print('Skips near data point {}: initial time {} not equal to ti-t1:{}'.format(i,data[index1]['totmin'], data[i]['totmin'] - t1)  )
                elif data[index2]['totmin'] != data[i]['totmin'] + t2:
                    if verbose: print('Skips near data point {}: final time {} not equal to ti+t2:{}'.format(i,data[index2]['totmin'], data[i]['totmin'] + t2)  )
                elif len(data[index1:i+1]['wind'])>0 and len(data[i+1:index2+1]['wind'])>0 and len(data[index1:i+1]['gust'])>0 and len(data[i+1:index2+1]['gust'])>0:
                    maxWindpast = min(vmax,max(data[index1:i+1]['wind']))
                    maxWindfut =  min(vmax,max(data[i+1:index2+1]['wind']))
                    nWindEvents[maxWindpast,maxWindfut] += 1
                    maxGustpast =  min(vmax,max(data[index1:i+1]['gust']))
                    maxGustfut =   min(vmax,max(data[i+1:index2+1]['gust']))
                    nGustEvents[maxGustpast,maxGustfut] += 1                      
        return nWindEvents,nGustEvents
    
    def probNextGTE(self,nMat,vmax):
        '''Probability that maximum speed in the next t_2 minutes will be >= v_ (index j), given that the last t1 min had a max of v (index i)'''
        prob = zeros((vmax+1,vmax+1),dtype = float32)
        rowCounts = zeros(vmax+1,dtype = int32)
        for i in range(vmax+1):
            rowCount = sum(nMat[i,:])
            if rowCount>0:
                rowCounts[i] = rowCount
                for j in range(vmax+1):
                    prob[i,j] = sum(nMat[i,j:])/float(rowCount)                
        return prob + self.probFloor,rowCounts
    
    def probNextGTECombine(self,nMat,rowsList,vmax):
        '''Probability that maximum speed in the next t_2 minutes will be >= v_ (index j), given that the last t1 min had a max of v (index i)'''
        prob = zeros((vmax+1),dtype = float32)
        rowsCount = 0
        for irow in rowsList:
            rowsCount += sum(nMat[irow,:])
            for j in range(vmax+1):
                prob[j] += sum(nMat[irow,j:]) 
        if rowsCount>0:             
            return prob/float(rowsCount) + self.probFloor, rowsCount
        else: 
            return zeros((vmax+1),dtype = float32)+self.probFloor, rowsCount

class plots:
    def __init__(self,path):
        self.iColor = 0 #counter for plots, so each variable has a different color
        self.path = path
        self.linew = 3.0
        self.colorsList = ['palevioletred', 'dodgerblue','green', 'darkorange', 'darkviolet','blue', 'red','orange', 
                   'limegreen', 'brown','mediumaquamarine',  'violet','lightcoral', 'olive','tomato','teal','peru','mediumorchid','slateblue','crimson']
        return 
    def xy(self,holdOpen,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None):
        '''To allow different time (x) arrays, we require the xs to be a list'''
        matplotlib.rcParams.update({'font.size': 14})
        fig, ax0 = subplots(figsize=(14,7))
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
        ax0.grid(color='lightgray', linestyle='-', linewidth=1)
        ax0.legend(loc ='lower right',framealpha=0.5)
#        legend(loc ='upper left')
        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
        ax0.set_ylim([ymin,ymax])
        #Set x limits
        if xmin is None or xmax is None:
            xmin = min(xs[0]); xmax = max(xs[0])
        ax0.set_xlim([xmin,xmax])
        #Draw zero line
        ax0.plot([xmin,xmax],[0,0],color='k',linewidth=self.linew)  
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        if holdOpen: print('Graphs ready...pausing after graph "{}"'.format(titlestr))
        show(block = holdOpen)
        return
    
    def matColor(self,holdOpen,mat,xlbl,ylbl,titlestr):
        matplotlib.rcParams.update({'font.size': 10})
        fig = figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(mat,origin='lower')
        fig.colorbar(cax)
        ax.set_xlabel(xlbl) 
        ax.set_ylabel(ylbl)
        ax.set_title(titlestr)
        ax.xaxis.tick_bottom()
        fig.savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
        if holdOpen: print('Graphs ready...pausing after graph "{}"'.format(titlestr))
        show(block = holdOpen)
        return


### read data ###  
close('all')         
# inPath = 'F:\\temp\\temp2'
# inPath = 'F:\\temp\\'
inPath ='F:\\gustsDataRaw\\'
outPath = inPath
t1 = 5 #min 
t2 = 5  #min
dt = 5 #min; the data sampling period
probFloor = 1e-9
vmax = 100 #kts
# verbose = True
# redoPast = True
verbose = False
redoPast = False
rdData = readData() #instance
corr = correlate(probFloor) #instance
# states = ['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','HI','IA','ID','IL','IN','KS','KY','LA','MA','MD','ME',
#           'MI','MN','MO','MS','MT','NC','ND','NE','NH','NJ','NM','NV','NY','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VA','VT',
#           'WA','WI','WV','WY']
# states = ['WY']
# states = ['UT','WY']
states = ['AK']
analysisDir = '{}\\analysis{},{}'.format(outPath,t1,t2)
if not os.path.exists(analysisDir):
    os.mkdir(analysisDir)
stationsAnalyzed = []
for state in states:
    print('==========',state,'=========='); sys.stdout.flush()
    analyzedFile = '{}\\{}_analyzed.dat'.format(analysisDir,state)
    if os.path.exists(analyzedFile):
        stationsAnalyzed += readfile(analyzedFile)
    nWindFile = '{}\\{}_WindEvents.dat'.format(analysisDir,state)
    nGustFile = '{}\\{}_GustEvents.dat'.format(analysisDir,state)
    if os.path.exists(nWindFile) and os.path.exists(nGustFile):
        nWindEvents = loadtxt(nWindFile, dtype=int32)
        nGustEvents = loadtxt(nGustFile, dtype=int32)
#         nWindEvents = [int(i) for i in readfile(nWindFile)]
#         nGustEvents = [int(i) for i in readfile(nGustFile)]
    else:
        nWindEvents = zeros((vmax+1,vmax+1),dtype = int32)
        nGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
    statePaths = rdData.readStatePaths(inPath,state)
    for stationPath in statePaths:
        station = stationPath.split('_')[1]
        if station not in stationsAnalyzed or redoPast:
            print(station,stationPath); sys.stdout.flush()
            nData,data = rdData.readStationData('{}\\{}'.format(inPath,stationPath),verbose)
            newWindEvents, newGustEvents = corr.nEventsCount(t1,t2,dt,nData,data,vmax,verbose)
            nWindEvents += newWindEvents
            nGustEvents += newGustEvents
            #record station as analyzed
            line = '{}\n'.format(station)
            if os.path.exists(analyzedFile):
                fd = open(analyzedFile,'a')
            else:
                fd = open(analyzedFile,'w') #creates file
            fd.write(line)
            fd.close() 
            #write state events counter to file
            savetxt(nWindFile,nWindEvents,fmt='%10d')
    ### correlations for each state ### 
    print('Number of wind events: {}'.format(sum(nWindEvents)))
    print('Number of gust events: {}'.format(sum(nGustEvents)))
#     for i in range(vmax+1):
#         for j in range(vmax+1):
#             print(i,j,'\t',nWindEvents[i,j],'\t',nGustEvents[i,j])
    if sum(nWindEvents)>0:
        nWindEventsDispl = nWindEvents
        nWindEventsDispl[0,0] = 0
        nGustEventsDispl = nGustEvents
        nGustEventsDispl[0,0] = 0
        ### correlated probability ###
        print('probabilities:' )
        #note if an (i,j) element in the end has no events in it, it's given probability of probFloor, a display floor useful when calculating logs
        probNextWindGTE,rowCountw = corr.probNextGTE(nWindEvents,vmax) 
        probNextGustGTE,rowCountg = corr.probNextGTE(nGustEvents,vmax)
        ### combined correlated probability 
        # Add initial max velocities from 0 to 15
        vPastCap = 15
        rows = range(vPastCap+1)
        probNextWindGETcomb,rowsCount = corr.probNextGTECombine(nWindEvents,rows,vmax)
        probNextGustGETcomb,rowsCount = corr.probNextGTECombine(nGustEvents,rows,vmax)
        # Independent probability...simply sum over all intial states
        rows = range(vmax+1)
        probNextWindInd,rowsCount = corr.probNextGTECombine(nWindEvents,rows,vmax)
        probNextGustInd,rowsCount = corr.probNextGTECombine(nGustEvents,rows,vmax)
        
        ### plots ###
        pl = plots(analysisDir) #instance
        ## nEvents
        #for display, log(1+arr) takes the log of 1+array to give contrast and avoid log(0)
        titleStr='{} Wind Correlation (log N+1 events)'.format(state)
        xlbl = 'Next {}-min max speed (kts)'.format(t2)
        ylbl = 'last {}-min max speed (kts)'.format(t1)     
        pl.matColor(False,log10(nWindEventsDispl+1),xlbl,ylbl,titleStr)
        titleStr='{} Gust Correlation (log N+1 events)'.format(state)
        pl.matColor(False,log10(nGustEventsDispl+1),xlbl,ylbl,titleStr)
        ## Probabilities
        titleStr='{} Wind probability (log10)'.format(state)
        pl.matColor(False,log10(probNextWindGTE),xlbl,ylbl,titleStr)
        titleStr='{} Gust probability (log10)'.format(state)
        pl.matColor(False,log10(probNextGustGTE),xlbl,ylbl,titleStr)
        # conditional
        xs = [range(vmax+1)]*3
        ys = [log10(probNextWindInd),log10(probNextWindGETcomb),log10(probNextGustInd),log10(probNextGustGETcomb)]
        titlestr = '{} Probabilities of wind max in next {} minutes'.format(state,t2)
        xlbl = 'Future wind speed (kts)'
        ylbl = 'Probability (log10)'
        strConditionalW = 'after {} min with winds <= {} kts'.format(t1,vPastCap)
        strConditionalG = 'after {} min with gusts <= {} kts'.format(t1,vPastCap)
        legendLabels = ['Wind: independent','Wind: {}'.format(strConditionalW),'Gust: independent', 'Gust: {}'.format(strConditionalG) ]
        pl.xy(True,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None)
        
        # plot(probNextWindGTE); title('Wind probability (independent)'),ylabel('Probability'.format(t1)),xlabel('Next {}-minutes, v max >= speed (kts)'.format(t2)); 
        # plot(probNextGustGTE); title('Gust probability (independent)'),ylabel('Probability'.format(t1)),xlabel('Next {}--minutes, v max >= speed (kts)'.format(t2)); 
        # plot(probNextWindGETcomb); title('Wind probability for {} kt mean intitial max v'.format(mean(rows))),ylabel('Probability'.format(t1)),xlabel('Next {}-minutes, v max >= speed (kts)'.format(t2)); 

print('Done')
