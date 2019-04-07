import os,sys,time
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,rint,where,\
    amin,amax,argmin,argmax,exp,mean,mod,int32,sum,log,log10,log1p,float32,transpose,savetxt,loadtxt

from datetime import datetime
import time

    
import matplotlib
from matplotlib.pyplot import ion,figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,\
    legend,title,grid,matshow,colorbar    
    
def readfile(filepath):
    with open(filepath) as f:
        lines = f.read().splitlines() #strips the lines of \n
    return lines

def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

def readAnalysisTask(path):
    ntry = 100
    itry = 1
    while itry <= ntry:
        try:
            lines = readfile('{}\\tasks.dat'.format(path))
            if len(lines) == 0:
                return 'NoTasks'
            f = open('{}\\tasks.dat'.format(path),"r+")
            lines = f.readlines() #strips the lines of \n
            f.seek(0) #pointer to file statrt
            f.writelines(lines[1:],) #write all but first line
            f.truncate() #file is smaller now
            f.close()
            task = lines[0].strip().split() 
            return task
        except:
            time.sleep(0.5) 
            itry += 1           
    return 'ReadFailed'

def writeDone(path,line):
    ntry = 100
    itry = 1
    while itry <= ntry:
        try:
            fd = open(path,'a')
            fd.write(line)
            fd.close()
            return 'OK'
        except:
            time.sleep(0.1) 
            itry += 1           
    return 'Failed'

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
        myEpoch = datetime(1950, 1, 1)
        for iline,line in enumerate(lines[nheaders:]):
            lineOK = True
            if mod(idata,10000)==0 and verbose:
                print(('{} ').format(idata),)
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
                diff = d - myEpoch
            except:
                if verbose:
                    print('Date {} could not be parsed'.format(iline),info)
                lineOK = False
            try:
                totmin = floor(int((diff.days* 24 * 3600 + diff.seconds))/60.0) #minutes since myEpoch
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
        nWindWindEvents = zeros((vmax+1,vmax+1),dtype = int32)
        nGustGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
        nWindGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
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
                    nWindWindEvents[maxWindpast,maxWindfut] += 1
                    maxGustpast =  min(vmax,max(data[index1:i+1]['gust']))
                    maxGustfut =   min(vmax,max(data[i+1:index2+1]['gust']))
                    nGustGustEvents[maxGustpast,maxGustfut] += 1 
                    nWindGustEvents[maxWindpast,maxGustfut] += 1                  
        return nWindWindEvents,nGustGustEvents,nWindGustEvents
    
    def probNextGTE(self,nMat,vmax):
        '''Probability that maximum speed in the next t_2 minutes will be >= v_ (index j), given that the last t1 min was v (index i)'''
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
    def xy(self,holdOpen,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None,ymin=None,ymax=None):
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
        ax0.legend(loc ='upper right',framealpha=0.5)
        #Set axis limits
        if ymin is None: ymin = min([min(y) for y in ys]); 
        if ymax is None: ymax = 1.1*max([max(y) for y in ys]);
        ax0.set_ylim([ymin,ymax])
        if xmin is None: xmin = min(xs[0]); 
        if xmax is None: xmax = max(xs[0])
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
