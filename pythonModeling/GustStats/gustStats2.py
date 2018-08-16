# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Statistical analysis of gust and wind speeds and extreme gust probabilities
'''
import os,sys
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,rint,where,\
    amin,amax,argmin,argmax,exp,mean,mod,int32,sum,log,log10,log1p,float32,transpose
import matplotlib
from matplotlib.pyplot import ion,figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,\
    legend,title,grid,matshow,colorbar
   
from datetime import datetime
import time

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

class readData:
    '''Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/awos/1min.php.  Collates them into one file with 
    the format: day (starting from first datapoint), time, wind(kts), gust(kts)
    '''
    
    def __init__(self):
        return 
    
    def readAll(self,inPaths):
        '''Reads Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/awos/1min.php. Best to download in pieces of 5 or 6 years.  Determines date, time and wind and gust speed into an array
        Format BNW,BOONE MUNI,1995-03-17 21:27,6,20,9, ... the 6,20,9 gives wind,direction,gust'''
        dtFormat = "%Y-%m-%d %H:%M"
        idata = 0
        #determine number of datapoints from all the files
        nLines = 0
        # use the first line to determine if the name has one or two parts, and the minutes baseline of the starting point
        info1 = readfile(inPaths[0])[1].split()
        iAdd = len(info1) - 7
        date1 = info1[2+iAdd]
        hrmin1 = info1[3+iAdd]        
        d1 = datetime.strptime(date1+' '+hrmin1, dtFormat)
        totmin1 = floor(int(time.mktime(d1.timetuple()))/60.0)        
        for inPath in inPaths:
            nLines += len(readfile(inPath))-1
        print 'Total number of lines in {} files: {}'.format(len(inPaths),nLines)
        data = zeros(nLines,dtype = [('year', int),('month', int),('day', int),('hour', int),('min', int),('totmin', int32),('wind',int),('direc', int),('gust', int)])
        for inPath in inPaths:
            print 'Reading from file {}'.format(inPath)
#             print 'Gust defined by max(wind,gust) reported'
            lines = readfile(inPath)
            for i,line in enumerate(lines[1:]):
                if mod(idata,10000)==0:
                    print '{} '.format(idata),
                if mod(idata,100000)==0:
                    print
                info = line.split()
                if len(info) == 7 + iAdd:
                    date = info[2+iAdd]
                    hrmin = info[3+iAdd]
                    wind = info[4+iAdd]
                    direc = info[5+iAdd]
                    gust = info[6+iAdd]
                    d = datetime.strptime(date+' '+hrmin, dtFormat)
                    totmin = floor(int(time.mktime(d.timetuple()))/60.0)
                    dateInfo = date.split('-')
                    data[idata]['year'] = dateInfo[0]
                    data[idata]['month'] = dateInfo[1]
                    data[idata]['day'] = dateInfo[2]
                    timeInfo = hrmin.split(':')
                    data[idata]['hour'] = timeInfo[0]
                    data[idata]['min'] = timeInfo[1]
                    data[idata]['totmin'] = totmin - totmin1 #since beginning of all data
                    data[idata]['wind'] = wind
                    data[idata]['direc'] = direc
#                     data[idata]['gust'] = max(wind,gust)
                    data[idata]['gust'] = gust
                    idata += 1 
                else:
                    print 'Data line {} has a length of only {}; skipped'.format(idata,len(info)),info
                   # sys.exit('Stop!')
        print
        return idata,data
    
class independent:
    def __init__(self,probFloor):
        self.probFloor = probFloor
        return    
    
    def nEventsCount(self,t2,dt,nData,data):
        '''For each data point find vmax during the time t2 past this point'''
        print 'Countingfuture wind and gust events'
        gustmax = max(data['gust']) #this is an integer in knots
        nWindEvents = zeros((gustmax+1),dtype = int32)
        nGustEvents = zeros((gustmax+1),dtype = int32)
        n2 = int(rint(t2/float(dt)))
        for i in range(nData):
#         for i in range(1000):
            if mod(i,10000)==0:
                print '{} '.format(i),
            if mod(i,100000)==0:
                print
            # check that there are no time skips in the region of ti+t2:
            index2 = i + n2
            if index2 <= nData - 1:
                if data[index2]['totmin'] != data[i]['totmin'] + t2:
                    print 'Skips near data point {}: final time {} not equal to ti+t2:{}'.format(i,data[index2]['totmin'], data[i]['totmin'] + t2)  
                else: #all OK
                    
                    maxWindfut =  max(data[i+1:index2]['wind'])
                    nWindEvents[maxWindfut] += 1
                    maxGustfut =  max(data[i+1:index2]['gust'])
                    nGustEvents[maxGustfut] += 1  
#                     print i, maxWindfut,maxGustfut               
        return nWindEvents,nGustEvents,gustmax
    
    def probNextGTE(self,nMat,gustmax):
        '''Probability that maximum speed in the next t_2 minutes will be >= vi, independent of previous velocity'''
        prob = zeros((gustmax+1),dtype = float32)
        arrayCount = sum(nMat)
        for i in range(gustmax+1):
            prob[i] = sum(nMat[i:])/float(arrayCount)                        
        return prob + self.probFloor

class correlate:   
    def __init__(self,probFloor):
        self.probFloor = probFloor
     
    def nEventsCount(self,t1,t2,dt,nData,data):
        '''For each data point with time greater than t1, find vmax during the past t1 (label _l).  
        Then during the time t2 past this point, find vmax (label _k)'''
        print 'Counting past and future wind and gust events'
        gustmax = max(data['gust']) #this is an integer in knots
        nWindEvents = zeros((gustmax+1,gustmax+1),dtype = int32)
        nGustEvents = zeros((gustmax+1,gustmax+1),dtype = int32)
        n1 = int(rint(t1/float(dt)))
        n2 = int(rint(t2/float(dt)))
        for i in range(nData):
            if mod(i,10000)==0:
                print '{} '.format(i),
            if mod(i,100000)==0:
                print
            # check that there are no time skips in the region of t-t1 to t+t2:
            index1 = i - n1
            index2 = i + n2
            if index1 >= 0 and index2 <= nData - 1:
                if data[index1]['totmin'] != data[i]['totmin'] - t1:
                    print 'Skips near data point {}: initial time {} not equal to ti-t1:{}'.format(i,data[index1]['totmin'], data[i]['totmin'] - t1)  
                elif data[index2]['totmin'] != data[i]['totmin'] + t2:
                    print 'Skips near data point {}: final time {} not equal to ti+t2:{}'.format(i,data[index2]['totmin'], data[i]['totmin'] + t2)  
                else: #all OK
                    maxWindpast = max(data[index1:i]['wind'])
                    maxWindfut =  max(data[i+1:index2]['wind'])
                    nWindEvents[maxWindpast,maxWindfut] += 1
                    maxGustpast = max(data[index1:i]['gust'])
                    maxGustfut =  max(data[i+1:index2]['gust'])
                    nGustEvents[maxGustpast,maxGustfut] += 1                   
        return nWindEvents,nGustEvents,gustmax
    
    def probNextGTE(self,nMat,gustmax):
        '''Probability that maximum speed in the next t_2 minutes will be >= v_ (index j), given that the last t1 min had a max of v (index i)'''
        prob = zeros((gustmax+1,gustmax+1),dtype = float32)
        rowCounts = zeros(gustmax+1,dtype = int32)
        for i in range(gustmax+1):
            rowCount = sum(nMat[i,:])
            if rowCount>0:
                rowCounts[i] = rowCount
                for j in range(gustmax+1):
                    prob[i,j] = sum(nMat[i,j:])/float(rowCount)                
        return prob + self.probFloor,rowCounts
    
    def probNextGTECombine(self,nMat,rowsList,gustmax):
        '''Probability that maximum speed in the next t_2 minutes will be >= v_ (index j), given that the last t1 min had a max of v (index i)'''
        prob = zeros((gustmax+1),dtype = float32)
        rowsCount = 0
        for irow in rowsList:
            rowsCount += sum(nMat[irow,:])
            for j in range(gustmax+1):
                prob[j] += sum(nMat[irow,j:])              
        return prob/float(rowsCount) + self.probFloor, rowsCount

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
        if holdOpen: print 'Graphs ready...pausing after graph "{}"'.format(titlestr)
        show(block = holdOpen)
#         show()

        return
### read data ###  
close('all')         
rdData = readData()
inPaths = ['C:\\Users\\owner\\Downloads\\knx95-00onemin.txt',
           'C:\\Users\\owner\\Downloads\\knx01-05onemin.txt',
            'C:\\Users\\owner\\Downloads\\knx06-10onemin.txt']
# inPaths = ['C:\\Users\\owner\\Downloads\\knxOneMonthTest.txt']
# inPaths = ['C:\\Users\\owner\\Downloads\\knx95-00onemin.txt']
outPath = 'C:\\Users\\owner\\Downloads\\'
nData,data = rdData.readAll(inPaths)
print 'number of complete datapoints',nData

t1 = 30 #min 
t2 = 5  #min
dt = 1 #min; the data sampling period
probFloor = 1e-9
### independent events probability ###
#count the events where gusts >= v occur in a time t2 ahead of ti
ind = independent(probFloor)
nWindEvents,nGustEvents,gustmax = ind.nEventsCount(t2,dt,nData,data)
probNextWindGTEind = ind.probNextGTE(nWindEvents,gustmax)
probNextGustGTEind = ind.probNextGTE(nGustEvents,gustmax)
show(block = False)
### correlations ### 
corr = correlate(probFloor) #instance
nWindEvents,nGustEvents,gustmax = corr.nEventsCount(t1,t2,dt,nData,data)
print 'Number of wind events {}; skipped {}'.format(sum(nWindEvents), nData - sum(nWindEvents))
print 'Number of gust events: {}'.format(sum(nWindEvents))
for i in range(gustmax+1):
    for j in range(gustmax+1):
        print i,j,'\t',nWindEvents[i,j],'\t',nGustEvents[i,j]

nWindEventsDispl = nWindEvents
nWindEventsDispl[0,0] = 0
nGustEventsDispl = nGustEvents
nGustEventsDispl[0,0] = 0
#for display, log(1+arr) takes the log of 1+array to give contrast and avoid log(0)
fig1 = matshow(log10(nWindEventsDispl+1));colorbar();title('Wind Correlation'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next -minutes, v max >= speed (kts)'.format(t2))
fig2 = matshow(log10(nGustEventsDispl+1));colorbar();title('Gust Correlation'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next -minutes, v max >= speed (kts)'.format(t2)) 

### correlated probability ###
print 'probabilities:' 
#note if a row (velocity) has no events in it, it's given probability of zero
probNextWindGTE,rowCountw = corr.probNextGTE(nWindEvents,gustmax) 
probNextGustGTE,rowCountg = corr.probNextGTE(nGustEvents,gustmax)
fig3 = matshow(log10(probNextWindGTE));colorbar();title('Wind probability (log10)'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next {}-min max speed (kts)'.format(t2))
fig4 = matshow(log10(probNextGustGTE));colorbar();title('Gust probability (log10)'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next {}-min max speed (kts)'.format(t2)) 
show(block = False)
### combined correlated probability 
# Add initial max velocities from 0 to 15
vPastCap = 15
rows = range(vPastCap+1)
probNextWindGETcomb,rowsCount = corr.probNextGTECombine(nWindEvents,rows,gustmax)
probNextGustGETcomb,rowsCount = corr.probNextGTECombine(nGustEvents,rows,gustmax)
# Plots
pl=plots(outPath) #instance
xs = [range(gustmax+1)]*3
ys = [log10(probNextWindGTEind),log10(probNextGustGTEind),log10(probNextWindGETcomb),log10(probNextGustGETcomb)]
titlestr = 'Probabilities of wind max in next {} minutes'.format(t2)
xlbl = 'Future wind speed (kts)'
ylbl = 'Probability (log10)'
strConditionalW = 'after {} min with winds <= {} kts'.format(t1,vPastCap)
strConditionalG = 'after {} min with gusts <= {} kts'.format(t1,vPastCap)
legendLabels = ['Wind: independent','Gust: independent', 'Wind: {}'.format(strConditionalW),'Gust: {}'.format(strConditionalG) ]
pl.xy(True,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None)

# plot(probNextWindGTE); title('Wind probability (independent)'),ylabel('Probability'.format(t1)),xlabel('Next {}-minutes, v max >= speed (kts)'.format(t2)); 
# plot(probNextGustGTE); title('Gust probability (independent)'),ylabel('Probability'.format(t1)),xlabel('Next {}--minutes, v max >= speed (kts)'.format(t2)); 
# plot(probNextWindGETcomb); title('Wind probability for {} kt mean intitial max v'.format(mean(rows))),ylabel('Probability'.format(t1)),xlabel('Next {}-minutes, v max >= speed (kts)'.format(t2)); 


print 'Done'
