# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Statistical analysis of gust and wind speeds and extreme gust probabilities
'''
import os,sys
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,rint,where,\
    amin,amax,argmin,argmax,exp,mean,mod,int32,sum,log,log10,log1p,float32,transpose
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,\
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
    def __init__(self):
        return    
    
    def nEventsCount(self,t2,dt,nData,data):
        '''For each data point find vmax during the time t2 past this point'''
        print 'Countingfuture wind and gust events'
        gustmax = max(data['gust']) #this is an integer in knots
        nWindEvents = zeros((gustmax+1),dtype = int32)
        nGustEvents = zeros((gustmax+1),dtype = int32)
        n2 = int(rint(t2/float(dt)))
        for i in range(nData):
            if mod(i,10000)==0:
                print '{} '.format(i),
            if mod(i,100000)==0:
                print
            # check that there are no time skips in the region of t-t1 to t+t2:
            index2 = i + n2
            if index2 <= nData - 1:
                if data[index2]['totmin'] != data[i]['totmin'] + t2:
                    print 'Skips near data point {}: indexed time {} not equal to ti+t2:{}'.format(i,data[index2]['totmin'], data[i]['totmin'] + t2)  
                else: #all OK
                    maxWindfut =  max(data[i+1:index2]['wind'])
                    nWindEvents[maxWindfut] += 1
                    maxGustfut =  max(data[i+1:index2]['gust'])
                    nGustEvents[maxGustfut] += 1                  
        return nWindEvents,nGustEvents,gustmax
    
    def probNextHigherThan(self,nMat,gustmax):
        '''Probability that maximum speed in the next t_2 minutes will be greater than v_ (index j), independent of previous velocity'''
        prob = zeros((gustmax+1),dtype = float32)
        arrayCount = sum(nMat)
        for i in range(gustmax+1):
            prob[i] = sum(nMat[i,j:])/float(arrayCount)                        
        return prob  
    
class correlate:
    
    def __init__(self):
        return
     
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
                    print 'Skips near data point {}: indexed time {} not equal to ti-t1:{}'.format(i,data[index1]['totmin'], data[i]['totmin'] - t1)  
                elif data[index2]['totmin'] != data[i]['totmin'] + t2:
                    print 'Skips near data point {}: indexed time {} not equal to ti+t2:{}'.format(i,data[index2]['totmin'], data[i]['totmin'] + t2)  
                else: #all OK
                    maxWindpast = max(data[index1:i]['wind'])
                    maxWindfut =  max(data[i+1:index2]['wind'])
                    nWindEvents[maxWindpast,maxWindfut] += 1
                    maxGustpast = max(data[index1:i]['gust'])
                    maxGustfut =  max(data[i+1:index2]['gust'])
                    nGustEvents[maxGustpast,maxGustfut] += 1                   
        return nWindEvents,nGustEvents,gustmax
    
    def probNextHigherThan(self,nMat,gustmax):
        '''Probability that maximum speed in the next t_2 minutes will be greater than v_ (index j), given that the last t1 min had a max of v (index i)'''
        prob = zeros((gustmax+1,gustmax+1),dtype = float32)
        rowCounts = zeros(gustmax+1,dtype = int32)
        for i in range(gustmax+1):
            rowCount = sum(nMat[i,:])
            if rowCount>0:
                rowCounts[i] = rowCount
                for j in range(gustmax+1):
                    prob[i,j] = sum(nMat[i,j:])/float(rowCount)                
        return prob,rowCounts

class plots:
    
    def __init__(self):
        return
### read data ###           
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
### independent events probability ###
#count the events where gusts greater than v occur in a time t2 ahead of ti
ind = independent()
nWindEvents,nGustEvents,gustmax = ind.nEventsCount(t2,dt,nData,data)
probNextWindHigherThan = ind.probNextHigherThan(nWindEvents,gustmax)


### correlations ### 
corr = correlate() #instance
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
close('all')
fig1 = matshow(log10(nWindEventsDispl+1));colorbar();title('Wind Correlation'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next {}-min max speed (kts)'.format(t2))
fig2 = matshow(log10(nGustEventsDispl+1));colorbar();title('Gust Correlation'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next {}-min max speed (kts)'.format(t2)) 
# show(fig1)
# show(fig2)

### correlated probability ###
print 'probabilities:' 
#note if a row (velocity) has no events in it, it's given probability of zero
probNextWindHigherThan,rowCountw = corr.probNextHigherThan(nWindEvents,gustmax)
probNextGustHigherThan,rowCountg = corr.probNextHigherThan(nGustEvents,gustmax)
factor1 = 100
minWprob = min(probNextWindHigherThan[where(probNextWindHigherThan>0)]) 
minGprob = min(probNextGustHigherThan[where(probNextGustHigherThan>0)]) 
relProbNextWindHigherThanDispl = probNextWindHigherThan/minWprob*factor1
relProbNextGustHigherThanDispl = probNextGustHigherThan/minGprob*factor1

for i in range(gustmax+1):
    for j in range(gustmax+1):
        if relProbNextWindHigherThanDispl[i,j] == 0:
            relProbNextWindHigherThanDispl[i,j] = 1 #set unknown relative probabilities (no data) to a floor factor1 below the lowest calculated probability
        if relProbNextGustHigherThanDispl[i,j] == 0:
            relProbNextGustHigherThanDispl[i,j] = 1 #set unknown relative probabilities (no data) to a floor factor1 below the lowest calculated probability
        print i,j,'\t',probNextWindHigherThan[i,j], '[{}]'.format(rowCountw[i]),'\t',probNextGustHigherThan[i,j], '[{}]'.format(rowCountg[i])
fig3 = matshow(log10(relProbNextWindHigherThanDispl));colorbar();title('Wind probability (relative)'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next {}-min max speed (kts)'.format(t2))
fig4 = matshow(log10(relProbNextGustHigherThanDispl));colorbar();title('Gust probability (relative)'),ylabel('last {}-min max speed (kts)'.format(t1)),xlabel('Next {}-min max speed (kts)'.format(t2)) 
show(fig1); show(fig2);show(fig3);show(fig4)




print 'Done'
