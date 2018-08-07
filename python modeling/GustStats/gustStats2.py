# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Statistical analysis of gust and wind speeds and extreme gust probabilities
'''
import os,sys
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,rint,where,\
    amin,amax,argmin,argmax,exp,mean,mod,int32
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend,title,grid
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
        # use the first line to determine if the name has one or two parts
        iAdd = len(readfile(inPaths[0])[1].split()) - 7
        for inPath in inPaths:
            nLines += len(readfile(inPath))-1
        print 'Total number of lines in {} files: {}'.format(len(inPaths),nLines)
        data = zeros(nLines,dtype = [('year', int),('month', int),('day', int),('hour', int),('min', int),('totmin', int32),('wind',int),('direc', int),('gust', int)])
        for inPath in inPaths:
            print 'Reading from file {}'.format(inPath)
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
                    data[idata]['totmin'] = totmin #since beginning of all data
                    data[idata]['wind'] = wind
                    data[idata]['direc'] = direc
                    data[idata]['gust'] = gust
                    idata += 1 
                else:
                    print 'Data line {} has a length of only {}; skipped'.format(idata,len(info)),info
                   # sys.exit('Stop!')
        print
        return idata,data
    
class correlate:
    '''
    '''
    
    def __init__(self):
        self.data = data
        self.ndata = data
        return
     
    def v_kv_lcounts(self,t1,t2,dt,nData,data):
        '''For each data point with time greater than t1, find vmax during the past t1 (label _l).  
        Then during the time t2 past this point, find vmax (label _k)'''
        gustmax = max(data['gust']) #this is an integer in knots
        nEvents = zeros((gustmax,gustmax),dtype = int32)
        n1 = int(rint(t1/float(dt)))
        n2 = int(rint(t2/float(dt)))
        for i in range(self.ndata):
            # determine if there are no time skips in the region of t-t1 to t+t2:
            index1 = i - n1 -1 
            index2 = i + n2
            OK = True
            if index1 >= 0 and index2 <= nData - 1:
                if data[index1]['totmin'] != data[i]['totmin'] - t1:
                    OK == False
                    print 'Skips near data point {}.  
                elif data[index2]['totmin'] != data[i]['totmin'] + t1:
                    OK == False
                    print 'Skips near data point {}.
                if OK:
                    maxv = max()
                    maxv` = max()
                    
                 
        
        return nEvents
           

rdData = readData()
# windPath = 'C:\\Users\\owner\\Downloads\\BooneWind1995-2011minute.txt'
# gustPath = 'C:\\Users\\owner\\Downloads\\BooneGust1995-2011minute.txt'
# outPath = 'C:\\Users\\owner\\Downloads\\BooneWindGust1995-2011minute.dat'
# windPath = 'C:\\Users\\owner\\Downloads\\wind10days.txt'
# gustPath = 'C:\\Users\\owner\\Downloads\\gust10days.txt'
# outPath = 'C:\\Users\\owner\\Downloads\\out10days.dat'
inPaths = ['C:\\Users\\owner\\Downloads\\knx95-00onemin.txt',
           'C:\\Users\\owner\\Downloads\\knx01-05onemin.txt',
            'C:\\Users\\owner\\Downloads\\knx06-10onemin.txt']
outPath = 'C:\\Users\\owner\\Downloads\\BooneWindGust1995-2010minute.dat'
nData,data,dt = rdData.readAll(inPaths)
print 'number of complete datapoints',nData
corr = correlate() #instance
t1 = 30 #min 
t2 = 5  #min
dt = 1 #min; the data sampling period
nEvents = corr.v_kv_lcounts(t1,t2,dt,nData,data)
print 'Done'
