# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Statistical analysis of gust and wind speeds and extreme gust probabilities
'''
import os,sys
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,where,\
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
    
    def readAll(self,inPath):
        '''Reads date, time and wind and gust speed into an array
        Format BNW,BOONE MUNI,1995-03-17 21:27,6,20,9, ... the 6,20,9 gives wind,direction,gust'''
        dtFormat = "%Y-%m-%d %H:%M"
        idata = 0
        lines = readfile(inPath)
        data = zeros(len(lines)-1,dtype = [('year', int),('month', int),('day', int),('hour', int),('min', int),('totsecs', int32),('wind',int),('direc', int),('gust', int)])
        print 'Reading {}s from raw data'.format(type)
        # use the first line to determine if the name has one or two parts
        iAdd = len(lines[1].split()) - 7
        for i,line in enumerate(lines[1:]):
            if mod(i,10000)==0:
                print '{} '.format(i),
            if mod(i,100000)==0:
                print
            info = line.split()
            if len(info) ==7:
                date = info[2+iAdd]
                hrmin = info[3+iAdd]
                wind = info[4+iAdd]
                direc = info[5+iAdd]
                gust = info[6+iAdd]
                d = datetime.strptime(date+' '+hrmin, dtFormat)
                secs = int(time.mktime(d.timetuple()))
                dateInfo = date.split('-')
                data[idata]['year'] = dateInfo[0]
                data[idata]['month'] = dateInfo[1]
                data[idata]['day'] = dateInfo[2]
                timeInfo = hrmin.split(':')
                data[idata]['hour'] = timeInfo[0]
                data[idata]['min'] = timeInfo[1]
                data[idata]['totsecs'] = secs
                data[idata]['wind'] = wind
                data[idata]['direc'] = direc
                data[idata]['gust'] = gust
                idata += 1 
            else:
                print 'Line {}  has a length of only {}.'.format(i,len(info)),info
                sys.exit('Stop!')
        print
        return idata
    
class correlate:
    '''Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/awos/1min.php.  Collates them into one file with 
    the format: day (starting from first datapoint), time, wind(m/s), gust(m/s)
    '''
    
    def __init__(self):
        return 


rdData = readData()
# windPath = 'C:\\Users\\owner\\Downloads\\BooneWind1995-2011minute.txt'
# gustPath = 'C:\\Users\\owner\\Downloads\\BooneGust1995-2011minute.txt'
# outPath = 'C:\\Users\\owner\\Downloads\\BooneWindGust1995-2011minute.dat'
# windPath = 'C:\\Users\\owner\\Downloads\\wind10days.txt'
# gustPath = 'C:\\Users\\owner\\Downloads\\gust10days.txt'
# outPath = 'C:\\Users\\owner\\Downloads\\out10days.dat'
inPath = 'C:\\Users\\owner\\Downloads\\knx95-00onemin.txt'
outPath = 'C:\\Users\\owner\\Downloads\\BooneWindGust1995-2011minute.dat'
nData = rdData.readAll(inPath)
print 'number of datapoints',nData
print 'Done'
        
#         'D:\\Winch launch physics\\results\\testSpy'

