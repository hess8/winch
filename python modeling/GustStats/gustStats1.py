# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Statistical analysis of gust and wind speeds and extreme gust probabilities
'''
import os,sys
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tan,tanh,ceil,floor,where,\
    amin,amax,argmin,argmax,exp,mean,mod
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

class extractToFile:
    '''Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/awos/1min.php.  Collates them into one file with 
    the format: day (starting from first datapoint), time, wind(m/s), gust(m/s)
    '''
    
    def __init__(self):
        return 
    
    def readVar(self,inPath):
        '''Reads date, time and wind or gust speed into an array
        Format BNW BOONE MUNI 1995-01-01 01:34 15'''
        dtFormat = "%Y-%m-%d %H:%M"
        out = []
        lines = readfile(inPath)
        for i,line in enumerate(lines[1:]):
            if mod(i,10000)==0:
                print '{} '.format(i),
            if mod(i,100000)==0:
                print
            info = line.split()
            if i>5120000: print 'info',i,info
            date = info[3]
            hrmin = info[4]
            var = info[5]
            d = datetime.strptime(date+' '+hrmin, dtFormat)
            secs = time.mktime(d.timetuple())
            out.append([date,hrmin,secs,var])
        print
        return out
    
    def collateWrite(self,windList,gustList,outPath):
        ''''''
        print 'Length of wind List:',len(windList) 
        print 'Length of gust List:',len(gustList)
        out = []
        for i,item in enumerate(windList):
            if item[:2] == gustList[i][:2]:
                date = item[0]
                hrmin = item[1] 
                secs = item[2]
                wind = item[3] 
                gust = gustList[i][3]   
                out.append('{}  {}  {}  {:4s}  {:4s}\n'.format(date,hrmin,secs,wind.rjust(4),gust.rjust(4)))
            else:
                print 'Line {} the date or times do not match:'.format(i+2,),item[:2],gustList[i][:2]                
        writefile(out,outPath)

class extractToFile:
    '''Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/awos/1min.php.  Collates them into one file with 
    the format: day (starting from first datapoint), time, wind(m/s), gust(m/s)
    '''
    
    def __init__(self):
        return 


exToFile = extractToFile()
windPath = 'C:\\Users\\owner\\Downloads\\BooneWind1995-2011minute.txt'
gustPath = 'C:\\Users\\owner\\Downloads\\BooneGust1995-2011minute.txt'
outPath = 'C:\\Users\\owner\\Downloads\\BooneWindGust1995-2011minute.dat'
# windPath = 'C:\\Users\\owner\\Downloads\\wind10days.txt'
# gustPath = 'C:\\Users\\owner\\Downloads\\gust10days.txt'
# outPath = 'C:\\Users\\owner\\Downloads\\out10days.dat'

windList = exToFile.readVar(windPath)
gustList = exToFile.readVar(gustPath)
exToFile.collateWrite(windList,gustList,outPath)



print 'Done'
        
#         'D:\\Winch launch physics\\results\\testSpy'

