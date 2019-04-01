# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com 
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
    
    def readVar(self,inPath,type):
        '''Reads date, time and wind *or* gust speed into an array
        Format BNW BOONE MUNI 1995-01-01 01:34 15'''
        dtFormat = "%Y-%m-%d %H:%M"
        out = []
        lines = readfile(inPath)
        print('Reading {}s from raw data'.format(type))
        for i,line in enumerate(lines[1:]):
            if mod(i,10000)==0:
                print('{} '.format(i),)
            if mod(i,100000)==0:
                print()
            info = line.split()
            if len(info) >=5:
#             if i>5120000: print('info',i,info)
                date = info[3]
                hrmin = info[4]
                if len(info) >= 6:
                    var = info[5]
                else:
                    var = ''
                    print('Line {} in type {} has a length of only {}.'.format(i,type,len(info)),info)
                d = datetime.strptime(date+' '+hrmin, dtFormat)
                secs = int(time.mktime(d.timetuple()))
                out.append([date,hrmin,secs,var])
            else:
                print('Line {} in type {} does not have complete time and date:'.format(i,type),info)
                sys.exit('Stop!')
        print()
        return out
    
    def collateData(self,windList,gustList,outPath=None):
        '''Writes data into structured array, and writes to a file if outPath exists'''
        print('Length of wind List:',len(windList) )
        print('Length of gust List:',len(gustList))
        if len(windList) != len(gustList):
            sys.exit('Stop: the two lengths are not equal!')
        data = zeros(len(windList),dtype = [('year', int),('month', int),('day', int),('hour', int),('min', int),('totsecs', int32),('wind',int),('gust', int)])
        out = []
        zeroSecs =  int(windList[0][2])
        idata = 0
        for i,item in enumerate(windList):
            if item[:2] == gustList[i][:2] and len(item[:2])>0 and len(gustList[i][:2]):
                date = item[0]
                hrmin = item[1] 
                secs = int(item[2]) - zeroSecs
                wind = item[3] 
                gust = gustList[i][3]
                dateInfo = date.split('-')
                data[idata]['year'] = dateInfo[0]
                data[idata]['month'] = dateInfo[1]
                data[idata]['day'] = dateInfo[2]
                timeInfo = hrmin.split(':')
                data[idata]['hour'] = timeInfo[0]
                data[idata]['min'] = timeInfo[1]
                data[idata]['totsecs'] = secs
                data[idata]['wind'] = wind
                data[idata]['gust'] = gust
                if not outPath is None:   
                    out.append('{}  {}  {}  {:4s}  {:4s}\n'.format(date,hrmin,secs,wind.rjust(4),gust.rjust(4)))  
            else:
                print('Line {} the date or times do not match:'.format(i+2,),item[:2],gustList[i][:2]  )
            idata += 1              
        if not outPath is None: writefile(out,outPath)
        
    def readCollatedData(self,inPath):
        '''Reads my saved data into structured array'''
        print('Reading from collated data file',inPath)
        if len(windList) != len(gustList):
            sys.exit('Stop: the two lengths are not equal!')
        data = zeros(len(windList),dtype = [('year', int),('month', int),('day', int),('hour', int),('min', int),('totsecs', int32),('wind',int),('gust', int)])
        out = []
        zeroSecs =  int(windList[0][2])
        idata = 0
        for i,item in enumerate(windList):
            if item[:2] == gustList[i][:2] and len(item[:2])>0 and len(gustList[i][:2]):
                date = item[0]
                hrmin = item[1] 
                secs = int(item[2]) - zeroSecs
                wind = item[3] 
                gust = gustList[i][3]
                dateInfo = date.split('-')
                data[idata]['year'] = dateInfo[0]
                data[idata]['month'] = dateInfo[1]
                data[idata]['day'] = dateInfo[2]
                timeInfo = hrmin.split(':')
                data[idata]['hour'] = timeInfo[0]
                data[idata]['min'] = timeInfo[1]
                data[idata]['totsecs'] = secs
                data[idata]['wind'] = wind
                data[idata]['gust'] = gust
                if not outPath is None:   
                    out.append('{}  {}  {}  {:4s}  {:4s}\n'.format(date,hrmin,secs,wind.rjust(4),gust.rjust(4)))  
            else:
                print('Line {} the date or times do not match:'.format(i+2,),item[:2],gustList[i][:2]  )
            idata += 1              
        if not outPath is None: writefile(out,outPath)


class correlate:
    '''Reads data from files downloaded from https://mesonet.agron.iastate.edu/request/awos/1min.php.  Collates them into one file with 
    the format: day (starting from first datapoint), time, wind(m/s), gust(m/s)
    '''
    
    def __init__(self):
        return 


rdData = readData()
windPath = 'C:\\Users\\owner\\Downloads\\BooneWind1995-2011minute.txt'
gustPath = 'C:\\Users\\owner\\Downloads\\BooneGust1995-2011minute.txt'
outPath = 'C:\\Users\\owner\\Downloads\\BooneWindGust1995-2011minute.dat'
# windPath = 'C:\\Users\\owner\\Downloads\\wind10days.txt'
# gustPath = 'C:\\Users\\owner\\Downloads\\gust10days.txt'
# outPath = 'C:\\Users\\owner\\Downloads\\out10days.dat'

windList = rdData.readVar(windPath,'wind')
gustList = rdData.readVar(gustPath, 'gust')
rdData.collateData(windList,gustList,outPath)

print('Done')
        
#         'D:\\Winch launch physics\\results\\testSpy'

