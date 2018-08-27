# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu

in one analysisDir, combines the events for wind and gusts into allWindEvents.dat and allGustEvents.dat, and makes plots

Read a task from a file tasks.dat (see taskGenerator.py), defined by state and parameters. 
0.  If a file read fails, wait a 0.5 sec...repeat x times or fail 
1.  Read the list of tasks. Start the first task. Delete the task from the list, and rewrite the list. 
2.  Task completion or failure is documented in each analysis folder (so there will be one tasks.dat file and multiple done.dat (which lists stations analyzed) or failed.dat
3.  Don't do the task if it is in done.dat
4.  done.dat is written to for each state analyzed for the same parameters (same analysis folder).  So we have to allow retries then.  

python D:\\WinchLaunchPhysics\\winchrep\\pythonModeling\\GustStats\\analysisTaskLooper.py

'''
import os,sys,time
from numpy import log10,savetxt,loadtxt,zeros,int32,sum
from matplotlib.pyplot import close
from analysisRoutines import correlate,plots,readfile,writefile

### Main script ###  
analysisDir = 'I:\\gustsData\\\\analysis30,5,15' 

params = analysisDir.split('analysis')[1].split(',')
t1 = int(params[0]); t2 = int(params[1]); vPastCap = int(params[2])
dt = 5 #min; the data sampling period
probFloor = 1e-9
vmax = 100 #kts
close('all') 
corr = correlate(probFloor) #instance
# get list of state files with analyzed events #  
windWindEventsPaths = []
gustGustEventsPaths = []
windGustEventsPaths = []
dirList = os.listdir(analysisDir)
for item in dirList:
    if 'windwindEvents.dat' in item:
        windWindEventsPaths.append('{}//{}'.format(analysisDir,item))
    elif 'gustgustEvents.dat' in item:
        gustGustEventsPaths.append('{}//{}'.format(analysisDir,item))
    elif 'windgustEvents.dat' in item:
        windGustEventsPaths.append('{}//{}'.format(analysisDir,item))
nFiles = len(windWindEventsPaths)
#collate events from state files
nWindWindEvents = zeros((vmax+1,vmax+1),dtype = int32)
nGustGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
nWindGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
for i,filename in enumerate(windWindEventsPaths):
    newwindWindEvents = loadtxt(filename, dtype=int32)
    newgustGustEvents = loadtxt(gustGustEventsPaths[i], dtype=int32)
    newwindGustEvents = loadtxt(windGustEventsPaths[i], dtype=int32)
    nWindWindEvents += newwindWindEvents
    nGustGustEvents += newgustGustEvents 
    nWindGustEvents += newwindGustEvents 
nWindTot = sum(nWindWindEvents)
print 'Total wind events:',nWindTot   
nGustTot = sum(nGustGustEvents)
print 'Total gust events:',nGustTot
nwindGustTot = sum(nWindGustEvents)
print 'Total wind-gust events:',nwindGustTot
savetxt('{}//allWindWindEvents.dat'.format(analysisDir),nWindWindEvents,fmt='%10d')
savetxt('{}//allGustGustEvents.dat'.format(analysisDir),nGustGustEvents,fmt='%10d')
savetxt('{}//allWindGustEvents.dat'.format(analysisDir),nWindGustEvents,fmt='%10d')
nWindWindEventsDispl = nWindWindEvents
nWindWindEventsDispl[0,0] = 0
nGustGustEventsDispl = nGustGustEvents
nGustGustEventsDispl[0,0] = 0
nWindGustEventsDispl = nWindGustEvents
nWindGustEventsDispl[0,0] = 0
### correlated probability ###
print 'probabilities:' 
#note if an (i,j) element in the end has no events in it, it's given probability of probFloor, a display floor used when calculating logs
probNextWindWindGTE,rowCountw = corr.probNextGTE(nWindWindEvents,vmax) 
probNextGustGustGTE,rowCountg = corr.probNextGTE(nGustGustEvents,vmax)
probNextWindGustGTE,rowCountg = corr.probNextGTE(nWindGustEvents,vmax)
### combined correlated probability 
# Add initial max velocities from 0 to 15
rows = range(vPastCap+1)
probNextWindWindGETcomb,rowsCount = corr.probNextGTECombine(nWindWindEvents,rows,vmax)
probNextGustGustGETcomb,rowsCount = corr.probNextGTECombine(nGustGustEvents,rows,vmax)
probNextWindGustGETcomb,rowsCount = corr.probNextGTECombine(nWindGustEvents,rows,vmax)
# Independent probability...simply sum over all intial states
rows = range(vmax+1)
probNextWindInd,rowsCount = corr.probNextGTECombine(nWindWindEvents,rows,vmax)
probNextGustInd,rowsCount = corr.probNextGTECombine(nGustGustEvents,rows,vmax)
### plots ###
pl = plots(analysisDir) #instance
## nEvents
#for display, log(1+arr) takes the log of 1+array to give contrast and avoid log(0)
titleStr='{} Wind Correlation (log N+1 events)'.format('All')
xlbl = 'Next {}-min max speed (kts)'.format(t2)
ylbl = 'last {}-min max speed (kts)'.format(t1)     
pl.matColor(False,log10(nWindWindEventsDispl+1),xlbl,ylbl,titleStr)
titleStr='{} Gust Correlation (log(N+1) events)'.format('All')
pl.matColor(False,log10(nGustGustEventsDispl+1),xlbl,ylbl,titleStr)
titleStr='{} Wind-Gust Correlation (log(N+1) events)'.format('All')
pl.matColor(False,log10(nWindGustEventsDispl+1),xlbl,ylbl,titleStr)
## Probabilities
titleStr='{} Wind probability (log10)'.format('All')
pl.matColor(False,log10(probNextWindWindGTE),xlbl,ylbl,titleStr)
titleStr='{} Gust probability (log10)'.format('All')
pl.matColor(False,log10(probNextGustGustGTE),xlbl,ylbl,titleStr)
titleStr='{} Wind-Gust probability (log10)'.format('All')
pl.matColor(False,log10(probNextWindGustGTE),xlbl,ylbl,titleStr)
# conditional
xs = [range(vmax+1)]*5
ys = [log10(probNextWindInd),log10(probNextWindWindGETcomb),log10(probNextGustInd),log10(probNextGustGustGETcomb),log10(probNextWindGustGETcomb)]
titlestr = '{} Probabilities of wind max in next {} minutes'.format('All',t2)
xlbl = 'Future wind speed (kts)'
ylbl = 'Probability (log10)'
strConditionalW = 'after {} min with winds <= {} kts'.format(t1,vPastCap)
strConditionalG = 'after {} min with gusts <= {} kts'.format(t1,vPastCap)
legendLabels = ['Wind: independent','Wind: {}'.format(strConditionalW),'Gust: independent', 'Gust: {}'.format(strConditionalG),'Gust: {}'.format(strConditionalW)  ]
#don't set hold to True, because it will stop looping. 
pl.xy(True,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None)
print 'Done'