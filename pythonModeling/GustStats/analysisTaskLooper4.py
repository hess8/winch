# -*- coding: utf-8 -*-
''' Bret Hess, bret.hess@gmail.com 

Read a task from a file tasks.dat (see taskGenerator.py), defined by state and parameters. 
0.  If a file read fails, wait...repeat x times or fail 
1.  Read the list of tasks. Start the first task. Delete the task from the list, and rewrite the list. 
2.  Task completion or failure is documented in each analysis folder (so there will be one tasks.dat file and multiple done.dat (which lists stations analyzed) or failed.dat
3.  Don't do the task if it is in done.dat
4.  done.dat is written for each state analyzed for the same parameters (same analysis folder).  So we have to allow retries then.  

python D:\\WinchLaunchPhysics\\winchrep\\pythonModeling\\GustStats\\analysisTaskLooper2.py

If we need more sophisticated file locking: https://github.com/dmfrey/FileLock/blob/master/filelock/filelock.py
'''

import os,sys,time
from numpy import log10,savetxt,loadtxt,zeros,int32,sum
from matplotlib.pyplot import close
from analysisRoutines import correlate,plots,readAnalysisTask,readfile,writefile,readData,writeDone
from copy import deepcopy

### Main script ###  
OK = True        
# inPath = 'F:\\temp\\temp2'
# inPath = 'F:\\temp\\'
inPath ='H:\\gustsData\\'
outPath = 'F:\\gustsAnalysis\\'
if not os.path.exists(outPath): os.mkdir(outPath)
#Analysis variables:
dt = 5 #min; the data sampling period
probFloor = 1e-9
vmax = 100 #kts
#end analysis variables:

loop = True
while loop:
    #find a state with task not done
    close('all') 
    taskOut = readAnalysisTask(outPath)
    if taskOut == 'ReadFailed':
        print('Failed to read task')
        time.sleep(1)
        OK = False
    elif taskOut == 'NoTasks':
        sys.exit('Done: no tasks remaining')
    else:
        state = taskOut[0]; t1 = int(taskOut[1]); t2 = int(taskOut[2]); vPastCap = int(taskOut[3]);
    #paths
    analysisDir = '{}\\analysis{},{},{}'.format(outPath,t1,t2,vPastCap) 
    if not os.path.exists(analysisDir): os.mkdir(analysisDir)
    statesDonePath = '{}\\statesDone.dat'.format(analysisDir)
    statesDone = [] #states that are done with all stations and plots
    if not os.path.exists(statesDonePath): 
        writefile([],statesDonePath) #write empty file
    else:
        statesDone = readfile(statesDonePath) #read stations done with this analysis
    #stations from any state that are done with analysis
    donePath = '{}\\done.dat'.format(analysisDir)
    done = [] 
    if not os.path.exists(donePath): 
        writefile([],donePath) #write empty file
    else:
        done = readfile(donePath)
    #start task
    if state not in statesDone:
        print('Running task {}'.format(taskOut))
        # verbose = True
        # redoPast = True
        verbose = False
        redoPast = False
        rdData = readData() #instance
        corr = correlate(probFloor) #instance
        nWindWindFile = '{}\\{}_windwindEvents.dat'.format(analysisDir,state)
        nGustGustFile = '{}\\{}_gustgustEvents.dat'.format(analysisDir,state)
        nWindGustFile = '{}\\{}_windgustEvents.dat'.format(analysisDir,state)
        if os.path.exists(nWindWindFile) and os.path.exists(nGustGustFile) and os.path.exists(nWindGustFile):
            nWindWindEvents = loadtxt(nWindWindFile, dtype=int32)
            nGustGustEvents = loadtxt(nGustGustFile, dtype=int32)
            nWindGustEvents = loadtxt(nWindGustFile, dtype=int32)
        else:
            nWindWindEvents = zeros((vmax+1,vmax+1),dtype = int32)
            nGustGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
            nWindGustEvents = zeros((vmax+1,vmax+1),dtype = int32)
        #----
        statePaths = rdData.readStatePaths(inPath,state)
        if len(statePaths) == 0:
            print('There are no data files matching state {}'.format(state))
        else:
            for stationPath in statePaths:
                station = stationPath.split('_')[1]
                if station not in done or redoPast:
                    print(station,stationPath); sys.stdout.flush()
                    nData,data = rdData.readStationData('{}\\{}'.format(inPath,stationPath),verbose)
                    newWindEvents, newGustEvents, newWindGustEvents = corr.nEventsCount(t1,t2,dt,nData,data,vmax,verbose)
                    nWindWindEvents += newWindEvents
                    nGustGustEvents += newGustEvents
                    nWindGustEvents += newWindGustEvents
                    #write state events counter to file
                    savetxt(nWindWindFile,nWindWindEvents,fmt='%10d')
                    savetxt(nGustGustFile,nGustGustEvents,fmt='%10d')
                    savetxt(nWindGustFile,nWindGustEvents,fmt='%10d')
                    #record station as analyzed
                    status = writeDone(donePath,'{}\n'.format(station))
                    if status == 'failed':
                        print('Failed to write to done.dat')
            ### correlations for this state ###
            print('Number of wind events: {}'.format(sum(nWindWindEvents)))
            print('Number of gust events: {}'.format(sum(nGustGustEvents)))
        #     for i in range(vmax+1):
        #         for j in range(vmax+1):
        #             print(i,j,'\t',nWindWindEvents[i,j],'\t',nGustGustEvents[i,j])
            if sum(nWindWindEvents)>0:
                nWindWindEventsDispl = deepcopy(nWindWindEvents)
                nWindWindEventsDispl[0,0] = 0
                nGustGustEventsDispl = deepcopy(nGustGustEvents)
                nGustGustEventsDispl[0,0] = 0
                nWindGustEventsDispl = deepcopy(nWindGustEvents)
                nWindGustEventsDispl[0,0] = 0
                ### correlated probability ###  
#                 print('probabilities:' )
                #note if an (i,j) element in the end has no events in it, it's given probability of probFloor, a display floor useful when calculating logs
#                 probNextWindWindGTE,rowCountw = corr.probNextGTE(nWindWindEvents,vmax) 
#                 probNextGustGustGTE,rowCountg = corr.probNextGTE(nGustGustEvents,vmax)
#                 probNextWindGustGTE,rowCountg = corr.probNextGTE(nWindGustEvents,vmax)
                ### combined correlated probability 
                # Add initial max velocities events from 0 to 15
                rows = range(vPastCap+1)
                probNextWindWindGTEcomb,rowsCount = corr.probNextGTECombine(nWindWindEvents,rows,vmax)
                probNextGustGustGTEcomb,rowsCount = corr.probNextGTECombine(nGustGustEvents,rows,vmax)
                probNextWindGustGTEcomb,rowsCount = corr.probNextGTECombine(nWindGustEvents,rows,vmax)
                # Independent probability...simply sum over all initial states
                rows = range(vmax+1)
                probNextWindInd,rowsCount = corr.probNextGTECombine(nWindWindEvents,rows,vmax)
                probNextGustInd,rowsCount = corr.probNextGTECombine(nGustGustEvents,rows,vmax)
                ### plots ###
                pl = plots(analysisDir) #instance
                ## nEvents
                #for display, log(1+arr) takes the log of 1+array to give contrast and avoid log(0)
                titleStr='{} Wind Correlation (log N+1 events)'.format(state)
                xlbl = 'Next {}-min max speed (kts)'.format(t2)
                ylbl = 'last {}-min max speed (kts)'.format(t1)     
                pl.matColor(False,log10(nWindWindEventsDispl+1),xlbl,ylbl,titleStr)
                titleStr='{} Gust Correlation (log(N+1) events)'.format(state)
                pl.matColor(False,log10(nGustGustEventsDispl+1),xlbl,ylbl,titleStr)
                titleStr='{} Wind-Gust Correlation (log(N+1) events)'.format(state)
                pl.matColor(False,log10(nWindGustEventsDispl+1),xlbl,ylbl,titleStr)
                ## Probabilities
#                 titleStr='{} Wind probability (log10)'.format(state)
#                 pl.matColor(False,log10(probNextWindWindGTE),xlbl,ylbl,titleStr)
#                 titleStr='{} Gust probability (log10)'.format(state)
#                 pl.matColor(False,log10(probNextGustGustGTE),xlbl,ylbl,titleStr)
#                 titleStr='{} Wind-Gust probability (log10)'.format(state)
#                 pl.matColor(False,log10(probNextWindGustGTE),xlbl,ylbl,titleStr)
                # conditional
                xs = [range(vmax+1)]*5
                ys = [log10(probNextWindInd),log10(probNextWindWindGTEcomb),log10(probNextGustInd),log10(probNextGustGustGTEcomb),log10(probNextWindGustGTEcomb)]
                titlestr = '{} Probabilities of wind max in next {} minutes'.format(state,t2)
                xlbl = 'Future wind speed (kts)'
                ylbl = 'Probability (log10)'
                strConditionalW = 'after {} min with winds <= {} kts'.format(t1,vPastCap)
                strConditionalG = 'after {} min with gusts <= {} kts'.format(t1,vPastCap)
                legendLabels = ['Wind: independent','Wind: {}'.format(strConditionalW),'Gust: independent', 'Gust: {}'.format(strConditionalG),'Gust: {}'.format(strConditionalW)  ]
                #don't set hold to True, because it will stop looping. 
                pl.xy(False,xs,ys,xlbl,ylbl,legendLabels,titlestr,xmin=None,xmax=None)
        #record state as analyzed
        status = writeDone(statesDonePath,'{}\n'.format(state))
        if status == 'failed':
            print('Failed to write to statesDone.dat'        )