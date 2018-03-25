# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 07:39:05 2018

@author: owner
"""
import os,sys
from matplotlib.pyplot import figure,plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim,legend,title,grid
from numpy import array
class plots:
    
    def __init__(self,path):
        self.i = 0 #counter for plots, so each variable has a different color
        self.path = path
        self.linew = 3.0
        self.colorsList = ['palevioletred', u'#fc4f30', u'#6d904f','darkorange', 'darkviolet', u'#8b8b8b',
        u'#348ABD', u'#e5ae38', u'#A60628', u'#7A68A6', 'mediumaquamarine', u'#D55E00', 'violet',
        u'#CC79A7',  u'#0072B2', u'#30a2da',u'#009E73','peru','slateblue'] # u'#F0E442',u'#467821','slateblue'      u'#56B4E9',
        return 
    
    def xy(self,xs,ys,xlbl,ylbl,legendLabels,titlestr):
        '''To allow different time (x) arrays, we require the xs to be a list'''
        if len(xs)<len(ys): #duplicate first x to every spot in list of correct length
            xs = [xs[0] for y in ys] 
        if len(legendLabels) != len(ys):
            sys.exit('Stop: number of legend labels and curves do not match for curve {}'.format(titlestr))
        figure(figsize=(20, 10))
        for iy,y in enumerate(ys):
            if len(xs[iy]) != len(y):
                sys.exit('Stop: curve {} on plot {} has x with dimensions different from y.'.format(iy+1,titlestr))            
            plot(xs[iy],y,'o-',color=self.colorsList[self.i],linewidth=self.linew,label=legendLabels[iy])
            self.i += 1
            if self.i > len(self.colorsList)-1:
                self.i = 0
        xlabel(xlbl)
        ylabel(ylbl)
        grid(color='k', linestyle='--', linewidth=1)
        legend(loc ='lower right',framealpha=0.5)
#        legend(loc ='upper left')
        ymin = min([min(y) for y in ys]); ymax = 1.1*max([max(y) for y in ys]);
        ylim([ymin,ymax])
        title(titlestr)
        savefig('{}{}{}.pdf'.format(self.path,os.sep,titlestr))
#        show(block = False)
        show()
        
Aasw27A =  [116]; Wasw27 = [70];
Aask21 =  [x]; Wask21 = [x];

Adg1000 =  [100]; Wasw27 = [81];
Adg505 =  [103]; Wdg505 = [76];
Aduo = [x]; Wduo = [x];
Apik20D =  [x]; Wpik20D = [x];

path = 'D:\\Winch launch physics\\maneuvering\\'  #for saving plots  
if not os.path.exists(path): os.mkdir(path)      
plts = plots(path)
plts.xy([Tinside,Toutside,Toutside2],[Vinside,Voutside,Voutside2],' Tension (lbs)','Signal (mV)',['Inside','Outside','Outside run 2'],'Amplfied tension signal 10 Mar 2018')

print 'Done'
        
        the ASW-27 has a Va = 116 kts and Vw = 70 kts.  

DG1000 has Va = 100 kts and Vw = 81 kts.  

DG505 has Va 103 kts and Vw = 76 kts. 