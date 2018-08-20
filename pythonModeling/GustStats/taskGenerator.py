''' Bret Hess, bret.hess@gmail.com or bret_hess@byu.edu
Generates the tasks from loops over parameters and writes to '''
def writefile(lines,filepath): #need to have \n's inserted already
    file1 = open(filepath,'w')
    file1.writelines(lines) 
    file1.close()
    return

class tasks:
    '''Generates the tasks and writes then to tasks.dat with the following format: first line: states.  Other lines t1,t2,vcap '''
    def __init__(self):
        return 

    def loopParams(self):
        return
    
    ### Main script ###
    inPath ='I:\\gustsDataRaw\\'
    outPath = inPath
    
    
    states = ['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','HI','IA','ID','IL','IN','KS','KY','LA','MA','MD','ME',
              'MI','MN','MO','MS','MT','NC','ND','NE','NH','NJ','NM','NV','NY','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VA','VT',
              'WA','WI','WV','WY']
    t1List = [5,30]
    t2List = [5]
    vcapList = [15] #kts
    lines = []
    for state in states:
        for t1 in t1List:
            for t2 in t2List:
                for vcap in vcapList:
                    lines.append('{:2s} {:8d} {:8d} {:8d}\n'.format(state,t1,t2,vcap))
    writefile(lines,'{}//tasks.dat'.format(outPath))
    print 'Done'