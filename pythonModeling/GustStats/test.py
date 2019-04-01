import json
import time,os,sys
import datetime
# Python 2 and 3: alternative 4
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

for i in range(10):
    print(i)
    
#     saveDir = 'C:\\Users\\owner\\Downloads\\'
saveDir = 'F:\\gustsDataRaw\\'
# os.chdir(saveDir)
# startts = datetime.datetime(2000, 1, 1)
# endts = datetime.datetime(2000, 1, 31)
# #     endts = datetime.datetime(2017, 1, 2)
# service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&"
# #     service = SERVICE + "data=sknt&data=drct&data=gust&tz=Etc/UTC&format=comma&latlon=no&"
# service += startts.strftime('year1=%Y&month1=%m&day1=%d&')
# service += endts.strftime('year2=%Y&month2=%m&day2=%d&')
states = ['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','HI','IA','ID','IL','IN','KS','KY','LA','MA','MD','ME',
          'MI','MN','MO','MS','MT','NC','ND','NE','NH','NJ','NM','NV','NY','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VA','VT',
          'WA','WI','WV','WY']
#     states = ['PA']
#     states = ['IA']
#     states = ['UT']
for state in states:
    print(state; sys.stdout.flush())
    stations = get_stations_from_network(state)
# stations = get_stations_from_filelist("mystations.txt")
    for station in stations:
        uri = '%s&station=%s' % (service, station)
        print('Downloading: %s' % (station, ))
        data = download_data(uri)
        outfn = '%s_%s_%s_%s.txt' % (state,station, startts.strftime("%Y%m%d%H%M"),
                                  endts.strftime("%Y%m%d%H%M"))
        out = open(outfn, 'w')
        out.write(data)
        out.close()
print("Done")
