"""
Modified example script that scrapes data from the IEM ASOS download service

D:\\WinchLaunchPhysics\\winchrep\\pythonModeling\\GustStats\\readASOSserver3.py

"""
# from __future__ import print_function
import json
import time,os,sys
import datetime
# Python 2 and 3: alternative 4
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

# Number of attempts to download data
MAX_ATTEMPTS = 6
# HTTPS here can be problematic for installs that don't have Lets Encrypt CA
SERVICE = "http://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"

def download_data(uri):
    """Fetch the data from the IEM
    The IEM download service has some protections in place to keep the number
    of inbound requests in check.  This function implements an exponential
    backoff to keep individual downloads from erroring.
    Args:
      uri (string): URL to fetch
    Returns:
      string data
    """
    attempt = 0
    while attempt < MAX_ATTEMPTS:
        try:
            data = urlopen(uri, timeout=300).read().decode('utf-8')
            if data is not None and not data.startswith('ERROR'):
                return data
        except Exception as exp:
            print("download_data({}) failed with{}".format(uri, exp))
            time.sleep(5)
        attempt += 1

    print("Exhausted attempts to download, returning empty data" )
    return ""


def get_stations_from_filelist(filename):
    """Build a listing of stations from a simple file listing the stations.
    The file should simply have one station per line.
    """
    stations = []
    for line in open(filename):
        stations.append(line.strip())
    return stations

def get_stations_from_network(network):
    """Build a station list"""
    stations = []
    # IEM quirk to have Iowa AWOS sites in its own labeled network
#     networks = ['AWOS']
    # Get metadata
    uri = ("https://mesonet.agron.iastate.edu/"
               "geojson/network/{}_ASOS.geojson".format(network))
    data = urlopen(uri)
    jdict = json.load(data)
    for site in jdict['features']:
        stations.append(site['properties']['sid'])
    return stations

def readStationsDone(outPath):
    '''Read all stations done from the filenames in the output directory'''
    dirList = os.listdir(outPath)
    stationsDone = []
    for dirTxt in dirList:
        if dirTxt[0].isalpha() and '_' in dirTxt and os.path.getsize(dirTxt) > 0:
            stationsDone.append(dirTxt.split('_')[1])
    return stationsDone

#########  Main script #############
#     outPath = 'C:\\Users\\owner\\Downloads\\' 
outPath = 'H:\gustsData'
if not os.path.exists(outPath): os.mkdir(outPath)
os.chdir(outPath)
startts = datetime.datetime(1960, 1, 1)
endts = datetime.datetime(2018, 12, 31)
#     endts = datetime.datetime(2017, 1, 2)
service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&"
#     service = SERVICE + "data=sknt&data=drct&data=gust&tz=Etc/UTC&format=comma&latlon=no&"
service += startts.strftime('year1=%Y&month1=%m&day1=%d&')
service += endts.strftime('year2=%Y&month2=%m&day2=%d&')
states = ['DE','FL','GA','HI','IA','ID','IL','IN','KS','KY','LA','MA','MD','ME',
          'MI','MN','MO','MS','MT','NC','ND','NE','NH','NJ','NM','NV','NY','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VA','VT',
          'WA','WI','WV','WY']
# states = ['PA']
# states = ['IA']
# states = ['UT','WY']
# states = ['WA','WI','WV','WY']
stationsDone = readStationsDone(outPath)
for state in states:
    print('==========',state,'=========='); sys.stdout.flush()
    stations = get_stations_from_network(state)
# stations = get_stations_from_filelist("mystations.txt")
    for station in stations:
        if station not in stationsDone:
            uri = '%s&station=%s' % (service, station)
            print('Downloading:', state, station); sys.stdout.flush()
            data = download_data(uri)
            if len(data)> 100:
                outfn = '%s_%s_%s_%s.txt' % (state,station, startts.strftime("%Y%m%d%H%M"),
                                          endts.strftime("%Y%m%d%H%M"))
                out = open(outfn, 'w')
                out.write(data)
                out.close()
print("Done")
