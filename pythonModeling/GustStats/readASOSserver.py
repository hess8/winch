"""
Example script that scrapes data from the IEM ASOS download service

D:\\WinchLaunchPhysics\\winchrep\\pythonModeling\\GustStats\\readASOSserver.py

"""
from __future__ import print_function
import json
import time,os
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
            print("download_data(%s) failed with %s" % (uri, exp))
            time.sleep(5)
        attempt += 1

    print("Exhausted attempts to download, returning empty data")
    return ""


def get_stations_from_filelist(filename):
    """Build a listing of stations from a simple file listing the stations.
    The file should simply have one station per line.
    """
    stations = []
    for line in open(filename):
        stations.append(line.strip())
    return stations


def get_stations_from_networks():
    """Build a station list by using a bunch of IEM networks."""
    stations = []
    stateList = []
#     states = """AK AL AR AZ CA CO CT DE FL GA HI IA ID IL IN KS KY LA MA MD ME
#      MI MN MO MS MT NC ND NE NH NJ NM NV NY OH OK OR PA RI SC SD TN TX UT VA VT
#      WA WI WV WY"""
#     states = """PA"""
#     states = """IA"""
    states = """UT"""
    # IEM quirk to have Iowa AWOS sites in its own labeled network
#     networks = ['AWOS']
    networks = []
    for state in states.split():
        print((state,),)
        networks.append("%s_ASOS" % (state,))

    for network in networks:
        # Get metadata
        uri = ("https://mesonet.agron.iastate.edu/"
               "geojson/network/%s.geojson") % (network,)
        data = urlopen(uri)
        jdict = json.load(data)
        for site in jdict['features']:
            stations.append(site['properties']['sid'])
            stateList.append(network.split('_')[0])
    return stations,stateList

def main():
    """Our main method"""
    # timestamps in UTC to request data for
#     saveDir = 'C:\\Users\\owner\\Downloads\\'
    saveDir = 'F:\\gustsDataRaw\\'
    os.chdir(saveDir)
    startts = datetime.datetime(2000, 1, 1)
    endts = datetime.datetime(2017, 12, 31)
#     endts = datetime.datetime(2017, 1, 2)
    service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&"
#     service = SERVICE + "data=sknt&data=drct&data=gust&tz=Etc/UTC&format=comma&latlon=no&"

    service += startts.strftime('year1=%Y&month1=%m&day1=%d&')
    service += endts.strftime('year2=%Y&month2=%m&day2=%d&')

    # Two examples of how to specify a list of stations
    stations,stateList = get_stations_from_networks()
    # stations = get_stations_from_filelist("mystations.txt")
    for i, station in enumerate(stations):
        uri = '%s&station=%s' % (service, station)
        print('Downloading: %s' % (station, ))
        data = download_data(uri)
        outfn = '%s_%s_%s_%s.txt' % (stateList[i],station, startts.strftime("%Y%m%d%H%M"),
                                  endts.strftime("%Y%m%d%H%M"))
        out = open(outfn, 'w')
        out.write(data)
        out.close()
    print("Done")

if __name__ == '__main__':
    main()