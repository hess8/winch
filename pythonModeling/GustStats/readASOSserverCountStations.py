"""
Just count the number of stations in the states chosen.

D:\\WinchLaunchPhysics\\winchrep\\pythonModeling\\GustStats\\readASOSserver2.py

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

#########  Main script #############
startts = datetime.datetime(1980, 1, 1)
endts = datetime.datetime(2017, 12, 31)
#     endts = datetime.datetime(2017, 1, 2)
states = ['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','HI','IA','ID','IL','IN','KS','KY','LA','MA','MD','ME',
          'MI','MN','MO','MS','MT','NC','ND','NE','NH','NJ','NM','NV','NY','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VA','VT',
          'WA','WI','WV','WY']
# states = ['PA']
# states = ['IA']
# states = ['UT']
nstations = 0
for state in states:
    print('==========',state,'=========='); sys.stdout.flush()
    stations = get_stations_from_network(state)
    print(len(stations))
    nstations += len(stations)
print('Total:', nstations)
print("Done")
