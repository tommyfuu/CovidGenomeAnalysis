import pandas as pd
import numpy as np
import csv
from geopy.exc import GeocoderTimedOut
from geopy.geocoders import Nominatim
from math import radians, cos, sin, asin, sqrt

# Store a list of file names
address1 = 'outputs1107/ORF1a0300000AlignmentScore.csv'
address2 = 'outputs1107/ORF1b0300000AlignmentScore.csv'
address3 = 'outputs1107/ORF3a0300000AlignmentScore.csv'
address4 = 'outputs1107/ORFN0300000AlignmentScore.csv'
address5 = 'outputs1107/ORFS0300000AlignmentScore.csv'
alignScores = [address1, address2, address3, address4, address5]

# Location of Wuhan, which we set as the origin point
originPoint = (30.5951051, 114.2999353)


def readLocation():
    arr = np.empty((0, 2), str)
    arr = np.append(arr,  np.array(
        [['Location', 'Distance from Wuhan']]), axis=0)
    print(arr)
    index = 0
    for fileName in alignScores:
        tempdf = pd.read_csv(fileName)
        outfileName = "updated" + fileName
        # List to store the latitude-longitude tuple
        geoLocation = []
        relativeDistance = []
        # for i in (tempdf["Location"]):
        for i in (tempdf["Location"]):
            print(index)
            index += 1
            arr = np.append(arr, np.array(
                [[i, distance(originPoint, findGeocode(i)[1])]]), axis=0)
            # row = np.array([i, distance(originPoint, findGeocode(i)[1])])
            # arr = np.append(arr, row, axis=0)
            # print(arr)
            # print(i)
            # print(distance(originPoint, findGeocode(i)[1]))
            np.savetxt("ReadingLocationOutput.csv",
                       arr, delimiter=',', fmt='%s')
            if findGeocode(i) != None:
                Loc = findGeocode(i)
                latitudeLongitude = (Loc.latitude, Loc.longitude)
                geoLocation.append(latitudeLongitude)
                if type(findGeocode(i)) != float:
                    relativeDistance.append(
                        distance(originPoint, latitudeLongitude))
                else:
                    relativeDistance.append(0)
            else:
                geoLocation.append(np.nan)
                relativeDistance.append(np.nan)

        # Now add the column to dataframe
        tempdf["geoLocation"] = geoLocation
        tempdf["relativeDistance"] = relativeDistance

        # And output the df to a new csv file.
        tempdf.to_csv(outfileName)

# Funcition to find the coordinate given a city name
# Reference:
# https://www.geeksforgeeks.org/how-to-find-longitude-and-latitude-for-a-list-of-regions-or-country-using-python/


def findGeocode(city):
    # try and catch is used to overcome
    # the exception thrown by geolocator
    # using geocodertimedout
    try:
        # Specify the user_agent as your
        # app name it should not be none
        geolocator = Nominatim(user_agent="your_app_name")
        return geolocator.geocode(city)

    except GeocoderTimedOut:
        return findGeocode(city)


def distance(latitudeLongitude1, latitudeLongitude2):
    lat1 = latitudeLongitude1[0]
    lon1 = latitudeLongitude1[1]
    lat2 = latitudeLongitude2[0]
    lon2 = latitudeLongitude2[1]

    # Convert from degrees to radians
    lon1 = radians(lon1)
    lon2 = radians(lon2)
    lat1 = radians(lat1)
    lat2 = radians(lat2)

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2

    c = 2 * asin(sqrt(a))

    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371

    # calculate the result
    return(c * r)
