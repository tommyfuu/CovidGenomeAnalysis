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
    memoizedDict = {}
    # index = 0
    for fileName in alignScores:
        tempdf = pd.read_csv(fileName)
        outfileName = "relativeDistAll" + fileName
        # initialize
        geoLocation = []
        relativeDistance = []
        ascensionNumL = []
        # get relative distance
        # for index in range(1000):
        for index in range(len(tempdf["Location"])):
            print(index)

            location = tempdf["Location"][index]
            ascensionNumL.append(tempdf["AscensionNum"][index])
            # check memo dictionary to see if we've already done this loc
            # note that if geolocation is none we just ignore it
            if not location in memoizedDict.keys():
                currentGeo = findGeocode(location)
                # save relative distance
                if currentGeo != None:
                    latitudeLongitude = (
                        currentGeo.latitude, currentGeo.longitude)
                    if type(currentGeo) != float:
                        memoizedDict[location] = (distance(
                            originPoint, latitudeLongitude), currentGeo)
                else:
                    memoizedDict[location] = (0, currentGeo)
            relativeDistance.append(memoizedDict[location][0])
            geoLocation.append(memoizedDict[location][1])

        # Now add the column to dataframe
        currentOutputDF = pd.DataFrame(
            {"geoLocation": geoLocation,
             "relativeDistance": relativeDistance,
             "AscensionNum": ascensionNumL
             })
        # And output the df to a new csv file.
        currentOutputDF.to_csv(outfileName)

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
