


# Import pandas package   
import pandas as pd  
import numpy as np 
from geopy.exc import GeocoderTimedOut 
from geopy.geocoders import Nominatim 

# Store a list of file names
alignScores = ["ORF1aAlignmentScore.csv", "ORF1bAlignmentScore.csv", "ORF3aAlignmentScore.csv"]

def readLocation():
    for fileName in alignScores:
        tempdf = pd.read_csv(fileName)
        outfileName = "updated" + fileName
        # List to store the latitude-longitude tuple 
        geoLocation = [] 
        for i in (tempdf["Location"]):
            if findGeocode(i) != None:
                Loc = findGeocode(i)
                geoLocation.append((Loc.latitude, Loc.longitude))
            else:
                geoLocation.append(np.nan)
        
        # Now add the column to dataframe
        tempdf["geoLocation"] = geoLocation

        # And output the df to a new csv file.
        tempdf.to_csv(outfileName)

## Funcition to find the coordinate given a city name
## Reference:
## https://www.geeksforgeeks.org/how-to-find-longitude-and-latitude-for-a-list-of-regions-or-country-using-python/
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
