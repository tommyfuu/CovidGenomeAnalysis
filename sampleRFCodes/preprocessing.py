"""
Author      : Tom Fu
Date        : 2020 June 2
FileName    : normaliseBINU.py (for Coarfa Lab)
Description : Transform and normalise input txt files such as EGS_pa
"""

import pandas as pd
import numpy as np
# from scipy.stats import zscore



def formXList(address):
    """
    Intro: Turn txt into a list form dataframe and invert it
    Para: address, a string address of a import file
    Output: a data frame with each row being a list of data for each
    sample/individual
    """
    df = pd.read_csv(address, delimiter='\t')
    Xdf = df.T
    X = Xdf.values.tolist()
    return X

                  
                  
def rowNormalise(row):
    """
    Intro: normalise a row of a dafaframe
    Para: row, a list
    Output: normalisedRow, a normalised list
    """
    normalisedRow = row
    featureName = normalisedRow.pop(0)
    maxRow = max(row) # the max number of the row
    minRow = min(row)
    rangeRow = maxRow - minRow

    # list comprehension to normalise each element
    normalisedRow = [(num-minRow)/rangeRow for num in normalisedRow]
    normalisedRow.insert(0, featureName)
    return normalisedRow


def BINUlabelDictsWGS(addressY):
    """
    Intro: turn a y address into three dicts where you input the name of the 
            variable to get the index, with which you can search for its 
            optimised ML parameters in either its binary or numeric dict. It
            also returns the 3 datasets matched with relevent datatypes
    Para: addressY, a string address of a import file for y
    Output: binaryLabelDict, a dict for binary variables
            numericLabelDict, a dict for numeric variables
            categoricLabelDict, a dict for categoric variables
            dfBinaryY, a pd dataframe containing all binary variable data for y
            dfCategoricY, a pd dataframe containing all categoric variable data for y
            dfNumericY, a pd dataframe containing all numeric variable data for y
    """
    # initialise
    binary_cols_Y = []
    catego_cols_Y = []
    numeri_cols_Y = []
    binaryLabelDict = {}
    numericLabelDict = {}
    categoricLabelDict = {}
    
    # recognise and pick out binary/categorical/numeric dataset
    binary_cols_Y = ['Depression_Remission', 'Anxiety_Remission', 'DepOutMod', 'AnxOutMod', 'SSRI',
                    'SNRI', 'ANTIDEPRESSANTOTHER', 'MOODSTABILIZER_OR_ANTICONVULSANT', 'ATYPICALANTIPSYCHOTIC', 
                    'BENZODIAZEPINE', 'ALTERNATIVE_OR_COMPLEMENTARY']
    catego_cols_Y = ['PHQDepressionCat', 'GAD7AnxietyCat', 'MSI']
    numeri_cols_Y = ['WHOATOBACCO_1', 'WHOAALCOHOL_1','WHOACANNIBIS_1','WHOACOCAINE_1','WHOAAMPHETAMINE_1',
                     'WHOASEDATIVE_1','WHOAOPIOID_1','WHOATOTAL_1','PHQDEPRESSION_1','GAD7ANXIETY_1',
                     'PHQSOMATIZATION_1','PHQDEPRESSION_D1','GAD7ANXIETY_D1','PHQSOMATIZATION_D1','SSRSIDEATIONMTH_1',
                     'SSRSIDEATIONMTH_D1','DERSTOTAL_1','DERSTOTAL_D1','AAQTOTAL_1','AAQTOTAL_D1','WHOFRAWSCORE_1',
                     'WHOFRAWSCORE_D1','WHODTOTAL_1','WHODTOTAL_D1','CSSRS_tot','CSSRS_F1','CSSRS_F2','PHQDEPRESSION_DIFF',
                     'PHQANXIETY_DIFF','PHQSOMATIZATION_DIFF','SSRSIDEATIONMTH_DIFF','DERSTOTAL_DIFF','AAQTOTAL_DIFF',
                     'WHOFRAWSCORE_DIFF','WHODTOTAL_DIFF']
    
    # get y data
    dfY = pd.read_excel(addressY)
    dfBinaryY = dfY[binary_cols_Y]
    dfCategoricY = dfY[catego_cols_Y]
    dfNumericY = dfY[numeri_cols_Y]
    
    # form labelDicts
    for binaryIndex in range(0, len(binary_cols_Y)):
        binaryLabelDict.update( {binary_cols_Y[binaryIndex] : binaryIndex})
    for categoIndex in range(0, len(catego_cols_Y)):
        categoricLabelDict.update( {catego_cols_Y[categoIndex] : categoIndex})
    for numericIndex in range(0, len(numeri_cols_Y)):
        numericLabelDict.update( {numeri_cols_Y[numericIndex] : numericIndex})
    
    return binaryLabelDict, categoricLabelDict, numericLabelDict, dfBinaryY, dfCategoricY, dfNumericY

def sampleDropper(X, currentY):
    """
    Intro: for a specific variable, drop samples in that y that is nan and make a respective copyY that doesn't have that sample
    Para: X, normalised X sample
          currentY, all y data belonging to the current variable
    Output: realX, realCurrentY, which are X and currentY with nan values dropped
    """
    nanIndices = [i for i,v in enumerate(currentY) if np.isnan(v)]
    realCurrentY = np.delete(currentY, nanIndices, axis=0)
    realCurrentX = np.delete(X, nanIndices, axis=0)
    return realCurrentX, realCurrentY

def XGetter(addressX):
    # get X
    dfList = formXList(addressX)
    featureList = dfList.pop(0)
    normalisedXList = [rowNormalise(row) for row in dfList]
    normalisedXdf = pd.DataFrame.from_records(normalisedXList)
    normalisedXdf = normalisedXdf.select_dtypes('number')
    # np arrayise
    X = np.array(normalisedXdf)
    return X, featureList

def sampleDropperAll(X, currentYType):
    nanIndices2D =np.argwhere(np.isnan(currentYType))
    nanIndices = [sample[0] for sample in nanIndices2D]
    nanIndices = list( dict.fromkeys(nanIndices) )
    realCurrentTypeY = np.delete(currentYType, nanIndices, axis=0)
    realCurrentX = np.delete(X, nanIndices, axis=0)
    return realCurrentX, realCurrentTypeY
    

################## overarching function ##################

def BINUnormalise(addressX, addressY):
    """
    Intro: normalise a dafaframe, a 2D list
    Para: addressX and addressY, two address corresponding to X and Y
    Output: normalised X, formatted y, and the respect label dicts
    """
    # get X
    X, allFeatureList = XGetter(addressX)
    # get Y and label dicts
    binaryLabelDict, categoricLabelDict, numericLabelDict, dfBinaryY, dfCategoricY, dfNumericY = BINUlabelDictsWGS(addressY)
    
    # for binary, change yes to 1 no to 0
    # mapping = {"Yes": 1, "No": 0}
    dfBinaryY = dfBinaryY.applymap(lambda v: 0 if v == "No" else (1 if v == "Yes" else v))
    
    
    # binaryY = dfBinaryY
    binaryY = np.array(dfBinaryY)
    # fill out missing values (TODO: MAY NEED TO DELETE THIS LATER)
    # binaryY = np.where(pd.isnull(binaryY), 0, binaryY)
    categoricalY = np.array(dfCategoricY)
    numericY = np.array(dfNumericY)
    


    return X, binaryY, categoricalY, numericY, binaryLabelDict, categoricLabelDict, numericLabelDict, allFeatureList 



    